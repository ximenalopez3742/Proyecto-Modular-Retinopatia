% Extracción de vasos de la retina humana
clc; 
clear; 
close all;

ruta_imagen='35.jpg';
imagen_rgb=imread(ruta_imagen);

[alto,ancho,canales]=size(imagen_rgb);
fprintf('Resolución Imagen: %d x %d px\n',ancho,alto);

ref=2000; % 2000 px de referencia a la imagen
img=max(alto,ancho);
escala=img/ref; % Factor de escala a la imagen de referencia
escala2=escala^2;

% Frangi escalas
sigma_min=max(0.5,0.5*escala); % Vasos más delgados
sigma_max=max(2.0,4.0*escala); % Vasos más gruesos
sigma_paso=max(0.5,0.5*escala);
frangi_escalas=sigma_min:sigma_paso:sigma_max;

% Parámetros
area_min_px=max(10,round(40*escala2));
area_min_geo=max(10,round(60*escala2));
longitud_min=max(3,round(8*escala));
exc_min=0.45;
area_grande=max(50,round(300*escala2));
ratio_max=8.0;

% Canal Verde + CLAHE
canal_verde=im2double(imagen_rgb(:,:,2));
clahe_cliplimit=max(0.005,min(0.025,0.08*std2(canal_verde)));
canal_clahe=adapthisteq(canal_verde,'NumTiles',[8 8],'ClipLimit',clahe_cliplimit,'Distribution','rayleigh');

% Máscara
gris=rgb2gray(imagen_rgb);
binaria=gris>15; % Covertir a gris casi negro
comp_conect=bwlabel(binaria);
area_comp=regionprops(comp_conect,'Area');
[~,indice]=max([area_comp.Area]);
mascara=imfill(comp_conect==indice,'holes');
radio_erosion=max(5,round(20*escala));
mascara=imerode(mascara,strel('disk',radio_erosion));

% Detección de exudados duros
R_ex=im2double(imagen_rgb(:,:,1));
G_ex=im2double(imagen_rgb(:,:,2));
B_ex=im2double(imagen_rgb(:,:,3));

radio_fondo=max(15,round(40*escala));
h_fondo=fspecial('disk',radio_fondo);
canal_det=G_ex-B_ex;
fondo_local=imfilter(canal_det,h_fondo,'replicate');
tophat_ex=canal_det-fondo_local;

med_ex=mean(tophat_ex(mascara));
desv_ex=std(tophat_ex(mascara));
k_ex=2.5;
umbral_ex=med_ex+k_ex*desv_ex;

candidatos_ex=(tophat_ex>umbral_ex) & mascara;
masc_azul=B_ex<0.55;
candidatos_ex=candidatos_ex & masc_azul;

radio_cierre_ex=max(2,round(5*escala));
area_min_ex=max(20,round(120*escala2));
area_max_ex=round(0.08*sum(mascara(:)));

candidatos_ex=imclose(candidatos_ex,strel('disk',radio_cierre_ex));
candidatos_ex=bwareaopen(candidatos_ex,area_min_ex);

comp_ex=bwlabel(candidatos_ex);
props_ex=regionprops(comp_ex,'Area','PixelIdxList');

for k_comp=1:length(props_ex)
    if props_ex(k_comp).Area > area_max_ex
        candidatos_ex(props_ex(k_comp).PixelIdxList) = false;
    end
end

radio_dilat_ex=max(2,round(3*escala));
mascara_exudados=imdilate(candidatos_ex,strel('disk',radio_dilat_ex));
mascara_exudados=mascara_exudados & mascara;

n_px_ex=sum(mascara_exudados(:));
n_ex=sum(mascara(:));

% Neutralizar zona de exudados
canal_frangi=imcomplement(canal_clahe);
valor_relleno=mean2(canal_frangi(mascara & ~mascara_exudados));
canal_frangi(mascara_exudados)=valor_relleno;

% Frangi multiescala
b_frangi=0.5;
g_frangi=15;
mapa_vasos=zeros(size(canal_frangi));

for sigma_f = frangi_escalas
    img_s=imgaussfilt(canal_frangi,sigma_f);
    % Hessiano 2D
    radio_h=round(3*sigma_f);
    [x_h,y_h]=meshgrid(-radio_h:radio_h,-radio_h:radio_h);
    g_h=exp(-(x_h.^2 + y_h.^2) / (2*sigma_f^2));
	Gxx_h=g_h.*(x_h.^2 / sigma_f^4 - 1 / sigma_f^2);
	Gxy_h=g_h.*(x_h.*y_h / sigma_f^4);
	Gyy_h=g_h.*(y_h.^2 / sigma_f^4-1 / sigma_f^2);
	e_h=sigma_f^2;
	Lxx=e_h*imfilter(img_s,Gxx_h,'replicate');
	Lxy=e_h*imfilter(img_s,Gxy_h,'replicate');
	Lyy=e_h*imfilter(img_s,Gyy_h,'replicate');
    % Valores propios
    D_h=sqrt((Lxx-Lyy).^2+ 4*Lxy.^2);
	tr_h=Lxx + Lyy;
	l1_h=(tr_h-D_h)/2;
	l2_h=(tr_h+D_h)/2;
    % Ordenar por valor absoluto
    sw_h=abs(l1_h) > abs(l2_h);
	l1_ord=l1_h;
	l2_ord=l2_h;
	l1_ord(sw_h)=l2_h(sw_h);
	l2_ord(sw_h)=l1_h(sw_h);

	Rb_f=l1_ord ./ (l2_ord+eps);
	S_f=sqrt(l1_ord.^2 + l2_ord.^2);
	v_f=zeros(size(canal_frangi));
	m_f=l2_ord < 0;
	v_f(m_f)=exp(-Rb_f(m_f).^2 / (2 * b_frangi^2)) .* (1-exp(-S_f(m_f).^2 / (2*g_frangi^2)));
	mapa_vasos=max(mapa_vasos, v_f);
end

mx_vasos=max(mapa_vasos(:));
if mx_vasos > 0
    mapa_vasos=mapa_vasos/mx_vasos;
end

% Umbralización
px_fov=mapa_vasos(mascara & ~mascara_exudados);
umbral_otsu=graythresh(px_fov);
umbral_pct=prctile(px_fov, 92);
umbral=max(min(umbral_otsu, umbral_pct), 0.05);
binaria=(mapa_vasos > umbral) & mascara & ~mascara_exudados;

% Limpieza
radio_open=max(1,round(0.5*escala));
radio_close=max(1,round(1.5*escala));
binaria=imopen(binaria,strel('disk',radio_open));
binaria=bwareaopen(binaria,area_min_px);
binaria=imclose(binaria,strel('disk',radio_close));
binaria=bwareaopen(binaria, area_min_px);
binaria=binaria & mascara;

% Reconectar fragmentos
radio_reconect=max(2,round(2.5*escala));
bin_dilatada=imdilate(binaria,strel('disk',radio_reconect));
bin_dilatada=bin_dilatada & mascara & ~mascara_exudados;
bin_erode=imerode(bin_dilatada,strel('disk',radio_reconect));
binaria_pre=(bin_erode|binaria) & mascara & ~mascara_exudados;

long_cierre_pre=max(3,round(5*escala));
binaria_pre=imclose(binaria_pre,strel('line',long_cierre_pre,0));
binaria_pre=imclose(binaria_pre,strel('line',long_cierre_pre,90));
binaria_pre=imclose(binaria_pre, strel('line',long_cierre_pre,45));
binaria_pre=imclose(binaria_pre, strel('line',long_cierre_pre,135));
binaria_pre=binaria_pre & mascara & ~mascara_exudados;

% Filtro geométrico
esqueleto_geo=bwskel(logical(binaria_pre),'MinBranchLength',0);
etiquetas_geo=bwlabel(binaria_pre);
props_geo=regionprops(etiquetas_geo,'Area','Eccentricity');
n_comp_geo=max(etiquetas_geo(:));
limpia=false(size(binaria_pre));

for k_geo=1:n_comp_geo
	mask_geo=(etiquetas_geo == k_geo);
	area_geo=props_geo(k_geo).Area;
	exc_geo=props_geo(k_geo).Eccentricity;
	long_skel=sum(esqueleto_geo(mask_geo));
	ratio_geo=area_geo / max(long_skel,1);
	conservar=false;

    if area_geo >= area_grande
        if long_skel >= longitud_min && ratio_geo <= ratio_max*4;
            conservar=true;
        end
    elseif area_geo >= area_min_geo
        if long_skel >= longitud_min && exc_geo >= exc_min && ratio_geo <= ratio_max
            conservar=true;
        end
    end

    if conservar
        limpia(mask_geo)=true;
    end
end

% Cierre
long_cierre=max(4,round(8*escala));
limpia=imclose(limpia,strel('line',long_cierre,0));
limpia=imclose(limpia,strel('line',long_cierre,90));
limpia=imclose(limpia,strel('line',long_cierre,45));
limpia=imclose(limpia,strel('line',long_cierre,135));
limpia=bwareaopen(limpia,area_min_geo);
limpia=limpia & mascara & ~mascara_exudados;

imagen_boxcounting=limpia;

% Conservación de vasos
diff_img=zeros(alto,ancho,3,'uint8');
diff_img(:,:,1)=uint8(binaria & ~imagen_boxcounting)*255;
diff_img(:,:,2)=uint8(imagen_boxcounting)*255;
% Superposición de vasos
ov=imagen_rgb;
r=ov(:,:,1);
g=ov(:,:,2);
b=ov(:,:,3);
r(imagen_boxcounting)=255;
g(imagen_boxcounting)=0;
b(imagen_boxcounting)=0;
ov(:,:,1)=r;
ov(:,:,2)=g;
ov(:,:,3)=b;

% Visualización de Resultados
figure();
subplot(2,3,1); imshow(imagen_rgb); title('Imagen Original');
subplot(2,3,2); imshow(canal_verde,[]); title('Canal Verde');
subplot(2,3,3); imshow(canal_clahe,[]); title('Canal Verde + CLAHE');
subplot(2,3,4); imshow(mapa_vasos,[]); title('Filtro Frangi');
subplot(2,3,5); imshow(imagen_boxcounting); title('Vasos Retinianos');
subplot(2,3,6); imshow(ov); title('Superposición de Vasos');
