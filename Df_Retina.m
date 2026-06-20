% Cálculo de la Dimensión Fractal – Método Box-Counting
clc;
clear all;
close all;

ruta_vasos = '35_vasos.png';
vasos = imread(ruta_vasos);

if ~islogical(vasos)
    vasos = vasos > 128;
end
if ndims(vasos) == 3
    vasos = vasos(:,:,1) > 0;
end

fprintf('Imagen cargada: %s\n', ruta_vasos);
fprintf('Resolución: %d x %d px\n', size(vasos,2), size(vasos,1));

[filas, columnas] = size(vasos);
vasos = logical(vasos);

num_cajas   = [2 4 8 16 32 64 128 256 512];
num_escalas = length(num_cajas);
N = zeros(1, num_escalas);
l = zeros(1, num_escalas);

for i = 1:num_escalas
    k        = num_cajas(i);
    l(i)     = 1/k;
    tam_caja = ceil(max(filas,columnas) / k);
    filas_pad    = k * tam_caja;
    columnas_pad = k * tam_caja;
    pad_filas    = filas_pad - filas;
    pad_columnas = columnas_pad - columnas;
    pad_sup = floor(pad_filas/2);
    pad_inf = ceil(pad_filas/2);
    pad_izq = floor(pad_columnas/2);
    pad_der = ceil(pad_columnas/2);
    imagen_pad = padarray(vasos, [pad_sup pad_izq], 0, 'pre');
    imagen_pad = padarray(imagen_pad, [pad_inf pad_der], 0, 'post');

    contador = 0;
    for r = 1:k
        for c = 1:k
            r_ini  = (r-1)*tam_caja + 1;
            r_fin  = r*tam_caja;
            c_ini  = (c-1)*tam_caja + 1;
            c_fin  = c*tam_caja;
            bloque = imagen_pad(r_ini:r_fin, c_ini:c_fin);
            if any(bloque(:))
                contador = contador + 1;
            end
        end
    end
    N(i) = contador;

    % Visualización de cuadrícula por escala
    figure;
    imshow(imagen_pad, []);
    hold on
    for r = 0:k
        y_lin = r * tam_caja + 0.5;
        line([1 columnas_pad], [y_lin y_lin], 'Color', 'y', 'LineWidth', 0.5);
    end
    for c = 0:k
        x_lin = c * tam_caja + 0.5;
        line([x_lin x_lin], [1 filas_pad], 'Color', 'y', 'LineWidth', 0.5);
    end
    title(['l = 1/', num2str(k), ',  N(l) = ', num2str(contador)])
    hold off
end

% Ajuste log-log
x = -log(l);
y =  log(N);

l_frac = strings(num_escalas, 1);
for i = 1:num_escalas
    l_frac(i) = "1/" + num2str(num_cajas(i));
end

fprintf('\nTabla de conteo de cajas (Box-Counting)\n');
fprintf('%6s %8s %10s %12s\n', 'l', 'N(l)', '-ln(l)', 'ln(N(l))');
fprintf('---------------------------------------------\n');
for i = 1:num_escalas
    fprintf('%6s %8d %10.4f %12.4f\n', l_frac{i}, N(i), x(i), y(i));
end
fprintf('---------------------------------------------\n\n');

% Ajuste por mínimos cuadrados
n  = num_escalas;
X  = x;
Y  = y;

q  = sum(X .* Y);
q1 = sum(X) * sum(Y);
z  = sum(X.^2);
z1 = (sum(X))^2;

a1   = (n*q - q1) / (n*z - z1); % pendiente = Df
a0   = mean(Y) - a1 * mean(X);  % Intercepto
Y_RL = a0 + a1 * X;             % Valores ajustados

% Error estándar de la pendiente (SE)
% Residuos del ajuste
residuos = Y - Y_RL;

MSE = sum(residuos.^2) / (n - 2);

% Error estándar de la pendiente
% SE(a1) = sqrt( MSE / sum((X - mean(X))^2) )
Sxx  = sum((X - mean(X)).^2);
SE_pendiente = sqrt(MSE / Sxx);

% Intervalo de confianza al 95% 
t_critico = tinv(0.975, n - 2);         
IC_inf = a1 - t_critico * SE_pendiente;
IC_sup = a1 + t_critico * SE_pendiente;

% Coeficiente de correlación de Pearson 
r_corr = (n*sum(X.*Y) - sum(X)*sum(Y)) / ...
         (sqrt(n*sum(X.^2) - (sum(X))^2) * ...
          sqrt(n*sum(Y.^2) - (sum(Y))^2));

% Resultados
fprintf('Dimensión fractal  Df  = %.4f\n', a1);
fprintf('Intercepto         a0  = %.4f\n', a0);
fprintf('Correlación de Pearson  r   = %.4f\n', r_corr);
fprintf('Error estándar de Df   SE  = %.4f\n', SE_pendiente);

% Gráfica log-log
figure;
plot(X, Y, 'bo', 'MarkerFaceColor', 'b')
hold on
plot(X, Y_RL, 'r', 'LineWidth', 1.5)


grid on
xlabel('-ln(l)')
ylabel('ln(N(l))')
title(['Método Box-Counting  |  Df = ', num2str(a1, '%.4f'), ...
       '  |  SE = ', num2str(SE_pendiente, '%.4f'), ...
       '  |  r = ', num2str(r_corr, '%.4f')])
legend('Datos', ...
       ['Pendiente: Df = ', num2str(a1, '%.4f')], ...
       'Location', 'NorthWest')
hold off