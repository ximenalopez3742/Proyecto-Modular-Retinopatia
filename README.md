# Aplicación de Métodos de Dinámica No Lineal (Dimensión Fractal) para Clasificar Imágenes de Retinopatía Diabética

**Proyecto Modular — Ing. en Electrónica y Computación**  
Centro Universitario de los Lagos, Universidad de Guadalajara  
Autora: Ximena López Guerra | Asesor: Dr. Rider Jaimes Reátegui  
Junio de 2026

---

## Descripción general

Este repositorio contiene los algoritmos desarrollados en MATLAB para la caracterización cuantitativa de la arquitectura vascular retiniana mediante métodos de dinámica no lineal. A partir de imágenes de fondo de ojo, el pipeline extrae la vasculatura, calcula la **dimensión fractal (D_f)** mediante el método de Box-Counting y la **lacunaridad media (Λ̄)** mediante el método de deslizamiento de cajas. Ambas métricas permiten discriminar entre retinas sanas y los diferentes estadios de la retinopatía diabética (RDNP, RDP-E1, RDP-E2, RDP-E3).

---

## Estructura del repositorio

```
/
├── Vasos_Retina.m       # Preprocesamiento y segmentación vascular (Filtro de Frangi)
├── Df_Retina.m          # Cálculo de la dimensión fractal (Box-Counting)
├── Lacunaridad.m        # Cálculo de la lacunaridad media
└── README.md
```

---

## Requisitos

- **MATLAB** R2021a o superior
- Toolboxes requeridos:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox (para `tinv` en `Df_Retina.m`)

---

## Pipeline general

```
Imagen RGB de fondo de ojo
        │
        ▼
  Vasos_Retina.m
  ┌─────────────────────────────────────┐
  │  Canal verde + CLAHE adaptativo     │
  │  Máscara del campo visual           │
  │  Detección y neutralización         │
  │  de exudados                        │
  │  Filtro de Frangi multiescala       │
  │  Umbralización adaptativa           │
  │  Refinamiento morfológico y         │
  │  filtro geométrico                  │
  └────────────────┬────────────────────┘
                   │  imagen binaria de vasos
          ┌────────┴────────┐
          ▼                 ▼
    Df_Retina.m       Lacunaridad.m
    Box-Counting      Deslizamiento
    → D_f, SE, r      → Λ por escala
                        → Λ̄ media
```

---

## Descripción de archivos

### `Vasos_Retina.m` — Segmentación de la arquitectura vascular

Procesa una imagen RGB de fondo de ojo y genera una imagen binaria con la vasculatura retiniana aislada. Esta imagen es la entrada para `Df_Retina.m` y `Lacunaridad.m`.

**Etapas del algoritmo:**

1. **Extracción del canal verde**  
   Se aísla el canal verde de la imagen RGB, que ofrece el mayor contraste natural entre los vasos sanguíneos y el fondo retiniano al presentar la menor absorción lumínica por los pigmentos de la retina.

2. **CLAHE adaptativo**  
   Se aplica ecualización adaptativa de histograma con limitación de contraste (`adapthisteq`). El parámetro `ClipLimit` se calcula dinámicamente a partir de la desviación estándar de cada imagen (`0.08 * std2(canal_verde)`), lo que evita amplificaciones inconsistentes entre estadios clínicos con distintas condiciones de captura.

3. **Generación de máscara del campo visual**  
   Se binariza la imagen en escala de grises con umbral de intensidad 15, se selecciona la componente conexa de mayor área y se refina mediante relleno de huecos (`imfill`) y erosión morfológica con un disco de radio proporcional a la escala de la imagen. Esto restringe todos los cálculos al área retiniana efectiva, excluyendo el fondo negro.

4. **Detección y neutralización de exudados**  
   Se aplica un filtro tophat sobre la diferencia entre los canales rojo y verde para identificar exudados brillantes mediante un umbral estadístico (`media + 2.5 * desviación estándar`). Las regiones detectadas se neutralizan antes del filtrado de Frangi para evitar falsos positivos en la segmentación vascular.

5. **Filtro de Frangi multiescala (Hessian-based)**  
   Se implementa el filtro de Frangi (Frangi et al., 1998) calculando el tensor Hessiano normalizado a múltiples escalas gaussianas (`sigma_min` a `sigma_max`, proporcionales a la resolución de la imagen respecto a una referencia de 2,000 px). Para cada escala se evalúa la función de *vesselness* 2D:

   ```
   V(σ) = exp(−R_B² / 2β²) · (1 − exp(−S² / 2c²)),   λ₂ < 0
   ```

   donde `R_B = λ₁/λ₂` es la medida de blobness y `S = sqrt(λ₁² + λ₂²)` la magnitud del Hessiano. El mapa de vasos final es el máximo de `V(σ)` sobre todas las escalas. Este método detecta simultáneamente vasos de diferente calibre sin requerir orientaciones predefinidas, superando las limitaciones de detectores de bordes clásicos (como Canny) para vasos finos.

6. **Umbralización adaptativa**  
   Se combina el umbral de Otsu con el percentil 92 de los píxeles dentro del campo visual, tomando el mínimo de ambos (con un mínimo absoluto de 0.05), para binarizar el mapa de *vesselness* de forma robusta ante variaciones de contraste entre estadios.

7. **Refinamiento morfológico y filtro geométrico**  
   Se aplican apertura y cierre morfológico con elementos estructurantes proporcionales a la escala. Se calcula el esqueleto de cada componente conexa y se descartan las que no cumplan criterios de longitud mínima de esqueleto, excentricidad mínima y razón área/longitud, conservando únicamente estructuras consistentes con vasos sanguíneos. Finalmente, se aplica cierre con elementos lineales en cuatro orientaciones (0°, 45°, 90°, 135°) para reconectar fragmentos vasculares discontinuos.

**Parámetros clave:**

| Parámetro | Descripción | Valor |
|---|---|---|
| `ref` | Resolución de referencia (px) | 2000 |
| `sigma_min / sigma_max` | Rango de escalas Frangi (proporcional a `escala`) | `0.5·escala` / `4.0·escala` |
| `b_frangi` | Parámetro de blobness β | 0.5 |
| `g_frangi` | Parámetro de magnitud c | 15 |
| `clahe_cliplimit` | Límite de contraste CLAHE | `0.08 · std2(canal_verde)` |
| `umbral_pct` | Percentil para umbralización | 92 |
| `k_ex` | Factor σ para detección de exudados | 2.5 |

**Salidas:**
- `imagen_boxcounting`: imagen binaria lógica con la arquitectura vascular segmentada (entrada para los otros dos scripts).
- Figura con 6 paneles: imagen original, canal verde, canal verde + CLAHE, mapa Frangi, vasos binarizados y superposición sobre la imagen RGB.

**Configuración — ruta de entrada:**
```matlab
ruta_imagen = 'D:\Dataset_Retina\...\imagen.png';
```

---

### `Df_Retina.m` — Dimensión Fractal por Box-Counting

Calcula la dimensión fractal `D_f` de la imagen binaria de vasos generada por `Vasos_Retina.m`, cuantificando cómo la estructura vascular ocupa el espacio bidimensional a través de diferentes escalas.

**Etapas del algoritmo:**

1. **Acondicionamiento (padding)**  
   La imagen binaria se somete a relleno simétrico para que sus dimensiones sean múltiplos exactos del número de cajas `k` por eje, evitando errores de truncamiento en los bordes de la cuadrícula.

2. **Conteo de cajas**  
   Se evalúan 9 escalas con `k ∈ {2, 4, 8, 16, 32, 64, 128, 256, 512}` (longitud de caja `l = 1/k`). Para cada escala, el algoritmo recorre sistemáticamente la cuadrícula y contabiliza el número de cajas ocupadas `N(l)`, definiendo una caja como ocupada si contiene al menos un píxel de vaso. El intervalo de escalas válido corresponde al rango donde la relación log-log es consistentemente lineal: por debajo de `k = 2` la cuadrícula es demasiado gruesa; por encima de `k = 512` el tamaño de caja se aproxima al ancho de los vasos más finos, introduciendo sesgo dependiente de la resolución.

3. **Ajuste lineal por mínimos cuadrados**  
   Se grafican los puntos `(−ln(l), ln(N(l)))` y se ajusta una recta por mínimos cuadrados. La pendiente representa la dimensión fractal `D_f`, con base en la definición de Mandelbrot:

   ```
   D_f = lim_{l→0}  ln(N(l)) / −ln(l)
   ```

4. **Validación estadística**  
   Se calculan:
   - **Coeficiente de correlación de Pearson (r):** un valor cercano a 1 confirma comportamiento fractal consistente.
   - **Error estándar de la pendiente (SE):** calculado a partir de la varianza residual del ajuste (`MSE / S_xx`), junto con el intervalo de confianza al 95% (t de Student, `gl = n − 2 = 7`).

**Resultados en consola:**

```
==============================================
   RESULTADOS DEL AJUSTE BOX-COUNTING
==============================================
Dimensión fractal  Df  = 1.4514
Intercepto         a0  = X.XXXX
Correlación de Pearson  r   = 0.999214
----------------------------------------------
Error estándar de Df   SE  = 0.0218
==============================================
```

**Salidas:**
- Tabla `(l, N(l), −ln(l), ln(N(l)))` impresa en consola.
- `D_f`, `a0`, `r`, `SE` y el intervalo de confianza al 95%.
- Figura log-log con datos experimentales y recta de mejor ajuste.
- Figuras de la cuadrícula superpuesta sobre la imagen para cada una de las 9 escalas.

**Configuración — ruta de entrada:**
```matlab
ruta_vasos = 'ruta\a\imagen_vasos.png';
```

---

### `Lacunaridad.m` — Lacunaridad media por deslizamiento de cajas

Calcula la lacunaridad media `Λ̄` de la imagen binaria de vasos como métrica complementaria a `D_f`. Mientras que `D_f` describe la complejidad global de la ramificación, la lacunaridad captura la **heterogeneidad espacial** de la red: dos estructuras con el mismo `D_f` pueden diferir significativamente en lacunaridad si sus vacíos están distribuidos de forma distinta.

**Etapas del algoritmo:**

1. **Acondicionamiento (padding)**  
   Idéntico al empleado en `Df_Retina.m`: relleno simétrico para múltiplos exactos de `k`.

2. **Deslizamiento de cajas y cálculo de masas**  
   Para cada escala `k ∈ {2, 4, 8, 16, 32, 64, 128, 256, 512}`, se divide la imagen en `k × k` cajas y se calcula la **masa** de cada caja (número de píxeles de vaso que contiene). A diferencia del Box-Counting —donde solo se registra si una caja está ocupada o vacía— aquí se cuantifica cuántos píxeles de vaso contiene cada caja, capturando la distribución espacial de la densidad vascular.

3. **Cálculo de la lacunaridad por escala**  
   A partir de la distribución de masas `{mᵢ}`:

   ```
   Λ(k) = σ²(k) / μ²(k) + 1
   ```

   donde `μ` y `σ²` son la media y varianza de las masas para la escala `k`. Un valor `Λ ≈ 1` indica distribución espacial homogénea; valores mayores reflejan heterogeneidad creciente (zonas de alta densidad vascular alternadas con zonas avasculares extensas).

4. **Lacunaridad media**  
   Se reporta `Λ̄ = (1/9) · Σ Λ(k)` como descriptor escalar por imagen, que integra la información de heterogeneidad en todos los niveles de escala analizados.

**Resultados en consola:**

```
Tabla de Lacunaridad
 k (cajas)       N(l)         Λ(r)
-------------------------------------
        2          4        X.XXXX
        4         14        X.XXXX
       ...
       512      13993        X.XXXX
-------------------------------------
Lacunaridad media = 8.6121
```

**Salidas:**
- Tabla `(k, N(l), Λ(k))` impresa en consola.
- `Λ̄` media impresa en consola.
- Figura de la curva de lacunaridad `Λ` en función de `log(1/l)`.

**Configuración — ruta de entrada:**
```matlab
ruta_vasos = 'ruta\a\imagen_vasos.png';
```

---

## Uso paso a paso

### 1. Segmentación vascular
Abrir `Vasos_Retina.m`, configurar la ruta de la imagen de fondo de ojo y ejecutar:
```matlab
ruta_imagen = 'ruta\a\imagen_original.png';
```
El script genera la variable `imagen_boxcounting` (imagen binaria de vasos) y la guarda o usa directamente en el workspace de MATLAB.

### 2. Guardar la imagen de vasos
Para usar la imagen binarizada como entrada en los siguientes scripts:
```matlab
imwrite(imagen_boxcounting, 'ruta\a\imagen_vasos.png');
```

### 3. Cálculo de la dimensión fractal
Abrir `Df_Retina.m`, configurar la ruta de la imagen de vasos y ejecutar:
```matlab
ruta_vasos = 'ruta\a\imagen_vasos.png';
```

### 4. Cálculo de la lacunaridad
Abrir `Lacunaridad.m`, configurar la misma ruta y ejecutar:
```matlab
ruta_vasos = 'ruta\a\imagen_vasos.png';
```

---

## Interpretación de resultados

| Métrica | Descripción | Valor típico (retina sana) |
|---|---|---|
| `D_f` | Complejidad global de la ramificación vascular | ~1.35 |
| `SE` | Error estándar de la pendiente del ajuste | < 0.02 (buen ajuste) |
| `r` | Coeficiente de correlación de Pearson | ≥ 0.99 |
| `Λ̄` | Heterogeneidad espacial de la red vascular | ~15.9 |

**Comportamiento esperado por estadio clínico:**

| Estadio | `D_f` (media) | `Λ̄` (media) | Interpretación |
|---|---|---|---|
| Sana | 1.3457 | 15.8861 | Red vascular regular y moderadamente compleja |
| RDNP | 1.4147 | 10.2654 | Mayor complejidad; red más densa y homogénea |
| RDP-E1 | 1.4119 | 10.5798 | Neovascularización incipiente |
| RDP-E2 | 1.4167 | 10.3460 | Neovascularización establecida |
| RDP-E3 | 1.3399 | 16.3185 | Pérdida de vasos finos; alta irregularidad espacial |

---

## Referencias

- Frangi, A. F., Niessen, W. J., Vincken, K. L., & Viergever, M. A. (1998). Multiscale vessel enhancement filtering. *MICCAI'98*, 130–137.
- Lynch, S. (2025). *Dynamical Systems with Applications using MATLAB®* (3rd ed.). Birkhäuser, Springer Nature.
- Masters, B. R. (2004). Fractal analysis of the vascular tree in the human retina. *Annual Review of Biomedical Engineering*, 6, 427–452.
- Liew, G., et al. (2008). The retinal vasculature as a fractal. *Ophthalmology*, 115(11), 1951–1956.
- Stosic, T., & Stosic, B. (2005). Multifractal analysis of human retinal vessels.
- Tolle, C. R., et al. (2003). Lacunarity definition for ramified data sets based on optimal cover. *Physica D*, 179(3–4), 129–152.
- The MathWorks, Inc. (2025). *adapthisteq — Contrast-limited adaptive histogram equalization (CLAHE)*. MathWorks.
