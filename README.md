# Aplicación de Métodos de Dinámica No Lineal (Dimensión Fractal) para Clasificar Imágenes de Retinopatía Diabética

**Proyecto Modular — Ing. en Electrónica y Computación**  
Centro Universitario de los Lagos, Universidad de Guadalajara  
Autora: Ximena López Guerra | Asesor: Dr. Rider Jaimes Reátegui  
Junio de 2026

---

## Descripción general

Este repositorio contiene los algoritmos desarrollados en MATLAB para la caracterización cuantitativa de la arquitectura vascular retiniana mediante métodos de dinámica no lineal. A partir de imágenes de fondo de ojo, el pipeline extrae la vasculatura, calcula la **dimensión fractal ($D_f$)** mediante el método de Box-Counting y la **lacunaridad media ($\bar{\Lambda}$)** mediante el método de deslizamiento de cajas. Ambas métricas permiten discriminar entre retinas sanas y los diferentes estadios de la retinopatía diabética (RDNP, RDP-E1, RDP-E2, RDP-E3).

---

## Pipeline general

```
Imagen RGB de fondo de ojo
        │
        ▼
  Vasos_Retina.m
  ┌───────────────────────────────────────────────┐
  │  Canal verde + CLAHE adaptativo               │
  │  Máscara del campo visual                     │
  │  Detección y neutralización de exudados       │
  │  Filtro de Frangi multiescala                 │
  │  Umbralización adaptativa                     │
  │  Refinamiento morfológico y filtro geométrico │
  └────────────────┬──────────────────────────────┘
                   │  Imagen binaria de vasos
          ┌────────┴────────┐
          ▼                 ▼
    Df_Retina.m       Lacunaridad.m
    Box-Counting      Deslizamiento
    → Df, SE, r      → Λ por escala
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
   Se aplica ecualización adaptativa de histograma con limitación de contraste (`adapthisteq`). El parámetro `ClipLimit` se calcula dinámicamente a partir de la desviación estándar de cada imagen, lo que evita amplificaciones inconsistentes entre estadios clínicos con distintas condiciones de captura.

3. **Generación de máscara del campo visual**  
   Se binariza la imagen en escala de grises con umbral de intensidad 15, se selecciona la componente conexa de mayor área y se refina mediante relleno de huecos (`imfill`) y erosión morfológica con un disco de radio proporcional a la escala de la imagen. Esto restringe todos los cálculos al área retiniana efectiva, excluyendo el fondo negro.

4. **Detección y neutralización de exudados**  
   Se aplica un filtro tophat sobre la diferencia entre los canales rojo y verde para identificar exudados brillantes mediante un umbral estadístico. Las regiones detectadas se neutralizan antes del filtrado de Frangi para evitar falsos positivos en la segmentación vascular.

5. **Filtro de Frangi multiescala**  
   Se implementa el filtro de Frangi calculando el tensor Hessiano normalizado a múltiples escalas gaussianas (`sigma_min` a `sigma_max`, proporcionales a la resolución de la imagen respecto a una referencia de 2,000 px). Para cada escala se evalúa la función de *vesselness* 2D. Este método detecta simultáneamente vasos de diferente calibre sin requerir orientaciones predefinidas, superando las limitaciones de detectores de bordes clásicos (como Canny) para vasos finos.

6. **Umbralización adaptativa**  
   Se combina el umbral de Otsu con el percentil 92 de los píxeles dentro del campo visual, tomando el mínimo de ambos (con un mínimo absoluto de 0.05), para binarizar el mapa de *vesselness* de forma robusta ante variaciones de contraste entre estadios.

7. **Refinamiento morfológico y filtro geométrico**  
   Se aplican apertura y cierre morfológico con elementos estructurantes proporcionales a la escala. Se calcula el esqueleto de cada componente conexa y se descartan las que no cumplan criterios de longitud mínima de esqueleto, excentricidad mínima y razón área/longitud, conservando únicamente estructuras consistentes con vasos sanguíneos. Finalmente, se aplica cierre con elementos lineales en cuatro orientaciones (0°, 45°, 90°, 135°) para reconectar fragmentos vasculares discontinuos.


**Salidas:**
- `imagen_boxcounting`: Imagen binaria lógica con la arquitectura vascular segmentada (entrada para los otros dos scripts).
- `Figura con 6 paneles`: Imagen original, canal verde, canal verde + CLAHE, filtro Frangi, vasos retinianos y superposición sobre la imagen RGB.

---

### `Df_Retina.m` — Dimensión Fractal por Box-Counting

Calcula la dimensión fractal $D_f$ de la imagen binaria de vasos generada por `Vasos_Retina.m`, cuantificando cómo la estructura vascular ocupa el espacio bidimensional a través de diferentes escalas.

**Etapas del algoritmo:**

1. **Acondicionamiento (padding)**  
   La imagen binaria se somete a relleno simétrico para que sus dimensiones sean múltiplos exactos del número de cajas `k` por eje, evitando errores de truncamiento en los bordes de la cuadrícula.

2. **Conteo de cajas**  
   Se evalúan 9 escalas con `k ∈ {2, 4, 8, 16, 32, 64, 128, 256, 512}` (longitud de caja `l = 1/k`). Para cada escala, el algoritmo recorre sistemáticamente la cuadrícula y contabiliza el número de cajas ocupadas `N(l)`, definiendo una caja como ocupada si contiene al menos un píxel de vaso. El intervalo de escalas válido corresponde al rango donde la relación log-log es consistentemente lineal: por debajo de `k = 2` la cuadrícula es demasiado gruesa; por encima de `k = 512` el tamaño de caja se aproxima al ancho de los vasos más finos, introduciendo sesgo dependiente de la resolución.

3. **Ajuste lineal por mínimos cuadrados**  
   Se grafican los puntos `(−ln(l), ln(N(l)))` y se ajusta una recta por mínimos cuadrados. La pendiente representa la dimensión fractal $D_f$.
   
4. **Validación estadística**  
   Se calculan:
   - **Coeficiente de correlación de Pearson (r):** un valor cercano a 1 confirma comportamiento fractal consistente.
   - **Error estándar de la pendiente (SE):** calculado a partir de la varianza residual del ajuste, junto con el intervalo de confianza al 95%.

**Salidas:**
- Tabla de $l$, $N(l)$, $−ln(l)$, $ln(N(l))$ impresa en consola.
- $D_f$, $a0$, $r$, $SE$.
- Figura $log-log$ con datos experimentales y recta de mejor ajuste.
- Figuras de la cuadrícula superpuesta sobre la imagen para cada una de las 9 escalas.

---

### `Lacunaridad.m` — Lacunaridad media por deslizamiento de cajas

Calcula la lacunaridad media $\bar{\Lambda}$ de la imagen binaria de vasos como métrica complementaria a $D_f$. Mientras que $D_f$ describe la complejidad global de la ramificación, la lacunaridad captura la **heterogeneidad espacial** de la red: dos estructuras con el mismo $D_f$ pueden diferir significativamente en lacunaridad si sus vacíos están distribuidos de forma distinta.

**Etapas del algoritmo:**

1. **Acondicionamiento (padding)**  
   Idéntico al empleado en `Df_Retina.m`: relleno simétrico para múltiplos exactos de `k`.

2. **Deslizamiento de cajas y cálculo de masas**  
   Para cada escala `k ∈ {2, 4, 8, 16, 32, 64, 128, 256, 512}`, se divide la imagen en `k × k` cajas y se calcula la **masa** de cada caja (número de píxeles de vaso que contiene). A diferencia del Box-Counting (donde solo se registra si una caja está ocupada o vacía) aquí se cuantifica cuántos píxeles de vaso contiene cada caja, capturando la distribución espacial de la densidad vascular.

3. **Cálculo de la lacunaridad por escala**  
   A partir de la distribución de masas `{mᵢ}`:

   ```
   Λ(k) = σ²(k) / μ²(k) + 1
   ```

   donde `μ` y `σ²` son la media y varianza de las masas para la escala `k`. Un valor `Λ ≈ 1` indica distribución espacial homogénea; valores mayores reflejan heterogeneidad creciente (zonas de alta densidad vascular alternadas con zonas avasculares extensas).

4. **Lacunaridad media**  
   Se reporta `Λ̄ = (1/9) · Σ Λ(k)` como descriptor escalar por imagen, que integra la información de heterogeneidad en todos los niveles de escala analizados.

**Salidas:**
- Tabla de $k$, $N(l)$, $Λ$ impresa en consola.
- $\bar{\Lambda}$ media impresa en consola.
- Figura de la curva de lacunaridad $\Lambda$ en función de $ln(1/l)$.
