# Guía de instalación y ejecución — Slip Rate v4

## ¿Qué es este kit?
Pequeño conjunto de archivos para estimar **tasas de deslizamiento** (pendiente lineal) a partir de eventos paleosísmicos con incertidumbre. 
Genera figuras por escenario y nivel de incertidumbre, mostrando la **recta** con la **tasa positiva hacia el presente** y un **eje superior** en **AP** (años antes de 2025).

**Archivos:**
- `events_sinteticos.csv` (plantilla de datos; edítalo para tus sitios/eventos).
- `sliprate_cli_v4.py` (script principal).

---

## Opción A — SIN instalar nada (Google Colab)
1. Abre https://colab.research.google.com/ e inicia sesión.
2. Crea un cuaderno nuevo (Archivo → Nuevo cuaderno).
3. Instala librerías en una celda:
   ```python
   !pip install numpy pandas matplotlib
   ```
4. Sube `events_sinteticos.csv` y `sliprate_cli_v4.py` (icono de carpeta → Subir).
5. Ejecuta (en una celda):
   ```python
   !python sliprate_cli_v4.py events_sinteticos.csv
   !python sliprate_cli_v4.py events_sinteticos.csv 2000
   !python sliprate_cli_v4.py events_sinteticos.csv 2000 3500
   !python sliprate_cli_v4.py events_sinteticos.csv 3000 3500 0.5,1.0,2.0
   ```
6. Descarga los PNG desde la carpeta `events_sinteticos_salidas_v4/` (navegador de archivos de Colab).

---

## Opción B — Instalación local (Windows/macOS/Linux)

### 1) Instalar Python 3.10+
- Windows/macOS: https://www.python.org/downloads/ (en Windows marca “Add Python to PATH”).
- Linux: suele venir instalado (`python3 --version`).

### 2) Instalar librerías necesarias
Con **pip**:
```bash
pip install numpy pandas matplotlib
# o (según tu sistema): pip3 install numpy pandas matplotlib
```
Con **conda** (opcional):
```bash
conda create -n sliprate python=3.10 numpy pandas matplotlib -y
conda activate sliprate
```

### 3) Ejecutar el script
Coloca `events_sinteticos.csv` y `sliprate_cli_v4.py` en la misma carpeta.
```bash
# Todo el registro
python3 sliprate_cli_v4.py events_sinteticos.csv

# Ventana de tasa 0–2000 BP
python3 sliprate_cli_v4.py events_sinteticos.csv 2000

# Registro truncado ≤ 3500 BP
python3 sliprate_cli_v4.py events_sinteticos.csv 2000 3500

# Incertidumbre (baja, media, alta)
python3 sliprate_cli_v4.py events_sinteticos.csv 3000 3500 0.5,1.0,2.0
```

---

## Parámetros clave
- **`rate_window_bp`** (2º argumento): Ventana temporal (0–X años BP) donde se estima la **pendiente**. Si se omite, usa todo el rango disponible.
- **`record_window_bp`** (3º): Trunca el registro a edades ≤ X BP (simula que el registro sedimentario **pierde eventos** más antiguos).
- **`sigma_scales`** (4º): Factores para escalar **todas** las σ (edad y desplazamiento). Por defecto: `0.5,1.0,2.0` (baja, media, alta).

---

## Lectura de resultados
- Carpeta de salida: `*_salidas_v4/`.
- **`timeline_eventos.png`**: rangos de edad (min–max) y edad preferida por escenario.
- **`<Escenario>_SR_<baja|media|alta>.png`**: 
  - Banda 5–95% del **acumulado** y curva mediana.
  - **Recta** con etiqueta: `pendiente = … mm/año; Σdisp = … m; N = …` (eventos usados en la ventana efectiva).
- **Ejes**:
  - Abajo: **BP** (0 BP = **1950**, más joven a la derecha).
  - Arriba: **AP** (Años antes de **2025**), con conversión **AP = BP + 75** (0 AP = 2025).

---

## Editar el CSV
Columnas por evento:
- `age_cal_yrBP`, `age_sigma_yr` (años BP y su 1σ).
- `displacement_m`, `disp_sigma_m` (m y su 1σ).
- `scenario`, `event_id` (agrupación).

Sugerencia: comienza con la plantilla y **sustituye** edades/incertidumbres por las del paper. 

---

## Problemas comunes
- **python no se reconoce** → prueba `python3` o `py -3` (Windows).
- **ModuleNotFoundError** → falta `numpy/pandas/matplotlib` (ver instalación).
- **FileNotFoundError** → verifica estar en la carpeta correcta o usa **ruta absoluta**.
- **Gráfica “inversa”** → recuerda: el acumulado sube al **presente**; con BP en X, la recta tiene pendiente **numérica negativa**, pero la **tasa física** mostrada es **positiva hacia el presente**.

---

## (cheat‑sheet)
```bash
python3 sliprate_cli_v4.py events_sinteticos.csv                 # todo
python3 sliprate_cli_v4.py events_sinteticos.csv 2000            # ventana de tasa
python3 sliprate_cli_v4.py events_sinteticos.csv 3000 3500       # + registro truncado
python3 sliprate_cli_v4.py events_sinteticos.csv 3000 3500 0.5,1.0,2.0  # + incertidumbre
```
