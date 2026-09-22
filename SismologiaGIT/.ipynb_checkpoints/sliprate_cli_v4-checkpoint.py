# -*- coding: utf-8 -*-
"""
Script: sliprate_cli_v4.py (BP/AP + tasas positivas + Σdisp & N)

Calcula tasas de deslizamiento "normales" (pendiente lineal) a partir de eventos con incertidumbre.
Sin tasas móviles. Eje X = años calibrados BP (0 BP = 1950; joven a la derecha). El acumulado (Y) crece hacia el presente.

USO
----
python sliprate_cli_v4.py <csv_file> [rate_window_bp] [record_window_bp] [sigma_scales]

- <csv_file>: CSV con columnas:
  scenario,event_id,age_min,age_max,age_pref,age_cal_yrBP,age_sigma_yr,displacement_m,disp_sigma_m
- rate_window_bp (opcional): ventana (en años BP) para estimar la pendiente en 0..rate_window_bp.
  Si se omite, se usa todo el rango temporal disponible del escenario (tras truncado si aplica).
- record_window_bp (opcional): recorta eventos a edades <= record_window_bp (simula registro incompleto).
- sigma_scales (opcional): factores separados por comas para escalar las sigma de edad y desplazamiento.
  Por defecto "0.5,1.0,2.0" → baja, media, alta.

SALIDAS
-------
- timeline_eventos.png
- Por escenario y por incertidumbre: <escenario>_SR_<baja|media|alta>.png
  (banda 5–95%, curva mediana y recta con la tasa **positiva hacia el presente**; label incluye Σdisp y N eventos usados)
"""
import sys, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PRESENT_BP_YEAR = 1950
CURRENT_YEAR = 2025
AP_SHIFT = CURRENT_YEAR - PRESENT_BP_YEAR  # 75 años: AP = BP + 75

def bp_to_ap(x):
    return x + AP_SHIFT

def ap_to_bp(x):
    return x - AP_SHIFT

def cumulative_bp_increasing_to_present(ages_bp, disps):
    """y(BP) = suma de desplazamientos de eventos con edad >= BP (sube hacia el presente)."""
    order = np.argsort(ages_bp)  # ascendente (joven->viejo)
    a = ages_bp[order]
    d = disps[order]
    s = np.zeros_like(d)
    run = 0.0
    for k in range(len(d)-1, -1, -1):
        run += d[k]
        s[k] = run
    return a, s

def simulate_curves(df_s, n_iter, sigma_scale=1.0, record_window_bp=None, rng=None):
    """Simula curvas de acumulado en un grid uniforme de BP. Retorna percentiles y eventos para marcadores."""
    if record_window_bp is not None:
        df_s = df_s[df_s["age_cal_yrBP"] <= float(record_window_bp)].copy()
    if df_s.empty or len(df_s) < 2:
        raise ValueError("Muy pocos eventos tras 'record_window_bp'; se requieren >=2.")
    rng = np.random.default_rng(42) if rng is None else rng

    ages = df_s["age_cal_yrBP"].to_numpy(float)
    sig_a = df_s["age_sigma_yr"].to_numpy(float) * float(sigma_scale)
    disp = df_s["displacement_m"].to_numpy(float)
    sig_d = df_s["disp_sigma_m"].to_numpy(float) * float(sigma_scale)

    ages_samp = np.clip(rng.normal(ages, sig_a, size=(n_iter, len(ages))), 1, None)
    disp_samp = np.clip(rng.normal(disp, sig_d, size=(n_iter, len(disp))), 0, None)

    t_bp_max = float(np.max(ages_samp))
    t_grid = np.linspace(0.0, t_bp_max, 400)
    curves = np.empty((n_iter, t_grid.size))

    for i in range(n_iter):
        a, s = cumulative_bp_increasing_to_present(ages_samp[i], disp_samp[i])
        idx = np.searchsorted(a, t_grid, side="left")
        s_ext = np.concatenate([s, [0.0]])
        curves[i] = s_ext[idx]

    p05 = np.percentile(curves, 5, axis=0)
    p50 = np.percentile(curves, 50, axis=0)
    p95 = np.percentile(curves, 95, axis=0)

    ages_pref = df_s["age_pref"].to_numpy(float) if "age_pref" in df_s.columns else ages
    order = np.argsort(ages_pref)
    a_pref = ages_pref[order]
    d_mean = df_s["displacement_m"].to_numpy(float)[order]
    s = np.zeros_like(d_mean)
    run = 0.0
    for k in range(len(d_mean)-1, -1, -1):
        run += d_mean[k]
        s[k] = run

    return dict(t_grid=t_grid, p05=p05, p50=p50, p95=p95, events_bp=a_pref, events_cum=s, curves=curves)

def slope_distribution_from_curves(t_grid, curves, rate_window_bp=None):
    """Distribución de pendientes por OLS para cada realización de curva (m/año).
    NOTA: con x=BP creciente, la pendiente numérica sale NEGATIVA; la tasa física hacia el presente es -pendiente."""
    ok = np.isfinite(curves).all(axis=1)
    curves = curves[ok]
    if curves.shape[0] == 0:
        return np.array([np.nan])
    if rate_window_bp is None:
        mask = np.ones_like(t_grid, dtype=bool)
    else:
        mask = (t_grid >= 0.0) & (t_grid <= float(rate_window_bp))
        if mask.sum() < 3:
            mask = np.ones_like(t_grid, dtype=bool)
    x = t_grid[mask]
    x_center = x - x.mean()
    varx = np.sum(x_center**2)
    if varx == 0:
        return np.array([np.nan])
    slopes = np.empty(curves.shape[0])
    for i in range(curves.shape[0]):
        y = curves[i, mask]
        y_center = y - y.mean()
        slopes[i] = np.dot(x_center, y_center)/varx  # m/año (NEGATIVA en BP)
    return slopes

def event_summary_for_window(sub_df, rate_window_bp, record_window_bp):
    """Cuenta eventos y suma desplazamientos dentro de la ventana efectiva (0..min(rate_window_bp, record_window_bp)).
    Si ambas ventanas son None, usa todos los eventos."""
    import numpy as np
    if rate_window_bp is None and record_window_bp is None:
        mask = np.ones(len(sub_df), dtype=bool)
    else:
        cutoff = float("inf") if rate_window_bp is None else float(rate_window_bp)
        if record_window_bp is not None:
            cutoff = min(cutoff, float(record_window_bp))
        ages_pref = sub_df["age_pref"].to_numpy(float) if "age_pref" in sub_df.columns else sub_df["age_cal_yrBP"].to_numpy(float)
        mask = ages_pref <= cutoff
    n_events = int(mask.sum())
    sum_disp = float(sub_df.loc[mask, "displacement_m"].sum()) if n_events > 0 else 0.0
    return n_events, sum_disp

def plot_timeline(df, scenarios, outdir):
    fig, ax = plt.subplots(figsize=(9,5))
    y_positions = {sc:i for i, sc in enumerate(reversed(scenarios))}
    for sc in scenarios:
        sub = df[df["scenario"]==sc]
        y0 = y_positions[sc]
        for _, row in sub.iterrows():
            xmin, xmax = float(row["age_min"]), float(row["age_max"])
            ax.hlines(y0, xmin, xmax, linewidth=3)
            ax.plot(float(row.get("age_pref", row["age_cal_yrBP"])), y0, marker="o")
    ax.set_yticks(list(range(len(scenarios)))[::-1])
    ax.set_yticklabels(list(scenarios))
    ax.set_xlabel("Edad (años calibrados BP) — 0 BP = 1950 (presente radiocarbónico)")
    ax.set_title("Línea de tiempo de eventos por escenario")
    ax.invert_xaxis()
    # Eje superior en AP (años antes de 2025)
    ax_top = ax.secondary_xaxis("top", functions=(bp_to_ap, ap_to_bp))
    ax_top.set_xlabel("Años antes de 2025 (AP) — 0 AP = 2025")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "timeline_eventos.png"), dpi=150)
    plt.close()

def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    csv_file = sys.argv[1]
    rate_window_bp = None if len(sys.argv) < 3 else float(sys.argv[2])
    record_window_bp = None if len(sys.argv) < 4 else float(sys.argv[3])
    sigma_scales = "0.5,1.0,2.0" if len(sys.argv) < 5 else sys.argv[4]
    scales = [float(s) for s in sigma_scales.split(",")]
    labels = ["baja" if s<=0.75 else "media" if s<=1.25 else "alta" for s in scales]

    df = pd.read_csv(csv_file).sort_values(["scenario","age_cal_yrBP","event_id"]).reset_index(drop=True)
    scenarios = sorted(df["scenario"].unique())
    outdir = os.path.splitext(csv_file)[0] + "_salidas_v4"
    os.makedirs(outdir, exist_ok=True)

    plot_timeline(df, scenarios, outdir)

    rng = np.random.default_rng(12345)
    for sc in scenarios:
        sub = df[df["scenario"]==sc].copy()
        for sscale, lab in zip(scales, labels):
            out = simulate_curves(sub, n_iter=4000, sigma_scale=sscale, record_window_bp=record_window_bp, rng=rng)
            slopes = slope_distribution_from_curves(out["t_grid"], out["curves"], rate_window_bp=rate_window_bp)

            # Tasas POSITIVAS hacia el presente
            slopes_pos = -slopes
            sr_p05 = np.nanpercentile(slopes_pos, 5)*1000.0
            sr_p50 = np.nanpercentile(slopes_pos, 50)*1000.0
            sr_p95 = np.nanpercentile(slopes_pos, 95)*1000.0

            # Recta representativa: pendiente positiva para la etiqueta, pero negativa en BP para que suba a la derecha
            mask = (out["t_grid"] >= 0.0) & (out["t_grid"] <= (rate_window_bp if rate_window_bp is not None else out["t_grid"].max()))
            x = out["t_grid"][mask]
            y = out["p50"][mask]
            m_pos = np.nanpercentile(slopes_pos, 50)/1000.0  # m/año positivo hacia el presente
            m_line = -m_pos                                   # en eje BP, línea con pendiente negativa
            b_med = np.nanmean(y) - m_line * np.nanmean(x)
            xline = np.array([x.min(), x.max()])
            yline = b_med + m_line * xline

            # Resumen de eventos usados y suma de desplazamientos
            n_used, sum_disp = event_summary_for_window(sub, rate_window_bp, record_window_bp)

            plt.figure(figsize=(8,5))
            plt.fill_between(out["t_grid"], out["p05"], out["p95"], alpha=0.3, label="Acumulado 5–95%")
            plt.plot(out["t_grid"], out["p50"], label="Acumulado mediano")
            plt.plot(xline, yline, linestyle="--", label=f"Ventana de tiempo (pendiente = {sr_p50:.3f} mm/año; Σdisp={sum_disp:.2f} m; N={n_used})")
            ax = plt.gca()
            ax.invert_xaxis()
            plt.xlabel("Edad (años calibrados BP) — 0 BP = 1950")
            # Eje superior en AP
            ax_top = ax.secondary_xaxis("top", functions=(bp_to_ap, ap_to_bp))
            ax_top.set_xlabel("Años antes de 2025 (AP) — 0 AP = 2025")
            plt.ylabel("Desplazamiento acumulado (m)")
            title = f"{sc} — incertidumbre {lab} (σ×{sscale:g})"
            if record_window_bp is not None:
                title += f" | registro ≤ {int(record_window_bp)} aBP"
            if rate_window_bp is not None:
                title += f" | tasa en 0–{int(rate_window_bp)} aBP"
            plt.title(title)
            txt = f"Tasa (mm/año, hacia el presente): {sr_p50:.3f}  (p05–p95: {sr_p05:.3f}–{sr_p95:.3f})"
            plt.text(0.02, 0.02, txt, transform=ax.transAxes, ha="left", va="bottom",
                     bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))
            plt.legend(loc="upper left")
            plt.tight_layout()
            fname = os.path.join(outdir, f"{sc}_SR_{lab}.png")
            plt.savefig(fname, dpi=150)
            plt.close()

    print("Listo. Salidas en:", outdir)

if __name__ == "__main__":
    main()
