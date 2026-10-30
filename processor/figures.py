"""
The report's own figures (matplotlib, PNG): the few that combine steps or summarise across
samples. Everything else is the upstream apps' own figures, copied as they are. One look
across them: a light chart surface, thin marks, recessive axes, a fixed categorical order
(never cycled; beyond eight series the rest fold into "Other"), a single-hue sequential
ramp for magnitudes and a two-hue diverging pair for signed values. Every function returns
the file name it wrote, or None when it has nothing to draw.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
OTHER = "#c3c2b7"
# Categorical slots, in order (the validated default palette: adjacent pairs pass CVD and
# normal-vision separation on the light surface).
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
# Diverging poles (cool = up / positive, warm = down / negative), neutral midpoint.
COOL, WARM = "#1c5cab", "#c4501f"
SEQ = LinearSegmentedColormap.from_list("seq", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DPI = 110

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 10, "axes.titleweight": "bold",
    "axes.titlecolor": INK, "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
})


def safe(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_") or "x"


def _save(fig, out_dir: Path, name: str) -> str:
    fig.savefig(out_dir / name, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return name


def _short(text: Any, n: int = 30) -> str:
    text = str(text)
    return text if len(text) <= n else text[: n - 1] + "…"


def condition_colors(levels: list[str]) -> dict[str, str]:
    return {lv: SLOTS[i] if i < len(SLOTS) else OTHER for i, lv in enumerate(levels)}


# ---------------------------------------------------------------- QC


def qc_overview(rows: list[dict[str, Any]], unit: str, out_dir: Path) -> str | None:
    df = pd.DataFrame(rows)
    if df.empty or "pct_removed" not in df:
        return None
    levels = list(dict.fromkeys(df["condition"].dropna().astype(str))) if "condition" in df else []
    df = df.iloc[::-1].reset_index(drop=True)  # first sample on top
    colors = condition_colors(levels)
    col = [colors.get(str(c), SLOTS[0]) for c in df["condition"]] if levels else [SLOTS[0]] * len(df)
    panels = [("pct_removed", f"% of {unit} removed by QC", False)]
    if "median_counts" in df and pd.to_numeric(df["median_counts"], errors="coerce").gt(0).any():
        panels.append(("median_counts", f"Median counts per {unit[:-1]} after QC (log)", True))
    n = len(df)
    fig, axes = plt.subplots(1, len(panels), figsize=(4.0 * len(panels) + 1.2, 0.26 * n + 1.3), sharey=True,
                             squeeze=False)
    y = np.arange(n)
    for ax, (key, title, log) in zip(axes[0], panels):
        vals = pd.to_numeric(df[key], errors="coerce").values
        ax.barh(y, np.nan_to_num(vals), color=col, height=0.66, edgecolor=SURFACE, linewidth=1)
        for yi, v in zip(y, vals):
            if not np.isnan(v):
                label = f" {v:,.0f}" if log else (f" {v:.1f}" if v >= 1 or v == 0 else f" {v:.2f}")
                ax.text(v, yi, label, va="center", fontsize=7, color=INK2)
        ax.set_title(title, loc="left")
        ax.xaxis.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        if log:
            ax.set_xscale("log")
            ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
            ax.set_xlim(np.nanmin(vals) / 1.6, np.nanmax(vals) * 3)
        else:
            ax.set_xlim(0, max(np.nanmax(vals), 0.5) * 1.3)
    axes[0][0].set_yticks(y, df["sample_id"], fontsize=7.5)
    if levels and len(levels) > 1:
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[lv]) for lv in levels]
        levels = [str(lv) for lv in levels]
        fig.legend(handles, levels, loc="lower center", ncols=min(len(levels), 6), bbox_to_anchor=(0.5, -0.02 - 0.12 / max(n, 1)),
                   fontsize=8)
    return _save(fig, out_dir, "qc_overview.png")


# ---------------------------------------------------------------- domains


def domain_coverage(cov: pd.DataFrame | None, cond: dict[str, str], unit: str, out_dir: Path) -> str | None:
    if cov is None or cov.empty:
        return None
    samples = list(cov.columns)
    if cond:
        levels = list(dict.fromkeys(cond.get(s, "") for s in samples))
        samples = sorted(samples, key=lambda s: (levels.index(cond.get(s, "")), samples.index(s)))
    m = cov[samples]
    nd, ns = m.shape
    fig, ax = plt.subplots(figsize=(0.55 * ns + 2.4, 0.3 * nd + 1.6))
    ax.imshow(m.values * 100, aspect="auto", cmap=SEQ, vmin=0, vmax=max(5.0, float(np.nanmax(m.values)) * 100),
              interpolation="nearest")
    for i in range(nd):
        for j in range(ns):
            v = m.values[i, j] * 100
            if v >= 0.5:
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6.5,
                        color=SURFACE if v > 0.6 * max(5.0, float(np.nanmax(m.values)) * 100) else INK2)
    ax.set_xticks(range(ns), samples, rotation=90, fontsize=7.5)
    ax.set_yticks(range(nd), m.index, fontsize=7.5)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    if cond and len(set(cond.get(s, "") for s in samples)) > 1:
        start = 0
        for j in range(1, ns + 1):
            if j == ns or cond.get(samples[j]) != cond.get(samples[start]):
                ax.annotate(cond.get(samples[start], ""), xy=((start + j - 1) / 2, -0.7), ha="center", va="bottom",
                            fontsize=8, fontweight="bold", color=INK, annotation_clip=False)
                if j < ns:
                    ax.axvline(j - 0.5, color=INK, linewidth=1)
                start = j
    ax.set_title(f"% of each sample's {unit} in each domain", loc="left", pad=18)
    return _save(fig, out_dir, "domain_coverage.png")


# ---------------------------------------------------------------- comparison


def de_counts(rows: list[dict[str, Any]], label: str, levels: tuple[str, str] | None, out_dir: Path) -> str | None:
    df = pd.DataFrame([r for r in rows if r.get("status") == "tested"])
    if df.empty or "n_up" not in df:
        return None
    df["n_up"] = pd.to_numeric(df["n_up"], errors="coerce").fillna(0)
    df["n_down"] = pd.to_numeric(df["n_down"], errors="coerce").fillna(0)
    if (df["n_up"] + df["n_down"]).sum() == 0:
        return None
    df["name"] = df["label"].astype(str) + (" · " + df["comparison"].astype(str) if df["comparison"].nunique() > 1 else "")
    df = df.iloc[::-1]
    n = len(df)
    fig, ax = plt.subplots(figsize=(6.2, 0.28 * n + 1.3))
    y = np.arange(n)
    ax.barh(y, df["n_up"], color=COOL, height=0.66, edgecolor=SURFACE, linewidth=1, label="up")
    ax.barh(y, -df["n_down"], color=WARM, height=0.66, edgecolor=SURFACE, linewidth=1, label="down")
    lim = float(max(df["n_up"].max(), df["n_down"].max())) * 1.3 + 1
    for yi, u, d in zip(y, df["n_up"], df["n_down"]):
        if u:
            ax.text(u, yi, f" {int(u):,}", va="center", fontsize=7, color=INK2)
        if d:
            ax.text(-d, yi, f"{int(d):,} ", va="center", ha="right", fontsize=7, color=INK2)
    ax.set_xlim(-lim, lim)
    ax.axvline(0, color=INK2, linewidth=0.8)
    ax.set_yticks(y, [_short(v, 34) for v in df["name"]], fontsize=7.5)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{abs(int(v)):,}"))
    lvl = f" in {levels[0]} vs {levels[1]}" if levels else ""
    ax.set_xlabel(f"significant genes: ← lower{lvl} | higher{lvl} →".replace("  ", " "))
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"Significant DE genes per {label}", loc="left")
    return _save(fig, out_dir, f"de_counts_{safe(label)}.png")


# ---------------------------------------------------------------- composition


def composition_bars(comp: pd.DataFrame | None, title: str, out_dir: Path, name: str, top: int = 7) -> str | None:
    if comp is None or comp.empty:
        return None
    cols = list(comp.columns)
    shown, rest = cols[:top], cols[top:]
    data = comp[shown].copy()
    if rest:
        data["Other"] = comp[rest].sum(axis=1)
    colors = SLOTS[: len(shown)] + ([OTHER] if rest else [])
    n = len(comp)
    fig, ax = plt.subplots(figsize=(max(4.5, 0.45 * n + 3.6), 3.8))
    x = np.arange(n)
    bottom = np.zeros(n)
    for c, color in zip(data.columns, colors):
        vals = data[c].values * 100
        ax.bar(x, vals, bottom=bottom, color=color, width=0.8, edgecolor=SURFACE, linewidth=1.2, label=_short(c, 32))
        bottom += vals
    ax.set_xticks(x, comp.index, rotation=90 if n > 4 else 0, fontsize=7.5)
    ax.set_ylim(0, 100)
    ax.set_ylabel("%")
    ax.yaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=7.5, reverse=True)
    ax.set_title(title, loc="left")
    return _save(fig, out_dir, name)


# ---------------------------------------------------------------- cell-cell communication


def ccc_spatial(spatial: pd.DataFrame | None, label: str, out_dir: Path, top: int = 15) -> str | None:
    if spatial is None or spatial.empty:
        return None
    df = spatial.copy()
    df["abs"] = df["mean_global_morans"].abs()
    df = df.sort_values("abs", ascending=False).head(top).sort_values("mean_global_morans")
    n = len(df)
    fig, ax = plt.subplots(figsize=(5.8, 0.27 * n + 1.3))
    y = np.arange(n)
    vals = df["mean_global_morans"].values
    ax.barh(y, vals, color=[COOL if v >= 0 else WARM for v in vals], height=0.66, edgecolor=SURFACE, linewidth=1)
    for yi, v, k, m in zip(y, vals, df["n_samples_significant"], df["n_samples"]):
        ax.text(v, yi, f" {v:+.2f} ({int(k)}/{int(m)})" if v >= 0 else f"{v:+.2f} ", va="center",
                ha="left" if v >= 0 else "right", fontsize=7, color=INK2)
    ax.axvline(0, color=INK2, linewidth=0.8)
    lim = float(np.nanmax(np.abs(vals))) * 1.45 + 0.02
    ax.set_xlim(-lim if (vals < 0).any() else -0.02, lim)
    ax.set_yticks(y, [f"{a} → {b}" for a, b in zip(df["ligand"], df["receptor"])], fontsize=7.5)
    ax.xaxis.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("bivariate Moran's R, mean over samples (in brackets: samples with R > 0, p ≤ 0.05)")
    ax.set_title("Spatial co-expression of ligand and receptor", loc="left")
    return _save(fig, out_dir, f"ccc_spatial_{safe(label)}.png")
