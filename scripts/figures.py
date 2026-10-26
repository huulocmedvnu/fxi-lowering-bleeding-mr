"""ATVB figures 2-4 (plotly_white, Plotly colours, Arial). 900 px wide = 174 mm print width, text 19 px (about 10.4 pt; AHA minimum 10 pt).
Inputs: outputs/results_all.tsv, cis-window summary statistics. Writes outputs/figures/ (png, tiff 600 dpi equivalent, svg)."""
import os, numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image
OUT = "outputs/figures"; os.makedirs(OUT, exist_ok=True)
BLUE, RED, GREEN, PURPLE, ORANGE, GREY, INK = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#8c96a8", "#111111"
FONT = dict(family="Arial", size=19, color=INK); WPX = 900
def save(fig, name, h, l=60, r=20, t=40, b=50):
    fig.update_layout(template="plotly_white", width=WPX, height=h, font=FONT, margin=dict(l=l, r=r, t=t, b=b), showlegend=False)
    fig.update_xaxes(showgrid=False, showline=True, linecolor=INK, ticks="outside")
    fig.update_yaxes(showgrid=False)
    fig.write_image(f"{OUT}/{name}.png", scale=4); fig.write_image(f"{OUT}/{name}.svg")
    im = Image.open(f"{OUT}/{name}.png").convert("RGB"); im.save(f"{OUT}/{name}.tiff", compression="tiff_lzw", dpi=(600, 600))
d = pd.read_csv("outputs/results_all.tsv", sep="\t")
def row(analysis, src, ep): return d[(d.analysis == analysis) & (d.source == src) & (d.outcome == ep)].iloc[0]
# ---------- Figure 2: forest by outcome and data source
items = [("Venous thromboembolism", None, None, None),
         ("  FinnGen", "IVW r2<0.1", "FinnGen", "I9_VTE"),
         ("  UK Biobank, pulmonary embolism", "IVW r2<0.1", "UKB", "icd10-I26-both_sexes"),
         ("  UK Biobank, thrombophlebitis", "IVW r2<0.1", "UKB", "icd10-I80-both_sexes"),
         ("  Million Veteran Program", "IVW r2<0.1", "MVP", "VTE"),
         ("Postpartum hemorrhage", None, None, None),
         ("  FinnGen", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HEAMORRH_ALLW"),
         ("  UK Biobank", "IVW r2<0.1", "UKB", "icd10-O72-both_sexes"),
         ("  Meta-analysis", "Fixed-effect meta FinnGen+UKB", "META", "PPH"),
         ("Heavy menstrual bleeding", None, None, None),
         ("  FinnGen", "IVW r2<0.1", "FinnGen", "N14_MENORRHAGIA"),
         ("  UK Biobank", "IVW r2<0.1", "UKB", "icd10-N92-both_sexes"),
         ("  Meta-analysis", "Fixed-effect meta FinnGen+UKB", "META", "HMB"),
         ("Other outcomes in FinnGen", None, None, None),
         ("  PPH, uterine atony", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HAEMORRH_ANATOMY"),
         ("  PPH, retained placenta", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HAEMORRH_PLACENTA"),
         ("  Antepartum hemorrhage", "IVW r2<0.1", "FinnGen", "O15_ANTEPART_HAEM_NAS"),
         ("  Other abnormal uterine bleeding", "IVW r2<0.1", "FinnGen", "N14_ABNORMALBLEED"),
         ("  Any bleeding", "IVW r2<0.1", "FinnGen", "BLEEDING"),
         ("  Intracerebral hemorrhage", "IVW r2<0.1", "FinnGen", "I9_ICH"),
         ("  Iron deficiency anemia, blood loss", "IVW r2<0.1", "FinnGen", "D3_ANAEMIA_IRONDEF_BLOODLOSS"),
         ("  Acute posthemorrhagic anemia", "IVW r2<0.1", "FinnGen", "D3_ACUTEPOSTBLEEDANAEMIA"),
         ("  Gestational hypertension", "IVW r2<0.1", "FinnGen", "O15_GESTAT_HYPERT")]
n = len(items); y = list(range(n, 0, -1))
fig = make_subplots(rows=1, cols=2, column_widths=[0.76, 0.24], horizontal_spacing=0.02, specs=[[{"type": "xy"}, {"type": "xy"}]])
labels, txt = [], []
for (lab, a, s, ep), yy in zip(items, y):
    labels.append(f"<b>{lab}</b>" if a is None else lab)
    if a is None: txt.append(("", "")); continue
    r = row(a, s, ep); col = BLUE if s == "MVP" or ep in ("icd10-I26-both_sexes", "icd10-I80-both_sexes", "I9_VTE") else (RED if s == "META" else INK)
    fig.add_trace(go.Scatter(x=[r.lo, r.hi], y=[yy, yy], mode="lines", line=dict(color=col, width=2.5)), 1, 2)
    fig.add_trace(go.Scatter(x=[r.OR], y=[yy], mode="markers", marker=dict(symbol="diamond" if s == "META" else "square", size=15 if s == "META" else 11, color=col)), 1, 2)
    p = f"{r.P:.1e}".replace("e-0", " × 10<sup>-").replace("e-", " × 10<sup>-") + "</sup>" if r.P < 0.001 else (f"{r.P:#.2g}" if r.P >= 0.01 else f"{r.P:.1g}")  # AHA: 2 significant digits, 1 below 0.01
    txt.append((f"{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})", p))
for lab, t_, yy in zip(labels, txt, y):
    fig.add_annotation(x=0, y=yy, xref="x1", yref="y1", text=lab, showarrow=False, xanchor="left", font=FONT)
    fig.add_annotation(x=0.585, y=yy, xref="x1", yref="y1", text=t_[0], showarrow=False, xanchor="left", font=FONT)
    fig.add_annotation(x=0.85, y=yy, xref="x1", yref="y1", text=t_[1], showarrow=False, xanchor="left", font=FONT)
fig.add_annotation(x=0.585, y=n + 1.2, xref="x1", yref="y1", text="<b>OR (95% CI)</b>", showarrow=False, xanchor="left", font=FONT)
fig.add_annotation(x=0.85, y=n + 1.2, xref="x1", yref="y1", text="<b><i>P</i></b>", showarrow=False, xanchor="left", font=FONT)
fig.add_vline(x=1, line=dict(color=GREY, dash="dot"), row=1, col=2)
fig.update_xaxes(range=[0, 1], visible=False, row=1, col=1); fig.update_yaxes(range=[0.3, n + 1.6], visible=False, row=1, col=1)
fig.update_yaxes(range=[0.3, n + 1.6], visible=False, row=1, col=2)
fig.update_xaxes(type="log", tickvals=[0.4, 0.7, 1, 1.4], ticktext=["0.4", "0.7", "1.0", "1.4"], range=[np.log10(0.38), np.log10(1.55)], title_text="OR per SD lower<br>factor XI", row=1, col=2)
save(fig, "Figure2", 1060, l=10, r=15, t=20, b=85)
# ---------- Figure 3: regional plot
def mlp(p): return -np.log10(np.clip(p.astype(float), 1e-300, 1))
h = open("data/pqtl/header.txt").read().strip().split("\t"); fx = pd.read_csv("data/pqtl/F11.cis.tsv", sep="\t", header=None, names=h)
hk = open("data/pqtl/KLKB1.header.txt").read().strip().split("\t"); kl = pd.read_csv("data/pqtl/KLKB1.cis.tsv", sep="\t", header=None, names=hk)
tracks = [("Plasma factor XI (UK Biobank)", fx.base_pair_location, mlp(fx.p_value), BLUE, ""),
          ("Plasma prekallikrein, KLKB1 (UK Biobank)", kl.base_pair_location, mlp(kl.p_value), GREY, "Shared with factor XI: PP4 0.000 (PP3 1.000)")]
for ep, lab, colr, pp in [("I9_VTE", "Venous thromboembolism (FinnGen)", PURPLE, "Shared with factor XI: PP4 0.999"),
                          ("O15_POSTPART_HEAMORRH_ALLW", "Postpartum hemorrhage (FinnGen)", RED, "Shared with factor XI: PP4 0.025"),
                          ("N14_MENORRHAGIA", "Heavy menstrual bleeding (FinnGen)", ORANGE, "Shared with factor XI: PP4 0.005")]:
    o = pd.read_csv(f"data/finngen/{ep}.cis.tsv", sep="\t"); tracks.append((lab, o.pos, mlp(o.pval), colr, pp))
fig = make_subplots(rows=len(tracks) + 1, cols=1, shared_xaxes=True, vertical_spacing=0.025, row_heights=[1] * len(tracks) + [0.35])
lead = 186264231
for i, (lab, x, yv, c, pp) in enumerate(tracks, 1):
    fig.add_trace(go.Scatter(x=x / 1e6, y=yv, mode="markers", marker=dict(size=5, color=c, opacity=0.8)), i, 1)
    fig.add_vline(x=lead / 1e6, line=dict(color=GREY, dash="dot", width=1), row=i, col=1)
    fig.update_yaxes(title_text="-log<sub>10</sub> P", title_font=dict(size=19), row=i, col=1)
    fig.add_annotation(x=0.01, y=0.98, xref=f"x{i} domain" if i > 1 else "x domain", yref=f"y{i} domain" if i > 1 else "y domain", text=f"<b>{lab}</b>" + (f"<br>{pp}" if pp else ""), showarrow=False, xanchor="left", yanchor="top", font=FONT, align="left", bgcolor="rgba(255,255,255,0.85)")
g = len(tracks) + 1
for name, s, e in [("KLKB1", 186227000, 186258471), ("F11", 186265524, 186289681)]:
    fig.add_trace(go.Scatter(x=[s / 1e6, e / 1e6], y=[0.5, 0.5], mode="lines", line=dict(color=INK, width=8)), g, 1)
    fig.add_annotation(x=(s if name == "KLKB1" else e) / 1e6, y=0.5, xref=f"x{g}", yref=f"y{g}", text=f"<i>{name}</i>", showarrow=False, xanchor="right" if name == "KLKB1" else "left", xshift=-4 if name == "KLKB1" else 4, font=FONT)
fig.update_yaxes(visible=False, range=[0, 1], row=g, col=1)
fig.update_xaxes(title_text="Chromosome 4 position, GRCh38 (Mb)", row=g, col=1)
save(fig, "Figure3", 1180, l=80, r=20, t=20, b=70)
# ---------- Figure 4: sensitivity analyses, one panel per outcome
meths = ["IVW r2<0.1", "IVW r2<0.01", "Weighted median r2<0.01", "MR-Egger r2<0.01", "Wald lead variant", "IVW excluding KLKB1-associated", "MVMR adjusted for KLKB1"]
mlab = ["IVW, r<sup>2</sup> < 0.1 (primary)", "IVW, r<sup>2</sup> < 0.01", "Weighted median, r<sup>2</sup> < 0.01", "MR-Egger, r<sup>2</sup> < 0.01", "Lead variant rs3756008", "Excluding prekallikrein variants", "Adjusted for prekallikrein"]
outs = [("I9_VTE", "A  Venous thromboembolism", BLUE, [0.38, 0.75]), ("O15_POSTPART_HEAMORRH_ALLW", "B  Postpartum hemorrhage", RED, [0.75, 1.45]), ("N14_MENORRHAGIA", "C  Heavy menstrual bleeding", ORANGE, [0.75, 1.45])]
fig = make_subplots(rows=3, cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.02, vertical_spacing=0.09)
yy = list(range(len(meths), 0, -1))
for k, (ep, title, c, rg) in enumerate(outs, 1):
    xa = "x" if k == 1 else f"x{2 * k - 1}"; ya = "y" if k == 1 else f"y{2 * k - 1}"
    for t, yv in zip(mlab, yy): fig.add_annotation(x=0, y=yv, xref=xa, yref=ya, text=t, showarrow=False, xanchor="left", font=FONT)
    fig.add_annotation(x=0, y=len(meths) + 1.1, xref=xa, yref=ya, text=f"<b>{title}</b>", showarrow=False, xanchor="left", font=FONT)
    for m_, yv in zip(meths, yy):
        r = row(m_, "FinnGen", ep)
        fig.add_trace(go.Scatter(x=[r.lo, r.hi], y=[yv, yv], mode="lines", line=dict(color=c, width=2.5)), k, 2)
        fig.add_trace(go.Scatter(x=[r.OR], y=[yv], mode="markers", marker=dict(symbol="square", size=11, color=c)), k, 2)
    fig.add_vline(x=1, line=dict(color=GREY, dash="dot"), row=k, col=2)
    fig.update_xaxes(type="log", range=[np.log10(rg[0]), np.log10(rg[1])], tickvals=[0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4], row=k, col=2)
    fig.update_xaxes(range=[0, 1], visible=False, row=k, col=1)
    for cc in (1, 2): fig.update_yaxes(range=[0.4, len(meths) + 1.5], visible=False, row=k, col=cc)
fig.update_xaxes(title_text="OR per SD lower plasma factor XI", row=3, col=2)
save(fig, "Figure4", 1150, l=10, r=20, t=15, b=75)
print("figures written")
