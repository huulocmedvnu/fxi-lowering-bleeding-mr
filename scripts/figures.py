"""ATVB figures 2-4 (plotly_white, Plotly colours, Arial). 900 px wide, text 13 px (about 7 pt at 174 mm print width).
Inputs: outputs/results_all.tsv, cis-window summary statistics. Writes outputs/figures/ (png, tiff 600 dpi equivalent, svg)."""
import os, numpy as np, pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from PIL import Image
OUT = "outputs/figures"; os.makedirs(OUT, exist_ok=True)
BLUE, RED, GREEN, PURPLE, ORANGE, GREY, INK = "#636EFA", "#EF553B", "#00CC96", "#AB63FA", "#FFA15A", "#8c96a8", "#111111"
FONT = dict(family="Arial", size=13, color=INK); WPX = 900
def save(fig, name, h, l=60, r=20, t=40, b=50):
    fig.update_layout(template="plotly_white", width=WPX, height=h, font=FONT, margin=dict(l=l, r=r, t=t, b=b), showlegend=False)
    fig.update_xaxes(showgrid=False, showline=True, linecolor=INK, ticks="outside")
    fig.update_yaxes(showgrid=False)
    fig.write_image(f"{OUT}/{name}.png", scale=4); fig.write_image(f"{OUT}/{name}.svg")
    im = Image.open(f"{OUT}/{name}.png").convert("RGB"); im.save(f"{OUT}/{name}.tiff", compression="tiff_lzw", dpi=(600, 600))
d = pd.read_csv("outputs/results_all.tsv", sep="\t")
def row(analysis, src, ep): return d[(d.analysis == analysis) & (d.source == src) & (d.outcome == ep)].iloc[0]
# ---------- Figure 2: forest by outcome and data source
items = [("Venous thromboembolism (positive control)", None, None, None),
         ("  FinnGen, venous thromboembolism", "IVW r2<0.1", "FinnGen", "I9_VTE"),
         ("  UK Biobank, pulmonary embolism", "IVW r2<0.1", "UKB", "icd10-I26-both_sexes"),
         ("  UK Biobank, phlebitis and thrombophlebitis", "IVW r2<0.1", "UKB", "icd10-I80-both_sexes"),
         ("  Million Veteran Program, venous thromboembolism", "IVW r2<0.1", "MVP", "VTE"),
         ("Postpartum hemorrhage", None, None, None),
         ("  FinnGen", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HEAMORRH_ALLW"),
         ("  UK Biobank", "IVW r2<0.1", "UKB", "icd10-O72-both_sexes"),
         ("  Meta-analysis", "Fixed-effect meta FinnGen+UKB", "META", "PPH"),
         ("Heavy menstrual bleeding", None, None, None),
         ("  FinnGen", "IVW r2<0.1", "FinnGen", "N14_MENORRHAGIA"),
         ("  UK Biobank", "IVW r2<0.1", "UKB", "icd10-N92-both_sexes"),
         ("  Meta-analysis", "Fixed-effect meta FinnGen+UKB", "META", "HMB"),
         ("Other outcomes, FinnGen", None, None, None),
         ("  Postpartum hemorrhage, uterine atony", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HAEMORRH_ANATOMY"),
         ("  Postpartum hemorrhage, retained placenta", "IVW r2<0.1", "FinnGen", "O15_POSTPART_HAEMORRH_PLACENTA"),
         ("  Antepartum hemorrhage", "IVW r2<0.1", "FinnGen", "O15_ANTEPART_HAEM_NAS"),
         ("  Other abnormal uterine bleeding", "IVW r2<0.1", "FinnGen", "N14_ABNORMALBLEED"),
         ("  Any bleeding", "IVW r2<0.1", "FinnGen", "BLEEDING"),
         ("  Intracerebral hemorrhage", "IVW r2<0.1", "FinnGen", "I9_ICH"),
         ("  Iron deficiency anemia due to blood loss", "IVW r2<0.1", "FinnGen", "D3_ANAEMIA_IRONDEF_BLOODLOSS"),
         ("  Acute posthemorrhagic anemia", "IVW r2<0.1", "FinnGen", "D3_ACUTEPOSTBLEEDANAEMIA"),
         ("  Gestational hypertension (negative control)", "IVW r2<0.1", "FinnGen", "O15_GESTAT_HYPERT")]
n = len(items); y = list(range(n, 0, -1))
fig = make_subplots(rows=1, cols=2, column_widths=[0.62, 0.38], horizontal_spacing=0.02, specs=[[{"type": "xy"}, {"type": "xy"}]])
labels, txt = [], []
for (lab, a, s, ep), yy in zip(items, y):
    labels.append(f"<b>{lab}</b>" if a is None else lab)
    if a is None: txt.append(""); continue
    r = row(a, s, ep); col = BLUE if "thrombo" in lab.lower() or s == "MVP" or ep in ("icd10-I26-both_sexes", "icd10-I80-both_sexes", "I9_VTE") else (RED if s == "META" else INK)
    fig.add_trace(go.Scatter(x=[r.lo, r.hi], y=[yy, yy], mode="lines", line=dict(color=col, width=2)), 1, 2)
    fig.add_trace(go.Scatter(x=[r.OR], y=[yy], mode="markers", marker=dict(symbol="diamond" if s == "META" else "square", size=11 if s == "META" else 8, color=col)), 1, 2)
    p = f"{r.P:.1e}".replace("e-0", " × 10<sup>-").replace("e-", " × 10<sup>-") + "</sup>" if r.P < 0.001 else f"{r.P:.2f}"
    txt.append(f"{r.OR:.2f} ({r.lo:.2f} to {r.hi:.2f})   P = {p}")
for lab, t_, yy in zip(labels, txt, y):
    fig.add_annotation(x=0, y=yy, xref="x1", yref="y1", text=lab, showarrow=False, xanchor="left", font=FONT)
    fig.add_annotation(x=0.66, y=yy, xref="x1", yref="y1", text=t_, showarrow=False, xanchor="left", font=FONT)
fig.add_annotation(x=0.66, y=n + 1.2, xref="x1", yref="y1", text="<b>OR (95% CI)</b>", showarrow=False, xanchor="left", font=FONT)
fig.add_vline(x=1, line=dict(color=GREY, dash="dot"), row=1, col=2)
fig.update_xaxes(range=[0, 1], visible=False, row=1, col=1); fig.update_yaxes(range=[0.3, n + 1.6], visible=False, row=1, col=1)
fig.update_yaxes(range=[0.3, n + 1.6], visible=False, row=1, col=2)
fig.update_xaxes(type="log", tickvals=[0.4, 0.5, 0.7, 1, 1.4], ticktext=["0.4", "0.5", "0.7", "1.0", "1.4"], range=[np.log10(0.38), np.log10(1.55)], title_text="Odds ratio per SD lower plasma factor XI", row=1, col=2)
save(fig, "Figure2", 760, l=10, r=20, t=20, b=55)
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
    fig.add_trace(go.Scattergl(x=x / 1e6, y=yv, mode="markers", marker=dict(size=4, color=c, opacity=0.8)), i, 1)
    fig.add_vline(x=lead / 1e6, line=dict(color=GREY, dash="dot", width=1), row=i, col=1)
    fig.update_yaxes(title_text="-log<sub>10</sub> P", title_font=dict(size=12), row=i, col=1)
    fig.add_annotation(x=0.01, y=0.98, xref=f"x{i} domain" if i > 1 else "x domain", yref=f"y{i} domain" if i > 1 else "y domain", text=f"<b>{lab}</b>" + (f"<br>{pp}" if pp else ""), showarrow=False, xanchor="left", yanchor="top", font=dict(size=12, family="Arial", color=INK), align="left")
g = len(tracks) + 1
for name, s, e in [("KLKB1", 186227000, 186258471), ("F11", 186265524, 186289681)]:
    fig.add_trace(go.Scatter(x=[s / 1e6, e / 1e6], y=[0.5, 0.5], mode="lines", line=dict(color=INK, width=8)), g, 1)
    fig.add_annotation(x=(s if name == "KLKB1" else e) / 1e6, y=0.5, xref=f"x{g}", yref=f"y{g}", text=f"<i>{name}</i>", showarrow=False, xanchor="right" if name == "KLKB1" else "left", xshift=-4 if name == "KLKB1" else 4, font=dict(size=12, family="Arial"))
fig.update_yaxes(visible=False, range=[0, 1], row=g, col=1)
fig.update_xaxes(title_text="Chromosome 4 position, GRCh38 (Mb)", row=g, col=1)
save(fig, "Figure3", 900, l=70, r=20, t=20, b=55)
# ---------- Figure 4: sensitivity analyses
meths = ["IVW r2<0.1", "IVW r2<0.01", "Weighted median r2<0.01", "MR-Egger r2<0.01", "Wald lead variant", "IVW excluding KLKB1-associated", "MVMR adjusted for KLKB1"]
mlab = ["IVW, instruments r<sup>2</sup> < 0.1 (primary)", "IVW, instruments r<sup>2</sup> < 0.01", "Weighted median", "MR-Egger", "Lead variant rs3756008", "IVW excluding variants associated with prekallikrein", "Multivariable, adjusted for prekallikrein"]
outs = [("I9_VTE", "A  Venous thromboembolism", BLUE, [0.35, 0.75]), ("O15_POSTPART_HEAMORRH_ALLW", "B  Postpartum hemorrhage", RED, [0.7, 1.45]), ("N14_MENORRHAGIA", "C  Heavy menstrual bleeding", ORANGE, [0.7, 1.45])]
fig = make_subplots(rows=1, cols=4, column_widths=[0.34, 0.22, 0.22, 0.22], horizontal_spacing=0.025)
yy = list(range(len(meths), 0, -1))
for t, yv in zip(mlab, yy): fig.add_annotation(x=0, y=yv, xref="x1", yref="y1", text=t, showarrow=False, xanchor="left", font=FONT)
for j, (ep, title, c, rg) in enumerate(outs, 2):
    for m_, yv in zip(meths, yy):
        r = row(m_, "FinnGen", ep)
        fig.add_trace(go.Scatter(x=[r.lo, r.hi], y=[yv, yv], mode="lines", line=dict(color=c, width=2)), 1, j)
        fig.add_trace(go.Scatter(x=[r.OR], y=[yv], mode="markers", marker=dict(symbol="square", size=8, color=c)), 1, j)
    fig.add_vline(x=1, line=dict(color=GREY, dash="dot"), row=1, col=j)
    fig.update_xaxes(type="log", range=[np.log10(rg[0]), np.log10(rg[1])], tickvals=[0.4, 0.5, 0.6, 0.8, 1, 1.2, 1.4], title_text="OR per SD lower factor XI", title_font=dict(size=12), row=1, col=j)
    fig.add_annotation(x=0.5, y=len(meths) + 0.9, xref=f"x{j} domain", yref=f"y{j}", text=f"<b>{title}</b>", showarrow=False, font=FONT)
    fig.update_yaxes(range=[0.4, len(meths) + 1.3], visible=False, row=1, col=j)
fig.update_xaxes(range=[0, 1], visible=False, row=1, col=1); fig.update_yaxes(range=[0.4, len(meths) + 1.3], visible=False, row=1, col=1)
save(fig, "Figure4", 380, l=10, r=20, t=15, b=60)
print("figures written")
