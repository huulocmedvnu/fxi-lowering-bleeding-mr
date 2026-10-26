"""Single source of all reported estimates (OR per SD lower plasma FXI). Writes outputs/results_all.tsv.
Instruments: cis F11 pQTL P<5e-8 clumped r2<0.1 (primary) or r2<0.01; positions deduplicated before matching."""
import os, subprocess, numpy as np, pandas as pd
from scipy import stats
PL = os.environ.get("PLINK", "plink"); BF = os.environ.get("REF_EUR", "ref/EUR")  # PLINK 1.9, 1000G EUR bfile prefix
bim = pd.read_csv(BF + ".bim", sep="\t", header=None, names=["c", "rsid", "cm", "pos37", "A1", "A2"])
def pq(name, hdrfile):
    h = open(hdrfile).read().strip().split("\t")
    d = pd.read_csv(f"data/pqtl/{name}.cis.tsv", sep="\t", header=None, names=h).rename(columns={"base_pair_location": "pos", "effect_allele": "ea", "other_allele": "oa", "beta": "b", "standard_error": "s", "p_value": "p"})
    return d[d.rsid.astype(str).str.startswith("rs")].drop_duplicates("pos").drop_duplicates("rsid")
fx = pq("F11", "data/pqtl/header.txt"); kl = pq("KLKB1", "data/pqtl/KLKB1.header.txt")
def instruments(r2):
    iv = pd.read_csv(f"outputs/clump_r{r2}.clumped", sep=r"\s+").SNP
    d = fx[fx.rsid.isin(iv)].copy(); fl = d.b > 0
    d["low"] = np.where(fl, d.oa, d.ea); d["oth"] = np.where(fl, d.ea, d.oa); d["bx"] = -d.b.abs(); d["sx"] = d.s
    k = kl.set_index("pos"); d["bk"] = [(k.loc[p, "b"] if k.loc[p, "ea"] == lw else -k.loc[p, "b"]) if p in k.index else np.nan for p, lw in zip(d.pos, d.low)]
    d["kp"] = [k.loc[p, "p"] if p in k.index else np.nan for p in d.pos]
    d = d.merge(bim[["rsid", "pos37"]], on="rsid", how="left")
    pd.Series(d.rsid).to_csv("outputs/tmp_iv.txt", index=False, header=False)
    subprocess.run([PL, "--bfile", BF, "--extract", "outputs/tmp_iv.txt", "--r", "square", "--keep-allele-order", "--write-snplist", "--out", "outputs/tmp_ld"], capture_output=True)
    s = open("outputs/tmp_ld.snplist").read().split(); L = np.loadtxt("outputs/tmp_ld.ld").reshape(len(s), len(s))
    b2 = bim.set_index("rsid").loc[s]; lo = dict(zip(d.rsid, d.low)); g = np.array([1 if b2.loc[x, "A1"] == lo[x] else -1 for x in s])
    d = d.set_index("rsid").loc[s].reset_index(); return d, s, L * np.outer(g, g)
def load_outcome(src, ep):
    if src == "FinnGen":
        o = pd.read_csv(f"data/finngen/{ep}.cis.tsv", sep="\t").rename(columns={"#chrom": "chr"}); return o.drop_duplicates("pos"), "pos", "ref", "alt", "beta", "sebeta"
    if src == "UKB":
        h = open(f"data/ukb/{ep}.header").read().strip().split("\t"); o = pd.read_csv(f"data/ukb/{ep}.cis.tsv", sep="\t", header=None, names=h)
        o = o[o.low_confidence_EUR.astype(str).str.lower() == "false"].rename(columns={"pos": "pos37"}).drop_duplicates("pos37"); return o, "pos37", "ref", "alt", "beta_EUR", "se_EUR"
    o = pd.read_csv("data/mvp/VTE.cis.tsv", sep="\t", na_values=["#NA"]).rename(columns={"base_pair_location": "pos"})
    o["beta"] = np.log(o.odds_ratio); o["se"] = (np.log(o.ci_upper) - np.log(o.ci_lower)) / 3.92
    return o.drop(columns=["rsid"]).dropna(subset=["beta", "se"]).drop_duplicates("pos"), "pos", "other_allele", "effect_allele", "beta", "se"
def matched(ins, src, ep):
    o, pk, ref, alt, be, se = load_outcome(src, ep); m = ins.merge(o, on=pk)
    sg = np.where((m[alt] == m.low) & (m[ref] == m.oth), 1, np.where((m[ref] == m.low) & (m[alt] == m.oth), -1, 0)); k = sg != 0
    m = m[k].copy(); m["by"] = m[be].values * sg[k]; m["sy"] = m[se].values; return m
def givw(m, snps, LD, X=None):
    idx = [snps.index(s) for s in m.rsid]; R = LD[np.ix_(idx, idx)]; Oi = np.linalg.pinv(np.outer(m.sy, m.sy) * R)
    X = m[["bx"]].values if X is None else X; A = X.T @ Oi @ X; b = np.linalg.solve(A, X.T @ Oi @ m.by.values); r = m.by.values - X @ b
    Q = r @ Oi @ r; phi = max(1, Q / (len(m) - X.shape[1])); se = np.sqrt(np.diag(np.linalg.inv(A)) * phi); return b[0], se[0], len(m), Q
def wm(m):
    est = (m.by / m.bx).values; w = (m.bx / m.sy).values ** 2; o = np.argsort(est); c = np.cumsum(w[o] / w.sum()) - 0.5 * w[o] / w.sum(); b = np.interp(0.5, c, est[o])
    rng = np.random.default_rng(1); bs = []
    for _ in range(2000):
        bx = rng.normal(m.bx, m.sx); by = rng.normal(m.by, m.sy); e = by / bx; ww = (bx / m.sy) ** 2; oo = np.argsort(e); cc = np.cumsum(ww.values[oo] / ww.sum()) - 0.5 * ww.values[oo] / ww.sum(); bs.append(np.interp(0.5, cc, e[oo]))
    return b, np.std(bs)
def egger(m):
    X = np.column_stack([np.ones(len(m)), -m.bx]); W = np.diag(1 / m.sy**2); y = -m.by.values
    beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ y); r = y - X @ beta; s2 = max(1, (r @ W @ r) / (len(m) - 2)); se = np.sqrt(np.diag(np.linalg.inv(X.T @ W @ X)) * s2)
    return beta[1], se[1], 2 * stats.norm.sf(abs(beta[0] / se[0]))
rows = []
def add(analysis, src, ep, b, se, n, extra=""):
    rows.append(dict(analysis=analysis, source=src, outcome=ep, n_iv=n, b=b, se=se, OR=np.exp(-b), lo=np.exp(-b - 1.96 * se), hi=np.exp(-b + 1.96 * se), P=2 * stats.norm.sf(abs(b / se)), note=extra))
fg = open("data/outcomes.txt").read().split()
ins1, s1, L1 = instruments("0.1"); ins2, s2, L2 = instruments("0.01")
for ep in fg:
    m = matched(ins1, "FinnGen", ep); b, se, n, Q = givw(m, s1, L1); add("IVW r2<0.1", "FinnGen", ep, b, se, n, f"Q_P={stats.chi2.sf(Q, n-1):.2g}")
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA"]:
    m2 = matched(ins2, "FinnGen", ep); b, se, n, Q = givw(m2, s2, L2); add("IVW r2<0.01", "FinnGen", ep, b, se, n)
    b, se = wm(m2); add("Weighted median r2<0.01", "FinnGen", ep, b, se, len(m2))
    b, se, ip = egger(m2); add("MR-Egger r2<0.01", "FinnGen", ep, b, se, len(m2), f"intercept_P={ip:.2g}")
    m = matched(ins1, "FinnGen", ep); j = m.bx.abs().div(m.sx).idxmax(); r = m.loc[j]; add("Wald lead variant", "FinnGen", ep, r.by / r.bx, r.sy / abs(r.bx), 1, r.rsid)
    keep = m[m.kp >= 1e-4]; b, se, n, _ = givw(keep, s1, L1); add("IVW excluding KLKB1-associated", "FinnGen", ep, b, se, n)
    mk = m.dropna(subset=["bk"]); b, se, n, _ = givw(mk, s1, L1, X=mk[["bx", "bk"]].values); add("MVMR adjusted for KLKB1", "FinnGen", ep, b, se, n)
for ep in ["icd10-O72-both_sexes", "icd10-N92-both_sexes", "icd10-I26-both_sexes", "icd10-I80-both_sexes"]:
    m = matched(ins1, "UKB", ep); b, se, n, _ = givw(m, s1, L1); add("IVW r2<0.1", "UKB", ep, b, se, n)
m = matched(ins1, "MVP", "VTE"); b, se, n, _ = givw(m, s1, L1); add("IVW r2<0.1", "MVP", "VTE", b, se, n)
d = pd.DataFrame(rows)
for a, bb, lab in [("O15_POSTPART_HEAMORRH_ALLW", "icd10-O72-both_sexes", "PPH"), ("N14_MENORRHAGIA", "icd10-N92-both_sexes", "HMB")]:
    x = d[(d.analysis == "IVW r2<0.1") & (d.outcome.isin([a, bb]))]; w = 1 / x.se**2; bm = (x.b * w).sum() / w.sum(); sm = np.sqrt(1 / w.sum())
    q = ((x.b - bm)**2 * w).sum(); add("Fixed-effect meta FinnGen+UKB", "META", lab, bm, sm, 0, f"het_P={stats.chi2.sf(q,1):.2g}")
d = pd.DataFrame(rows); d.to_csv("outputs/results_all.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(d[["analysis", "source", "outcome", "n_iv", "OR", "lo", "hi", "P", "note"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
