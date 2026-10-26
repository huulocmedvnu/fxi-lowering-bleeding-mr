"""Follow-up: coloc ABF, KLKB1 (adjacent gene) checks, replication exposure, weighted median / Egger.
OR per SD lower FXI throughout. Run from repo root after mr_fxi.py."""
import subprocess, numpy as np, pandas as pd
from scipy import stats
PL = "/mnt/hdd2/loc-tran/miniconda3/envs/ldsc3/bin/plink"; BF = "/mnt/hdd2/loc-tran/birth_canal_gen/ref/EUR"
def pqtl(name):
    h = open(f"data/pqtl/{'header.txt' if name=='F11' else name+'.header.txt'}").read().strip().split("\t")
    d = pd.read_csv(f"data/pqtl/{name}.cis.tsv", sep="\t", header=None, names=h)
    d = d.rename(columns={"base_pair_location": "pos", "effect_allele": "ea", "other_allele": "oa", "beta": "b", "standard_error": "s", "p_value": "p"})
    return d[d.rsid.astype(str).str.startswith("rs")][["rsid", "pos", "ea", "oa", "b", "s", "p"]].drop_duplicates("pos")
def outcome(ep):
    d = pd.read_csv(f"data/finngen/{ep}.cis.tsv", sep="\t").rename(columns={"#chrom": "chr"})
    return d[["pos", "ref", "alt", "beta", "sebeta", "pval"]].drop_duplicates("pos")
def align(a, o, ea="ea", oa="oa"):
    """effect of o on a's effect allele"""
    m = a.merge(o, on="pos")
    sgn = np.where((m.alt == m[ea]) & (m.ref == m[oa]), 1, np.where((m.ref == m[ea]) & (m.alt == m[oa]), -1, 0))
    m = m[sgn != 0].copy(); m["by"] = m.beta * sgn[sgn != 0]; return m
def labf(z, V, W): r = W / (V + W); return 0.5 * np.log(1 - r) + 0.5 * z * z * r
def coloc(b1, s1, W1, b2, s2, W2, p1=1e-4, p2=1e-4, p12=1e-5):
    l1, l2 = labf(b1 / s1, s1**2, W1), labf(b2 / s2, s2**2, W2)
    lse = lambda x: x.max() + np.log(np.exp(x - x.max()).sum())
    l_1, l_2, l_12 = lse(l1), lse(l2), lse(l1 + l2)
    h = np.array([0, np.log(p1) + l_1, np.log(p2) + l_2, np.log(p1) + np.log(p2) + (l_1 + l_2 + np.log1p(-np.exp(l_12 - l_1 - l_2)) if l_12 - l_1 - l_2 < 0 else l_1 + l_2 - 50), np.log(p12) + l_12])
    pp = np.exp(h - lse(h)); return pp
fx, kl, rep = pqtl("F11"), pqtl("KLKB1"), pqtl("F11rep")
print("== coloc (cis window, prior p12 1e-5): PP3 different variants, PP4 shared")
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA", "O15_GESTAT_HYPERT"]:
    m = align(fx, outcome(ep)); pp = coloc(m.b.values, m.s.values, 0.15**2, m.by.values, m.sebeta.values, 0.2**2)
    print(f"  FXI x {ep:28s} n={len(m)} PP0-4 " + " ".join(f"{x:.3f}" for x in pp))
m = fx.merge(kl, on=["pos"], suffixes=("", "_k")); m = m[(m.ea == m.ea_k) & (m.oa == m.oa_k) | (m.ea == m.oa_k) & (m.oa == m.ea_k)].copy()
m["bk"] = np.where(m.ea == m.ea_k, m.b_k, -m.b_k)
pp = coloc(m.b.values, m.s.values, 0.15**2, m.bk.values, m.s_k.values, 0.15**2); print("  FXI x KLKB1 PP0-4 " + " ".join(f"{x:.3f}" for x in pp))
# instruments (r2<0.1 clump from mr_fxi.py), oriented to FXI-lowering allele
iv = pd.read_csv("outputs/clump_r0.1.clumped", sep=r"\s+").SNP
ins = fx[fx.rsid.isin(iv)].copy(); fl = ins.b > 0
ins["low"] = np.where(fl, ins.oa, ins.ea); ins["oth"] = np.where(fl, ins.ea, ins.oa); ins["bx"] = -ins.b.abs(); ins["sx"] = ins.s
k = ins.merge(kl[["pos", "ea", "oa", "b", "s", "p"]].rename(columns={"ea": "kea", "oa": "koa", "b": "kb", "s": "ks", "p": "kp"}), on="pos")
k["bk"] = np.where(k.kea == k.low, k.kb, -k.kb)
print(f"== KLKB1: instruments with KLKB1 P<1e-4: {(k.kp<1e-4).sum()}/{len(k)}; corr(bx, bk) {np.corrcoef(k.bx, k.bk)[0,1]:.2f}")
# LD
def ldmat(snps, low):
    pd.Series(snps).to_csv("outputs/tmp_iv.txt", index=False, header=False)
    subprocess.run([PL, "--bfile", BF, "--extract", "outputs/tmp_iv.txt", "--r", "square", "--keep-allele-order", "--write-snplist", "--out", "outputs/tmp_ld"], capture_output=True)
    s = open("outputs/tmp_ld.snplist").read().split(); L = np.loadtxt("outputs/tmp_ld.ld").reshape(len(s), len(s))
    bim = pd.read_csv(BF + ".bim", sep="\t", header=None, names=["c", "SNP", "cm", "bp", "A1", "A2"]).set_index("SNP").loc[s]
    lo = dict(zip(snps, low)); g = np.array([1 if bim.loc[x, "A1"] == lo[x] else -1 for x in s]); return s, L * np.outer(g, g)
def givw(bx, by, sy, R):
    Oi = np.linalg.pinv(np.outer(sy, sy) * R); b = (bx @ Oi @ by) / (bx @ Oi @ bx); se = np.sqrt(1 / (bx @ Oi @ bx))
    r = by - b * bx; Q = r @ Oi @ r; phi = max(1, Q / (len(bx) - 1)); return b, se * np.sqrt(phi)
def mvivw(X, by, sy, R):
    Oi = np.linalg.pinv(np.outer(sy, sy) * R); A = X.T @ Oi @ X; b = np.linalg.solve(A, X.T @ Oi @ by); V = np.linalg.inv(A)
    r = by - X @ b; Q = r @ Oi @ r; phi = max(1, Q / (len(by) - X.shape[1])); return b, np.sqrt(np.diag(V) * phi)
fmt = lambda b, se: f"OR {np.exp(-b):.2f} ({np.exp(-b-1.96*se):.2f}-{np.exp(-b+1.96*se):.2f}) P {2*stats.norm.sf(abs(b/se)):.2g}"
snps, LD = ldmat(k.rsid.tolist(), k.low.tolist()); k = k.set_index("rsid").loc[snps].reset_index()
def wm(bx, by, sy):
    est = by / bx; w = (bx / sy) ** 2; o = np.argsort(est); est, w = est[o], w[o] / w.sum(); c = np.cumsum(w) - 0.5 * w
    b = np.interp(0.5, c, est)
    bs = []; rng = np.random.default_rng(1)
    for _ in range(1000):
        bb = rng.normal(bx, k.sx.values) ; yy = rng.normal(by, sy); e = yy / bb; o2 = np.argsort(e); ww = ((bb/sy)**2)[o2]; ww /= ww.sum(); cc = np.cumsum(ww) - 0.5 * ww; bs.append(np.interp(0.5, cc, e[o2]))
    return b, np.std(bs)
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA"]:
    o = outcome(ep); m = k.merge(o, on="pos")
    sg = np.where((m.alt == m.low) & (m.ref == m.oth), 1, np.where((m.ref == m.low) & (m.alt == m.oth), -1, 0)); m = m[sg != 0]; by = m.beta.values * sg[sg != 0]
    idx = [snps.index(s) for s in m.rsid]; R = LD[np.ix_(idx, idx)]
    b, se = givw(m.bx.values, by, m.sebeta.values, R)
    keep = (m.kp >= 1e-4).values
    b2, se2 = givw(m.bx.values[keep], by[keep], m.sebeta.values[keep], R[np.ix_(keep, keep)])
    bm, sem = mvivw(np.column_stack([m.bx.values, m.bk.values]), by, m.sebeta.values, R)
    print(f"== {ep}\n  main IVW          {fmt(b, se)}\n  excl KLKB1-assoc  {fmt(b2, se2)} (n={keep.sum()})\n  MVMR FXI|KLKB1    {fmt(bm[0], sem[0])}")
# weighted median and Egger on near-independent set (r2<0.01)
iv2 = pd.read_csv("outputs/clump_r0.01.clumped", sep=r"\s+").SNP
j = ins[ins.rsid.isin(iv2)]
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA"]:
    o = outcome(ep); m = j.merge(o, on="pos"); sg = np.where((m.alt == m.low) & (m.ref == m.oth), 1, np.where((m.ref == m.low) & (m.alt == m.oth), -1, 0)); m = m[sg != 0]; by = m.beta.values * sg[sg != 0]
    bx, sy = m.bx.values, m.sebeta.values
    est = by / bx; w = (bx / sy) ** 2; o_ = np.argsort(est); c = np.cumsum(w[o_] / w.sum()) - 0.5 * w[o_] / w.sum(); bwm = np.interp(0.5, c, est[o_])
    # Egger (orient bx positive per convention: flip to FXI-raising)
    X = np.column_stack([np.ones(len(bx)), -bx]); W = np.diag(1 / sy**2); beg = np.linalg.solve(X.T @ W @ X, X.T @ W @ (-by)); res = -by - X @ beg
    s2 = max(1, (res @ W @ res) / (len(bx) - 2)); seg = np.sqrt(np.diag(np.linalg.inv(X.T @ W @ X)) * s2)
    print(f"  {ep:28s} WM OR {np.exp(-bwm):.2f} | Egger slope OR {np.exp(-beg[1]):.2f} ({np.exp(-beg[1]-1.96*seg[1]):.2f}-{np.exp(-beg[1]+1.96*seg[1]):.2f}) intercept P {2*stats.norm.sf(abs(beg[0]/seg[0])):.2g}")
# replication exposure (SomaScan): primary instruments with replication pQTL effects
r = ins.merge(rep[["pos", "ea", "oa", "b", "s", "p"]].rename(columns={"ea": "rea", "oa": "roa", "b": "rb", "s": "rs_", "p": "rp"}), on="pos")
r["rbx"] = np.where(r.rea == r.low, r.rb, np.where(r.roa == r.low, -r.rb, np.nan)); r = r.dropna(subset=["rbx"])
print(f"== replication pQTL: {len(r)} instruments; sign concordance {(np.sign(r.rbx)==np.sign(r.bx)).mean():.2f}; corr {np.corrcoef(r.rbx, r.bx)[0,1]:.2f}; replication-P<5e-8 instruments {(r.rp<5e-8).sum()}")
top = r.loc[r.rp.idxmin()]
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA"]:
    o = outcome(ep).set_index("pos").loc[top.pos]; sg = 1 if (o.alt == top.low) else -1; by = o.beta * sg
    w = by / top.rbx; se = o.sebeta / abs(top.rbx); print(f"  rep Wald {top.rsid} {ep:28s} {fmt(w, se)} (per SD lower FXI, SomaScan scale)")
