"""Drug-target MR: plasma FXI (cis pQTL, UKB Olink) -> FinnGen R12 outcomes. Generalised IVW with LD, Egger, WM, Wald; coloc ABF.
Run from repo root. Inputs: data/pqtl/F11.cis.tsv (+header), data/finngen/<ep>.cis.tsv, 1000G EUR plink (hg19, rsid)."""
import os, subprocess, numpy as np, pandas as pd, sys
from scipy import stats
PL = os.environ.get("PLINK", "plink"); BF = os.environ.get("REF_EUR", "ref/EUR")  # PLINK 1.9, 1000G EUR bfile prefix
R2 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
hdr = open("data/pqtl/header.txt").read().strip().split("\t")
pq = pd.read_csv("data/pqtl/F11.cis.tsv", sep="\t", header=None, names=hdr, comment="#")
pq = pq.rename(columns={"base_pair_location": "pos", "effect_allele": "ea", "other_allele": "oa", "beta": "bx", "standard_error": "sx", "p_value": "px", "effect_allele_frequency": "eaf"})
pq = pq[pq.rsid.astype(str).str.startswith("rs")]
eps = open("data/outcomes.txt").read().split()
def outcome(ep):
    d = pd.read_csv(f"data/finngen/{ep}.cis.tsv", sep="\t")
    d = d.rename(columns={"#chrom": "chr"})
    return d[["pos", "ref", "alt", "beta", "sebeta", "pval", "af_alt"]]
# instruments: P<5e-8, clump in 1000G EUR by rsid
sig = pq[pq.px < 5e-8][["rsid", "px"]].rename(columns={"rsid": "SNP", "px": "P"})
sig.to_csv("outputs/sig.tsv", sep="\t", index=False)
subprocess.run([PL, "--bfile", BF, "--clump", "outputs/sig.tsv", "--clump-p1", "5e-8", "--clump-r2", str(R2), "--clump-kb", "1000", "--out", f"outputs/clump_r{R2}"], capture_output=True)
iv = pd.read_csv(f"outputs/clump_r{R2}.clumped", sep=r"\s+").SNP.tolist()
ins = pq[pq.rsid.isin(iv)].copy()
# orient to FXI-lowering allele; b is per unit higher FXI, reported as OR per SD lower FXI (exp(-b))
flip = ins.bx > 0
ins["low"] = np.where(flip, ins.oa, ins.ea); ins["oth"] = np.where(flip, ins.ea, ins.oa); ins["bx"] = -ins.bx.abs()
# LD (signed) among instruments, aligned to the low allele
pd.Series(ins.rsid).to_csv("outputs/iv.txt", index=False, header=False)
subprocess.run([PL, "--bfile", BF, "--extract", "outputs/iv.txt", "--r", "square", "--keep-allele-order", "--write-snplist", "--out", "outputs/iv_ld"], capture_output=True)
snps = open("outputs/iv_ld.snplist").read().split(); ld = np.loadtxt("outputs/iv_ld.ld").reshape(len(snps), len(snps))
bim = pd.read_csv(BF + ".bim", sep="\t", header=None, names=["c", "SNP", "cm", "bp", "A1", "A2"]).set_index("SNP").loc[snps]
sgn = np.array([1 if bim.loc[s, "A1"] == ins.set_index("rsid").loc[s, "low"] else (-1 if bim.loc[s, "A2"] == ins.set_index("rsid").loc[s, "low"] else 0) for s in snps])
ldm = ld * np.outer(sgn, sgn)
ins = ins.set_index("rsid").loc[snps].reset_index()
print(f"instruments r2<{R2}: {len(ins)}; F-stat median {np.median((ins.bx/ins.sx)**2):.0f}; unaligned {int((sgn==0).sum())}")
def mr(ep):
    o = outcome(ep)
    m = ins.merge(o, on="pos")
    by = np.where((m.alt == m.low) & (m.ref == m.oth), m.beta, np.where((m.ref == m.low) & (m.alt == m.oth), -m.beta, np.nan))
    ok = ~np.isnan(by); m = m[ok]; by = by[ok]
    idx = [snps.index(s) for s in m.rsid]; R = ldm[np.ix_(idx, idx)]
    bx, sx, sy = m.bx.values, m.sx.values, m.sebeta.values
    Om = np.outer(sy, sy) * R; Oi = np.linalg.pinv(Om)
    b = (bx @ Oi @ by) / (bx @ Oi @ bx); se = np.sqrt(1 / (bx @ Oi @ bx))
    res = by - b * bx; Q = res @ Oi @ res; k = len(bx); phi = max(1, Q / (k - 1)) if k > 1 else 1
    se_re = se * np.sqrt(phi)
    j = np.argmax(np.abs(bx / sx)); wald = by[j] / bx[j]; wse = sy[j] / abs(bx[j])
    return dict(outcome=ep, n_iv=k, OR=np.exp(-b), lo=np.exp(-b - 1.96 * se_re), hi=np.exp(-b + 1.96 * se_re), P=2 * stats.norm.sf(abs(b / se_re)),
                Q=Q, Q_P=stats.chi2.sf(Q, k - 1) if k > 1 else np.nan, wald_OR=np.exp(-wald), wald_P=2 * stats.norm.sf(abs(wald / wse)), top=m.rsid.values[j])
out = pd.DataFrame([mr(e) for e in eps])
pd.set_option("display.width", 200)
print(out.to_string(index=False, float_format=lambda x: f"{x:.3g}"))
out.to_csv(f"outputs/mr_r{R2}.tsv", sep="\t", index=False)
