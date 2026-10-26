"""Outcome replication (amendment 2026-10-26a): UKB (Pan-UKBB EUR) and MVP VTE; FinnGen + UKB fixed-effect meta.
OR per SD lower FXI, generalised IVW with LD (r2<0.1 instruments), multiplicative random effects."""
import subprocess, numpy as np, pandas as pd
from scipy import stats
PL = "/mnt/hdd2/loc-tran/miniconda3/envs/ldsc3/bin/plink"; BF = "/mnt/hdd2/loc-tran/birth_canal_gen/ref/EUR"
hdr = open("data/pqtl/header.txt").read().strip().split("\t")
fx = pd.read_csv("data/pqtl/F11.cis.tsv", sep="\t", header=None, names=hdr).rename(columns={"base_pair_location": "pos", "effect_allele": "ea", "other_allele": "oa", "beta": "b", "standard_error": "s"})
iv = pd.read_csv("outputs/clump_r0.1.clumped", sep=r"\s+").SNP
ins = fx[fx.rsid.isin(iv)].drop_duplicates("rsid").copy(); fl = ins.b > 0
ins["low"] = np.where(fl, ins.oa, ins.ea); ins["oth"] = np.where(fl, ins.ea, ins.oa); ins["bx"] = -ins.b.abs()
bim = pd.read_csv(BF + ".bim", sep="\t", header=None, names=["c", "rsid", "cm", "pos37", "A1", "A2"])
ins = ins.merge(bim[["rsid", "pos37"]], on="rsid", how="left")
pd.Series(ins.rsid).to_csv("outputs/tmp_iv.txt", index=False, header=False)
subprocess.run([PL, "--bfile", BF, "--extract", "outputs/tmp_iv.txt", "--r", "square", "--keep-allele-order", "--write-snplist", "--out", "outputs/tmp_ld"], capture_output=True)
snps = open("outputs/tmp_ld.snplist").read().split(); L = np.loadtxt("outputs/tmp_ld.ld").reshape(len(snps), len(snps))
b2 = bim.set_index("rsid").loc[snps]; lo = dict(zip(ins.rsid, ins.low)); g = np.array([1 if b2.loc[x, "A1"] == lo[x] else -1 for x in snps]); LD = L * np.outer(g, g)
def givw(m, by, sy):
    idx = [snps.index(s) for s in m.rsid]; R = LD[np.ix_(idx, idx)]; bx = m.bx.values
    Oi = np.linalg.pinv(np.outer(sy, sy) * R); b = (bx @ Oi @ by) / (bx @ Oi @ bx); se = np.sqrt(1 / (bx @ Oi @ bx))
    r = by - b * bx; Q = r @ Oi @ r; return b, se * np.sqrt(max(1, Q / (len(bx) - 1))), len(bx)
def align(m, ref, alt, beta):
    s = np.where((m[alt] == m.low) & (m[ref] == m.oth), 1, np.where((m[ref] == m.low) & (m[alt] == m.oth), -1, 0)); k = s != 0
    return m[k], (m[beta].values * s)[k]
res = {}
def report(name, b, se, n):
    res[name] = (b, se); print(f"{name:34s} n_iv={n:3d} OR {np.exp(-b):.2f} ({np.exp(-b-1.96*se):.2f}-{np.exp(-b+1.96*se):.2f}) P {2*stats.norm.sf(abs(b/se)):.2g}")
for ep in ["I9_VTE", "O15_POSTPART_HEAMORRH_ALLW", "N14_MENORRHAGIA"]:
    o = pd.read_csv(f"data/finngen/{ep}.cis.tsv", sep="\t").rename(columns={"#chrom": "chr"}).drop_duplicates("pos")
    m, by = align(ins.merge(o, on="pos"), "ref", "alt", "beta"); report("FinnGen " + ep, *givw(m, by, m.sebeta.values))
for f, lab in [("icd10-O72-both_sexes", "UKB O72 PPH"), ("icd10-N92-both_sexes", "UKB N92 heavy/irregular menses"), ("icd10-I26-both_sexes", "UKB I26 pulmonary embolism"), ("icd10-I80-both_sexes", "UKB I80 phlebitis/thrombophlebitis")]:
    h = open(f"data/ukb/{f}.header").read().strip().split("\t")
    o = pd.read_csv(f"data/ukb/{f}.cis.tsv", sep="\t", header=None, names=h)
    o = o[(o.low_confidence_EUR.astype(str) == "false") | (o.low_confidence_EUR.astype(str) == "False")].rename(columns={"pos": "pos37"}).drop_duplicates("pos37")
    m, by = align(ins.merge(o[["pos37", "ref", "alt", "beta_EUR", "se_EUR"]], on="pos37"), "ref", "alt", "beta_EUR"); report(lab, *givw(m, by, m.se_EUR.values))
o = pd.read_csv("data/mvp/VTE.cis.tsv", sep="\t", na_values=["#NA"]).rename(columns={"base_pair_location": "pos"})
o["beta"] = np.log(o.odds_ratio); o["se"] = (np.log(o.ci_upper) - np.log(o.ci_lower)) / 3.92
o = o.dropna(subset=["beta", "se"]).drop_duplicates("pos")
m, by = align(ins.merge(o[["pos", "effect_allele", "other_allele", "beta", "se"]], on="pos"), "other_allele", "effect_allele", "beta"); report("MVP VTE", *givw(m, by, m.se.values))
def meta(a, b, lab):
    (b1, s1), (b2, s2) = res[a], res[b]; w1, w2 = 1 / s1**2, 1 / s2**2; bm = (b1 * w1 + b2 * w2) / (w1 + w2); sm = np.sqrt(1 / (w1 + w2))
    q = (b1 - bm)**2 * w1 + (b2 - bm)**2 * w2; report(lab, bm, sm, 0); print(f"{'':34s} heterogeneity P {stats.chi2.sf(q,1):.2g}")
print("== fixed-effect meta FinnGen + UKB")
meta("FinnGen O15_POSTPART_HEAMORRH_ALLW", "UKB O72 PPH", "META PPH")
meta("FinnGen N14_MENORRHAGIA", "UKB N92 heavy/irregular menses", "META heavy menstrual bleeding")
