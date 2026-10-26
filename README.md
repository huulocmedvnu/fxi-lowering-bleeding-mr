# Factor XI lowering and bleeding in women: drug-target Mendelian randomization

Analysis code for "Genetically Proxied Factor XI Lowering and Bleeding in Women:
A Drug-Target Mendelian Randomization Study" (T.H.L. Tran).

## Analysis plan
`docs/PLAN.md` was committed (288f467) before any outcome association in the
F11 region was extracted. The amendment adding outcome replication in the UK
Biobank and the Million Veteran Program was committed in 4c09b42 before those
data were extracted. The git history records the order of these steps.

## Data (all public summary statistics)
- Exposure: plasma factor XI and prekallikrein, UK Biobank Olink
  (GWAS Catalog GCST90469165, GCST90469708), harmonised files.
- Outcomes: FinnGen release 12 endpoints in `data/outcomes.txt`,
  Pan-UK Biobank ICD-10 O72, N92, I26, I80 (both sexes, EUR),
  Million Veteran Program VTE (GCST90476007).
- The cis window (F11 +/- 500 kb, GRCh38) was read from the remote
  bgzipped files with `tabix` into `data/<source>/<endpoint>.cis.tsv`.
- LD reference: 1000 Genomes Phase 3 EUR in PLINK format.

## Run
Set `PLINK` (PLINK 1.9 binary) and `REF_EUR` (1000G EUR bfile prefix), then
from the repo root:

    python3 scripts/mr_fxi.py 0.1      # primary, r2 < 0.1
    python3 scripts/mr_fxi.py 0.01     # sensitivity, r2 < 0.01
    python3 scripts/followup.py        # coloc, KLKB1, weighted median, Egger
    python3 scripts/replicate.py       # UK Biobank and MVP replication
    python3 scripts/final_tables.py    # outputs/results_all.tsv
    python3 scripts/figures.py         # Figures 2-4

Python 3.9 with numpy, pandas, scipy, plotly, kaleido, Pillow.
`outputs/results_all.tsv` holds every estimate reported in the manuscript and
`outputs/TableS1_instruments.tsv` lists the instruments.

## License
MIT
