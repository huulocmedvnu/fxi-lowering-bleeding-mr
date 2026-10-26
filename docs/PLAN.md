# PLAN: genetically proxied factor XI lowering and bleeding in women (drug-target Mendelian randomization)

Registered 2026-10-26 before any outcome association in the F11 region was extracted for this project.
Disclosure: in the earlier project vte_pph_gen (submitted, JTHA-S-26-02148) two F11 variants associated with
venous thromboembolism (VTE) were tested against postpartum haemorrhage (PPH) in FinnGen (P 0.07 and 0.39,
risk allele direction towards lower PPH). No pQTL-based analysis of F11 was done.

## Question
Factor XI (FXI) inhibitors (abelacimab, milvexian, asundexian) are in phase 3 trials as anticoagulants with less
bleeding. Women of reproductive age on anticoagulants have heavy menstrual bleeding, and pregnant women are
excluded from trials. Does lifelong genetically lower plasma FXI raise the risk of postpartum haemorrhage or heavy
menstrual bleeding, while lowering VTE?

## Exposure
Plasma F11 protein (Olink, UK Biobank, n 47,745; GCST90469165, hg38). Instruments: variants within 500 kb of F11
(hg38 chr4:186,265,524-186,289,681) with P < 5e-8, clumped at r2 < 0.1 (primary) and r2 < 0.01 (sensitivity) in
1000 Genomes EUR. Effects oriented to the allele lowering FXI; estimates per SD lower FXI.
Independent exposure data for replication: GCST90424248 (n 3,506) if instruments are available.

## Outcomes (FinnGen R12, no overlap with UK Biobank)
- Positive control: VTE (I9_VTE). Expected: lower risk with lower FXI.
- Primary safety: PPH (O15_POSTPART_HEAMORRH_ALLW), heavy menstrual bleeding (N14_MENORRHAGIA).
- Secondary: PPH due to atony and due to retained placenta, antepartum haemorrhage, gastrointestinal bleeding,
  intracerebral haemorrhage, epistaxis, iron deficiency anaemia (where available in R12 with >= 1,000 cases).
- Negative control outcome: gestational hypertension.

## Analysis
Inverse-variance weighted MR accounting for residual correlation between instruments (generalised IVW with LD
matrix from 1000G EUR); MR-Egger and weighted median as sensitivity; single strongest cis variant (Wald ratio).
Colocalisation (coloc, default priors) of FXI pQTL with each outcome in the cis window; also with plasma
prekallikrein (KLKB1, adjacent gene) to flag shared signals.
Positive-control requirement: the analysis is interpreted only if FXI lowering associates with lower VTE (P < 0.05).
Safety outcomes: report odds ratio per SD lower FXI with 95% CI; an increase in PPH or heavy menstrual bleeding with
P < 0.05/2 is called a signal; a null result is reported with the CI bound to state the largest increase compatible
with the data.

## Computation
cis-window extraction of summary statistics in Rust (or tabix); timed on one file first.

## Amendments
### 2026-10-26a (after primary and follow-up results; before replication outcome data were extracted)
Target journal ATVB requires replication of MR findings in an independent sample. Outcome replication:
UK Biobank (Pan-UKBB, EUR, GRCh37): O72 postpartum haemorrhage, N92 excessive/frequent/irregular menstruation;
VTE positive control: I26 pulmonary embolism, I80 phlebitis and thrombophlebitis; and MVP VTE (GCST90476007, no UKB).
Same instruments (r2<0.1), matched by rsID through 1000G EUR (GRCh37). The FXI pQTL sample (47,745) is a subset of
UK Biobank, so UKB outcome estimates have partial sample overlap (about 10%); this is stated, and the SomaScan
exposure (non-UKB) Wald estimate is also given for UKB outcomes. FinnGen and UKB estimates combined by fixed-effect
inverse-variance meta-analysis. Endpoint deviations from the original list (GI bleeding and epistaxis not available
in FinnGen R12; composite BLEEDING and anaemia endpoints used) are recorded here.
