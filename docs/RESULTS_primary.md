# Primary results 2026-10-26 (OR per SD lower plasma FXI; generalised IVW, multiplicative random effects)
Instruments: cis F11 pQTL (UKB Olink n 47,745) P<5e-8; r2<0.1 primary (65 clumped, 57 in FinnGen, median F 59); r2<0.01 sensitivity (18/16).
Outcome endpoints chosen from R12 by availability: GI bleeding and epistaxis not available; BLEEDING composite, iron deficiency
anaemia due to blood loss, acute posthaemorrhagic anaemia and other abnormal uterine bleeding used instead (deviation logged).
Positive control VTE: 0.54 (0.48-0.60), P 7e-25 (sensitivity 0.53). PASS. Heterogeneity high (Q P 5e-31).
PPH: 1.03 (0.95-1.11), P 0.50; sensitivity 1.09 (0.98-1.21), P 0.12; lead variant rs3756008 Wald 1.09, P 0.08.
Heavy menstrual bleeding: 0.96 (0.90-1.03), P 0.23; sensitivity 0.97 (0.91-1.04).
PPH atony 1.03 (0.90-1.17); retained placenta 1.01 (0.88-1.15); antepartum 0.96 (0.83-1.10); any bleeding 1.00 (0.97-1.05);
ICH 1.07 (0.95-1.20); IDA blood loss 0.92 (0.83-1.03); other abnormal uterine bleeding 1.00 (0.93-1.09);
acute posthaemorrhagic anaemia 1.11 (0.86-1.43); gestational hypertension (negative control) 1.00 (0.92-1.07).
No bleeding signal at P<0.025. Largest increase compatible (upper CI): PPH +11% (primary), +21% (sensitivity); HMB +3%.

# Follow-up 2026-10-26 (scripts/followup.py, outputs/followup.txt)
Coloc (p12 1e-5): FXI x VTE PP4 0.999; FXI x PPH PP1 0.66, PP4 0.025; FXI x HMB PP1 0.80, PP4 0.005;
  FXI x gestational hypertension PP4 0.007; FXI x KLKB1 PP3 1.000 (distinct causal variants).
KLKB1: 46/65 instruments have KLKB1 P<1e-4 (regional LD) but effect correlation 0.04. Excluding them (n=17):
  VTE 0.54 (0.45-0.65), PPH 0.93 (0.81-1.07), HMB 0.92 (0.84-1.02). MVMR FXI|KLKB1: VTE 0.53, PPH 1.03 (0.95-1.11), HMB 0.96 (0.90-1.02).
Weighted median / MR-Egger (r2<0.01): VTE 0.52 / 0.49 (0.44-0.55); PPH 1.10 / 1.16 (0.98-1.37), intercept P 0.40;
  HMB 0.96 / 0.97 (0.87-1.09).
Replication exposure (SomaScan, n 3,506, GCST90424248): 43 shared instruments, sign concordance 0.81, r 0.42,
  none P<5e-8; Wald at rs3775297: VTE 0.47 (0.38-0.57), PPH 0.94 (0.70-1.28), HMB 0.92 (0.74-1.14).
Pregnancy-associated VTE: FinnGen endpoints O15_PREG_VENOUS_HAEMORRH (685) and O15_VENOUS_COMPLIC (389) include
  haemorrhoids and are below 1,000 cases; not analysed (limitation).
