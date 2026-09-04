# PRINCE Mito3 survival sensitivity analysis

The release includes a post hoc exploratory OS/PFS sensitivity analysis for the
38 nivolumab-exposed patients used in the response association. The clinical
fields are `clinical.observation.os` and `clinical.observation.pfs` with their
corresponding boolean event fields. Mito3 is entered per one sample standard
deviation in treatment-arm-adjusted Cox models. Kaplan–Meier curves use a
cohort-median split only for visualization; the cutoff was not optimized.

The aggregate results are in `results/PRINCE_Mito3_survival_summary.tsv` and
the curves are in `figures/Supplementary_FigureS11_PRINCE_Mito3_survival.pdf`.
Run the script with a locally obtained, permitted Mito3 score table:

```text
python scripts/analyze_prince_mito3_survival.py --scores <local-score-table.tsv> --output results/PRINCE_Mito3_survival_summary.tsv --figure figures/Supplementary_FigureS11_PRINCE_Mito3_survival.pdf
```

The survival analysis is outside the seven-module response-screen BH family and
does not establish prognostic value, predictive treatment interaction, or a
clinical cutoff.
