# PRISM mitochondrial-state-selective drug screen

The screen matched 31 PDAC cell models and 4686 compounds. Cell lines were split at the median tumor `mitochondrial_oxidative` expression score; lower PRISM log-fold-change indicates stronger growth inhibition.

## Result

0 compounds met the exploratory ranking rule (low-state minus high-state median LFC < -0.25 and FDR < 0.20). The top candidates and all prespecified-pathway-related compounds are reported in machine-readable tables.

No compound from this screen can be called RICTOR-specific, ICI-synergistic, or clinically effective. The analysis is a cell-line state-association screen with a median split, not a randomized treatment comparison. FDR is reported across the complete tested compound set.

## Top-ranked rows

    compound                                          moa                                                                                                                                                    target  delta_lfc_low_minus_high  fdr
 digoxigenin                                      steroid                                                                                                                                                       NaN                 -1.341622  1.0
  mevastatin                              HMGCR inhibitor                                                                                                                                                     HMGCR                 -1.321054  1.0
     AT-7519                                CDK inhibitor                                                                                                                        CDK1, CDK2, CDK4, CDK5, CDK6, CDK9                 -1.237698  1.0
pevonedistat             nedd activating enzyme inhibitor                                                                                                                                                NAE1, UBA3                 -1.103843  1.0
     AZD5438                                CDK inhibitor                                                                                                                                                     KCNH2                 -1.018172  1.0
       M-344                               HDAC inhibitor                                                                                                                                                       NaN                 -0.963213  1.0
 doxorubicin                      topoisomerase inhibitor                                                                                                                                                     TOP2A                 -0.959835  1.0
  bortezomib NFkB pathway inhibitor, proteasome inhibitor PSMA1, PSMA2, PSMA3, PSMA4, PSMA5, PSMA6, PSMA7, PSMA8, PSMB1, PSMB10, PSMB11, PSMB2, PSMB3, PSMB4, PSMB5, PSMB6, PSMB7, PSMB8, PSMB9, PSMD1, PSMD2, RELA                 -0.930245  1.0
        AZ20                         ATR kinase inhibitor                                                                                                                                                 ATR, MTOR                 -0.913149  1.0
  hesperadin                      Aurora kinase inhibitor                                                                                                                                                     AURKB                 -0.901141  1.0

## Outputs

- `all_compounds_state_selectivity.tsv`
- `candidate_state_selective_compounds.tsv`
- `prespecified_pathway_related_compounds.tsv`