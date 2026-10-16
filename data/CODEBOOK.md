# Codebook

## Adjudication codes

These codes appear in the `human_verdict` / `verdict` columns of the census tables
(`data/census_adjudication_v2_r3.tsv`, `data/dataset_uses_v2.tsv`) and in
Supplementary Table S1.

| Code | Meaning |
|---|---|
| C_AD / C_BM / C_UC | external-validation candidate; dominant source adipose tissue / bone marrow / umbilical cord |
| D | discovery-cohort study (development layer) |
| MULTI_NEW | additional multi-source study (more than one source measured within the study) |
| ORTHO_SCRNA | orthopaedic score / single-cell study kept for the scRNA arm |
| PASSAGE_SERIES | passage-series study (passage-stability analysis) |
| X_NOT_MSC | not an MSC / stromal-cell culture profile |
| X_NOT_CULTURE | primary tissue or dissociated cells, not culture-expanded MSCs |
| X_NOT_HUMAN | non-human material |
| X_NOT_EXPR | no expression matrix |
| X_CELLLINE_IPSC | cell line or iPSC-derived material |
| X_TREAT_NO_CTRL | treatment design without a usable control arm |
| X_SOURCE_OOC | other tissue of origin (out-of-class; application scope only) |
| X_FETAL_EMBRYO | fetal or embryonic material |
| X_NO_MATRIX | matrix not available |
| X_SMALL_N | below the required minimum sample count |
| X_DUP | duplicate series |
| X_EXPLORATORY_PHASE | carried over from the earlier exploratory phase of the project |
| UNCLEAR | insufficient information; not carried into any analysis |

## data/dataset_uses_v2.tsv — column meanings

| Column | Meaning |
|---|---|
| gse | GEO series accession |
| uses_en | cohorts and analyses that use the study (baseline scope and application scope) |
| n_samples_baseline | samples entering the baseline (primary) analyses |
| n_samples_application | samples entering the application-scope analyses (non-baseline arms, out-of-class and mixture analyses) |
| n_samples_total | total number of samples used |
| cohort_membership | cohort memberships, combined with `+` when a study belongs to more than one (D, C, P, OOC, ARM) |
| verdict | final adjudication code (see the table above) |
| note | free-text curation remark |

## Study-group codes

Study groups are labelled with the neutral codes G01-G92 used throughout the
repository and in Supplementary Table S1. The accession-to-group table is
`curation/team_numbering.tsv`.
