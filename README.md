# MSC-OriginPairs

Code and classifier for the paper:

> Zhen Zhang, Tianzhu Xie, Chunqi Xu, Sisi Deng, Yubin Sun, Xi Wu, Zeyu Huang,
> Ning Chen, Meirong Qin and Ping Wang. Multi-study transcriptomic analysis
> identifies reproducible tissue-of-origin signatures in human mesenchymal
> stromal cells. *STEM CELLS* (submitted).

Archived release: [Zenodo DOI]

MSC-OriginPairs assigns the tissue of origin of culture-expanded human mesenchymal
stromal cells (bone marrow, adipose tissue or umbilical cord) from 30 within-sample
gene-pair rules. Each sample is classified on its own, without normalization against
other samples. This repository also contains the analysis pipeline used in the paper.

## Repository contents

| Folder | Contents |
|---|---|
| `classifier/` | Final model, rule table, prediction script and worked example |
| `analysis/` | R analysis pipeline, numbered in run order |
| `R/` | Shared R functions used by the pipeline |
| `figures/` | Python scripts that draw the figures from result tables |
| `config/` | Analysis settings (thresholds, seeds, class definitions) |
| `data/` | Dataset list, download plan, result tables and codebook |
| `curation/` | Curated input lists (reported markers, HOX genes, study lists) |
| `tests/` | Unit tests on simulated data |
| `environment/` | Full R session information |

## System requirements

- Tested on Windows 11 (x64) with R 4.5.2 and Python 3.12. The code uses no
  OS-specific features and should run on Linux and macOS.
- Classifier: R >= 4.1 with `multiclassPairs` (0.4.3), `switchBox` (1.44.0) and `readr`.
- Full pipeline: the R packages listed in `environment/sessionInfo_R.txt`, including
  limma 3.64.3, edgeR 4.6.3, glmnet 4.1.10, lme4 1.1.38 and variancePartition 1.38.1.
- Figures: Python packages in `requirements.txt`.
- No special hardware is required; the full pipeline was run with 16 GB RAM.

## Installation

```r
install.packages(c("readr", "multiclassPairs"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("switchBox")
```

Then clone this repository. Typical installation time: about 3 minutes on a
standard desktop computer.

## Demo

```
Rscript classifier/predict_msc_origin.R classifier/example_input.tsv
```

The example contains three samples from the external validation cohort, one per
source. The expected output is `classifier/example_expected_output.tsv`: each sample is
assigned to its labeled source, and the umbilical cord sample is flagged as a
low-confidence call. Expected run time: under 10 seconds.

## Using the classifier on your own data

- **Input:** a tab-separated file. The first column, named `gene`, holds HGNC symbols;
  each further column holds one sample.
- **Scale:** any scale on which higher values mean higher expression and the
  within-sample ordering is preserved: log-scale expression or within-sample ranks.
  For qPCR, use negative Ct values (−Ct), because lower Ct means higher expression.
- **Missing genes:** rules involving a missing gene are skipped; samples are scored
  independently.
- **Output:** the call (AD, BM or UC), three source scores, the confidence (the margin
  between the two highest scores) and a low-confidence flag (threshold 0.2494).
- **Scope:** the classifier is designed to verify samples labeled as BM, AD or UC. It
  always returns one of these three sources. Cells from other tissues (for example
  bone or dental pulp) may be called BM, and admixture below about 20% is not reliably
  detected. Low-confidence calls should be checked against other information.

## Reproducing the analyses

1. Download the public expression data listed in `data/download_plan_v2.tsv`
   (`analysis/10_download_verify.R` checks that downloads are complete).
2. Run the scripts in `analysis/` in numbered order from `10_prepare_plan.R`. Scripts
   `01`–`05` document the dataset census and are provided for reference only; their
   curated output is `data/dataset_uses_v2.tsv` and Supplementary Table S1.
3. Run the scripts in `figures/` to draw all figures.
4. Run `Rscript tests/run_tests.R` to check the shared functions on simulated data.

Expected run time for the full pipeline: several hours to about a day on a standard
desktop computer (most of it is downloading and assembling the public expression data).
The result tables underlying the figures are included in this repository (`data/`).

Study groups are labeled with the codes G01–G92 used in Supplementary Table S1.

## License

Code is released under the MIT License. Result tables in `data/` are released
under CC BY 4.0.

## Citation

Please cite the paper above and the archived release ([Zenodo DOI]). Citation metadata
is in `CITATION.cff`.

## Contact

Zhen Zhang, zhangzhen7812@live.com
