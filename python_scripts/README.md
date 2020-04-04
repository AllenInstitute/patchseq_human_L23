# Python scripts
Code for reproducing the analyses presented in "Human cortical expansion involves diversification and specialization of  supragranular intratelencephalic-projecting neurons".

## Prerequisites

Install the necessary packages in a Python 3.7 environment via `pip install -r requirements.txt`

## Scripts

All high-level analysis and figure code is presented in jupyter notebooks as follows:

`pathology_analysis.ipynb`: Pathology and electrophysiology analysis from Figure 2.

`morpho_electro_feature_analysis.ipynb`: All remaining depth, electrophysiology, and morphology analysis from main paper, Figures 4-6.

`classifier_and_species_analysis.ipynb`: Additional Extended Data results including random forest classifier for t-types and cross-species comparison of electrophysiology.

Additionally, code to calculate gene-depth correlations (very compute-intensive) is in `gene_depth_corr_calculation.py`. For convenience, results are included in data so there is no need to re-run this analysis.
