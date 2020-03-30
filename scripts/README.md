# Overview
This folder contains the analysis scripts.  Prior to running these scripts, download all of them into your **main working directory** (*not* a scripts subdirector).  Then run the scripts in order and follow any additional directions indicated below or within the scripts.  

# List of scripts

This document is divided into two sections based on which programming language was used for the scripts.  The R scripts should be run prior to running the python scripts, as the python scripts required cluster calls for the human Patch-seq cells (which are geenrated using R).  If only the Python scripts are to be run, we provide cluster calls for Patch-seq cells in the `data` directory.

## R scripts

To replicate the figures in the paper, please use the code below:  

1. **Code 1. Download and prepare the dissociation data** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/Code1_prepareComparisonDataSets.nb.html)  This script converts data downloaded from the [Allen Institute Cell Types Database](http://celltypes.brain-map.org/rnaseq) into a format compatible for use as comparison data to human and mouse Patch-seq.  
2. **Code 2. Complete preparation for dissociation and Patch-seq data** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/Code2_prepare_data_L23exc.nb.html).  This script prepares the data for comparing mouse VISp and human MTG using dissociated and Patch-seq in superficial cortical layers, and outputs the results as R objects for remaining scripts.  
3. **Code 3. Re-analysis and visualization of FACS data** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/Code3_FREM3_crossSpeciesComparisons.nb.html).  This script performs some statistical analyses and makes some UMAP plots showing glutamatergic neurons in supragranular layers of human MTG and mouse VISp for Figure 1.  One purpose is to quantify how cell type diversity, variation, and specialization of human types compares with mouse.  We also plot gene expression of some genes for Figure 6.  
4. **Code 4. Visualization of mouse Patch-seq data** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/*****CORRECT_LINK*******.nb.html).  This script compares mouse FACs to Patch-seq data and visualizes the results for Figure 3.  It also outputs some plots and tables for later figures.  
5. **Code 5. Mapping and visualization of human Patch-seq data** [(LINK TO SCRIPT)](http://htmlpreview.github.io/?https://github.com/AllenInstitute/patchseq_human_L23/blob/master/scripts/*****CORRECT_LINK*******.nb.html).  This script first performs Seurat mapping on the human Patch-seq cells (using FACS MTG data as references) and then compares and visualizes these results for Figure 3.  It also outputs some plots and tables (including mapping results) for later figures.  
6. **Code 6&7. Soma size/density quantification**: `Code6_cortical_depth_code_human.r` and `Code7_cortical_depth_code_mouse.r` together perform quantifications of neuron density and soma size in mouse and human that are presented in Figure 1a-c.  **Before running these two scripts** data in the `data\cell_sizes` folder needs to be downloaded and unzipped in a specific way:
* Download everything in `data\cell_sizes` and transfer to your working directory
* Unzip `mouse_cell_sizes.zip` (make sure the directory is in the format `[working_directory]/mouse_cell_sizes/case_set_1` rather than `[working_directory]/mouse_cell_sizes/mouse_cell_sizes/case_set_1`)
* Go into the `human_cell_sizes.zip` and uzip each of the five `case_###.zip` in that folder.  Again, check the file structure (e.g., avoid `[working_directory]/human_cell_sizes/case_072/case_072/`)
* Delete the six zip files that you unzipped.  **Do not unzip or delete the `L2/3.zip` or `lines.zip` files in any of the case folders!**
* Run Code6 and Code7 from the working directory in order.  

All R code was run in Windows, but should produce similar results in UNIX.  See session info in above scripts for details.  


## Python scripts

***DESCRIPTION OF PYTHON SCRIPTS GOES HERE***
