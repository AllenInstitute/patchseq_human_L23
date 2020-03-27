# Overview
This folder contains the analysis scripts.  Prior to running these scripts, download all of them into your **main working directory** (*not* a scripts subdirector).  Then run the scripts in order and follow any additional directions indicated below or within the scripts.  

# List of scripts

This document is divided into two sections based on which programming language was used for the scripts.  The R scripts should be run prior to running the python scripts, as the python scripts required cluster calls for the human Patch-seq cells (which are geenrated using R).  If only the Python scripts are to be run, we provide cluster calls for Patch-seq cells in the `data` directory.

## R scripts

### (Placeholder for main R scripts, Rmd files).

`Code6_cortical_depth_code_human.r` and `Code7_cortical_depth_code_mouse.r` together perform quantifications of neuron density and soma size in mouse and human that are presented in Figure 1a-c.  **Before running these two scripts** data in the `data\cell_sizes` folder needs to be downloaded and unzipped in a specific way:
1. Download everything in `data\cell_sizes` and transfer to your working directory
2. Unzip `mouse_cell_sizes.zip` (make sure the directory is in the format `[working_directory]/mouse_cell_sizes/case_set_1` rather than `[working_directory]/mouse_cell_sizes/mouse_cell_sizes/case_set_1`)
3. Go into the `human_cell_sizes.zip` and uzip each of the five `case_###.zip` in that folder.  Again, check the file structure (e.g., avoid `[working_directory]/human_cell_sizes/case_072/case_072/`)
4. Delete the six zip files that you unzipped.  **Do not unzip or delete the `L2/3.zip` or `lines.zip` files in any of the case folders!**
5. Run Code6 and Code7 from the working directory in order.  

## Python scripts

***DESCRIPTION OF PYTHON SCRIPTS GOES HERE***
