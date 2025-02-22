# Identifying Regenerative Organizing Cells in Xenopus laevis Using Multi-Cluster Analysis and Machine Learning

This repository contains all materials for our project on discovering and validating Regenerative Organizing Cells (ROCs) in *Xenopus laevis* tadpole tails, following up on the dataset and findings from [Aztekin *et al.*, *Science* **364**, 653–658 (2019)](https://doi.org/10.1126/science.aav9996).

## Repository Structure

1. **Identifying Regenerative Organizing Cells in Xenopus laevis Using Multi-Cluster Analysis and Machine Learning.pdf**  
   - Our primary paper (in PDF format). It describes the motivation, data processing pipeline, clustering and marker analyses, results, and discussion. This document provides a comprehensive overview of the entire project.

2. **aav9996_tables3.xlsx**  
   - The supplementary ROC gene list (Table 3) from the original publication, used for overlap comparisons with the cluster-derived marker genes.

3. **ratio=0, res=1.0, km n=10 / ratio=0, res=2.0, km n=45 / ratio=0.2, res=1.0, km n=10 / ratio=0.2, res=2.0, km n=45**  
   - Each folder represents one configuration of filtering ratio (0 or 0.2), Leiden/Louvain resolution (1.0 vs. 2.0/3.0), and K-means cluster number (10 vs. 45).  
   - Inside each, you will find subfolders named `0_leiden`, `0_louvain`, `0_kmeans`, etc., for the five time points (`0`, `1`, `2`, `3`, `all`) combined with each clustering method.  
   - Each sub-subfolder contains top marker results for **logreg**, **xgboost**, and **svm** (e.g., `0_leiden_logreg_marker.csv`), plus UMAP images.

4. **FINAL_ANALYSIS**  
   - A folder containing aggregated results after merging all configurations. This includes `all_clusters_overlap.csv`, best overlap count/ratio CSVs, bar plots of overlaps, “strong” ROC genes per `(time, method)`, and final “very strong” or “super strong” ROIs.  

5. **data analysis.py**  
   - A Python script that **re-processes** or **analyzes** the marker files generated in the above subfolders.  
   - Iterates through the four root folders (each config), collects cluster–gene overlaps with the known ROC gene set, and computes summary statistics (best overlap count, best overlap ratio, bar plots by time/method/model).  
   - Also merges “best clusters” from each configuration to identify “strong” ROC genes and saves them in subdirectories.

6. **data analysis2.py**  
   - A second Python script focusing on a variant or extended analysis flow. It again scans the same root directories and produces additional consolidated outputs in a final folder. This includes “very strong” genes that appear in many groups across multiple cluster methods and filtering conditions.  

7. **data process2.py**  
   - A script demonstrating the main single-cell pipeline for each time point, performing QC, clustering (Leiden/Louvain/K-means), and marker selection (logreg, XGBoost, SVM). It also includes the optional 0 vs. 1,2,3 expression filter for identifying potential ROC clusters prior to overlap with Table 3.


## Reference

- **Aztekin, C., et al.** *Identification of a regeneration-organizing cell in the Xenopus tail.* *Science* 364, 653–658 (2019).  
  [DOI: 10.1126/science.aav9996](https://doi.org/10.1126/science.aav9996)

## License

This project is for academic purposes. Please cite appropriately if you adapt or reuse any portion of this code or methodology.

## Contact

For questions or collaborations, please contact [Zhengze Zhang](mailto:zhangzhengze2018@gmail.com) or open an issue in this repository. 
