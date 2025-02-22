'''READ ME'''
#README
#My folder is too large so I upload everything to my github repositories, please read my code and paper there.
# https://github.com/Garthzzz/Identifying-Regenerative-Cells-in-Xenopus-laevis-Using-Machine-Learning
'''READ ME'''

import os
import pandas as pd
import numpy as np

import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.metrics import silhouette_score, adjusted_rand_score, rand_score
from sklearn.model_selection import train_test_split
import xgboost as xgb
from sklearn.svm import LinearSVC
from sklearn.cluster import KMeans

# ---------------------------------------------------------------------
# User-defined constants
# ---------------------------------------------------------------------
TIMEPOINTS = ["0","1","2","3","all"]  
CLUSTER_METHODS = ["leiden","louvain","kmeans"]
MODELS = ["logreg","xgboost","svm"]  
TOP_N = 75  
EXP_THRESHOLD = 1.0  
PORTION_REQUIRED = 0.0   # changed from 0.6
RESULT_DIR = "results_new"  # main output directory

# Provide the absolute path to your Table3 Excel file
TABLE3_PATH = r"E:\5243\proj1\aav9996_tables3.xlsx"
TABLE3_SHEET = "ROC markers"

# ---------------------------------------------------------------------
# Function: get_subset_anndata
# ---------------------------------------------------------------------
def get_subset_anndata(adata, tp):
    """
    If tp == 'all', return the entire dataset (a copy).
    Otherwise, filter by adata.obs['DaysPostAmputation'] == tp.
    """
    if tp == "all":
        return adata.copy()
    else:
        return adata[adata.obs['DaysPostAmputation'].astype(str) == tp].copy()

# ---------------------------------------------------------------------
# Function: do_clustering_and_umap
# ---------------------------------------------------------------------
def do_clustering_and_umap(adata_sub, time_label):
    """
    Perform normalization, HVG, PCA on adata_sub.
    Then run Leiden, Louvain, and K-means.
    For each method, create individual UMAP plots and save them separately.
    """
    sc.pp.normalize_total(adata_sub, target_sum=1e4)
    sc.pp.log1p(adata_sub)
    sc.pp.highly_variable_genes(adata_sub, flavor='seurat', n_top_genes=2000)
    adata_sub = adata_sub[:, adata_sub.var['highly_variable']].copy()

    sc.tl.pca(adata_sub, svd_solver='arpack')
    sc.pp.neighbors(adata_sub, n_neighbors=10, n_pcs=20)
    sc.tl.umap(adata_sub)

    # 1) Leiden
    leiden_key = f"leiden_{time_label}"
    sc.tl.leiden(adata_sub, resolution=2.0, key_added=leiden_key)
    # Save UMAP figure
    outdir = os.path.join(RESULT_DIR, f"{time_label}")
    os.makedirs(outdir, exist_ok=True)
    sc.pl.umap(
        adata_sub,
        color=[leiden_key],
        title=[f"{time_label} Leiden"],
        show=False
    )
    plt.savefig(os.path.join(outdir, f"{time_label}_Leiden_umap.png"))
    plt.close()

    # 2) Louvain
    louvain_key = f"louvain_{time_label}"
    sc.tl.louvain(adata_sub, resolution=3.0, key_added=louvain_key)
    sc.pl.umap(
        adata_sub,
        color=[louvain_key],
        title=[f"{time_label} Louvain"],
        show=False
    )
    plt.savefig(os.path.join(outdir, f"{time_label}_Louvain_umap.png"))
    plt.close()

    # 3) K-means
    X_pca = adata_sub.obsm['X_pca']
    kmeans = KMeans(n_clusters=45, random_state=42)
    kmeans_labels = kmeans.fit_predict(X_pca)
    kmeans_key = f"kmeans_{time_label}"
    adata_sub.obs[kmeans_key] = kmeans_labels.astype(str)
    sc.pl.umap(
        adata_sub,
        color=[kmeans_key],
        title=[f"{time_label} K-means"],
        show=False
    )
    plt.savefig(os.path.join(outdir, f"{time_label}_Kmeans_umap.png"))
    plt.close()

    return adata_sub

# ---------------------------------------------------------------------
# Function: check_0_vs_123_expression
# ---------------------------------------------------------------------
def check_0_vs_123_expression(
    adata_0, adata_123, gene_list,
    threshold=EXP_THRESHOLD,
    portion=PORTION_REQUIRED
):
    """
    For each gene in gene_list, compute mean expression in adata_0 vs adata_123.
    If (mean123 - mean0) > threshold, we call it "abnormal".
    If at least portion (e.g. 0.2) of genes are abnormal, return True, else False.
    """
    genes_in_both = [g for g in gene_list if (g in adata_0.var_names) and (g in adata_123.var_names)]
    if len(genes_in_both) == 0:
        return False

    idx0 = [np.where(adata_0.var_names==g)[0][0] for g in genes_in_both]
    idx123 = [np.where(adata_123.var_names==g)[0][0] for g in genes_in_both]

    X0 = adata_0.X
    if not isinstance(X0, np.ndarray):
        X0 = X0.toarray()
    mean0 = X0[:, idx0].mean(axis=0)

    X123 = adata_123.X
    if not isinstance(X123, np.ndarray):
        X123 = X123.toarray()
    mean123 = X123[:, idx123].mean(axis=0)

    diff = mean123 - mean0
    count_abnormal = np.sum(diff > threshold)
    ratio = count_abnormal / len(genes_in_both)
    return (ratio >= portion)

# ---------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------
def main():
    # 1) Read the original data and do QC
    raw_path = r"E:\5243\proj1\cleaned_processed_frogtail.h5ad"
    adata_raw = sc.read_h5ad(raw_path)
    print(adata_raw)

    sc.pp.filter_cells(adata_raw, min_genes=200)
    sc.pp.filter_genes(adata_raw, min_cells=3)
    adata_raw.var["mt"] = adata_raw.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata_raw, qc_vars=['mt'], inplace=True)
    adata_raw = adata_raw[adata_raw.obs['pct_counts_mt']<20].copy()

    # Create result root
    os.makedirs(RESULT_DIR, exist_ok=True)

    # Prepare 0 vs 123 for expression check
    adata_0 = adata_raw[adata_raw.obs['DaysPostAmputation'].astype(str)=='0'].copy()
    adata_123 = adata_raw[adata_raw.obs['DaysPostAmputation'].astype(str).isin(['1','2','3'])].copy()

    # We'll store some cluster-level info in a list of dictionaries
    potential_roc_list = []

    # 2) For each of the five timepoints, do clustering, then for each method => markers
    for tp in TIMEPOINTS:
        adata_sub = get_subset_anndata(adata_raw, tp)
        if adata_sub.n_obs < 50:
            print(f"[{tp}] has only {adata_sub.n_obs} cells, skip.")
            continue

        # Clustering + separate UMAP plots
        adata_clustered = do_clustering_and_umap(adata_sub, tp)

        # Now do three cluster methods: we have columns like "leiden_{tp}", "louvain_{tp}", "kmeans_{tp}"
        for method in CLUSTER_METHODS:
            cluster_col = f"{method}_{tp}"
            if cluster_col not in adata_clustered.obs.columns:
                print(f"Warning: no col {cluster_col} in obs, skip.")
                continue

            # We'll store marker results in a file
            # We'll do 3 models: logreg (Scanpy), xgboost, svm
            out_dir = os.path.join(RESULT_DIR, f"{tp}_{method}")
            os.makedirs(out_dir, exist_ok=True)

            for model in MODELS:
                out_marker_file = os.path.join(out_dir, f"{tp}_{method}_{model}_marker.csv")

                if model == "logreg":
                    # use rank_genes_groups with method="logreg"
                    sc.tl.rank_genes_groups(adata_clustered, groupby=cluster_col,
                                            method="logreg", n_genes=500, max_iter=1000)
                    rgg = adata_clustered.uns['rank_genes_groups']
                    clusters = rgg['names'].dtype.names

                    recs = []
                    for clv in clusters:
                        # get top75
                        genes_all = rgg['names'][clv]
                        scores_all = rgg['scores'][clv]
                        n_take = min(TOP_N, len(genes_all))
                        genes_all = genes_all[:n_take]
                        scores_all = scores_all[:n_take]
                        for rank_i, (g, scv) in enumerate(zip(genes_all, scores_all)):
                            recs.append({
                                'time': tp,
                                'method': method,
                                'model': model,
                                'cluster': clv,
                                'gene': g,
                                'score': float(scv),
                                'rank_in_cluster': rank_i+1
                            })
                    df_marker = pd.DataFrame(recs)
                    df_marker.to_csv(out_marker_file, index=False)

                else:
                    # xgboost or svm => manual 1 vs rest
                    recs = []
                    X = adata_clustered.X
                    if not isinstance(X, np.ndarray):
                        X = X.toarray()
                    cluster_vals = adata_clustered.obs[cluster_col].unique()

                    for clv in cluster_vals:
                        labels = (adata_clustered.obs[cluster_col].values == clv).astype(int)
                        X_train, X_test, y_train, y_test = train_test_split(X, labels, test_size=0.2, random_state=42)
                        if model=="xgboost":
                            clf = xgb.XGBClassifier(n_estimators=100, max_depth=3,
                                                    random_state=42, eval_metric='logloss')
                            clf.fit(X_train, y_train)
                            importances = clf.feature_importances_
                            idx_sorted = np.argsort(importances)[::-1]
                            n_take = min(TOP_N, len(idx_sorted))
                            top_idx = idx_sorted[:n_take]
                            for rank_i, idxg in enumerate(top_idx):
                                recs.append({
                                    'time': tp,
                                    'method': method,
                                    'model': model,
                                    'cluster': clv,
                                    'gene': adata_clustered.var_names[idxg],
                                    'score': float(importances[idxg]),
                                    'rank_in_cluster': rank_i+1
                                })
                        else:
                            # SVM
                            svc = LinearSVC(C=1.0, max_iter=10000, random_state=42)
                            svc.fit(X_train, y_train)
                            coefs = svc.coef_.flatten()
                            idx_sorted = np.argsort(np.abs(coefs))[::-1]
                            n_take = min(TOP_N, len(idx_sorted))
                            top_idx = idx_sorted[:n_take]
                            for rank_i, idxg in enumerate(top_idx):
                                recs.append({
                                    'time': tp,
                                    'method': method,
                                    'model': model,
                                    'cluster': clv,
                                    'gene': adata_clustered.var_names[idxg],
                                    'score': float(coefs[idxg]),
                                    'rank_in_cluster': rank_i+1
                                })

                    df_marker = pd.DataFrame(recs)
                    df_marker.to_csv(out_marker_file, index=False)

                print(f"Saved marker => {out_marker_file}")

    # Next: do the "0 vs 123" expression filter
    print("\n=== Step: 0 vs 123 expression filter ===")
    adata_0 = adata_raw[adata_raw.obs['DaysPostAmputation'].astype(str)=='0'].copy()
    adata_123 = adata_raw[adata_raw.obs['DaysPostAmputation'].astype(str).isin(['1','2','3'])].copy()

    potential_roc_list = []

    for tp in TIMEPOINTS:
        for method in CLUSTER_METHODS:
            cluster_dir = os.path.join(RESULT_DIR, f"{tp}_{method}")
            if not os.path.isdir(cluster_dir):
                continue
            for model in MODELS:
                marker_file = os.path.join(cluster_dir, f"{tp}_{method}_{model}_marker.csv")
                if not os.path.isfile(marker_file):
                    continue
                df_marker = pd.read_csv(marker_file)
                # group by cluster
                for clv, subdf in df_marker.groupby('cluster'):
                    gene_list = subdf['gene'].unique().tolist()
                    is_potential = check_0_vs_123_expression(
                        adata_0, adata_123, gene_list,
                        threshold=EXP_THRESHOLD,
                        portion=PORTION_REQUIRED
                    )
                    if is_potential:
                        potential_roc_list.append({
                            'time': tp,
                            'method': method,
                            'model': model,
                            'cluster': clv,
                            'marker_count': len(gene_list)
                        })

    df_potential = pd.DataFrame(potential_roc_list)
    df_potential_path = os.path.join(RESULT_DIR, "potential_roc_clusters.csv")
    df_potential.to_csv(df_potential_path, index=False)
    print(f"[0 vs 123 filter] potential clusters => {df_potential_path}")

    # Next: read table3
    df_table3 = pd.read_excel(TABLE3_PATH, sheet_name=TABLE3_SHEET)
    roc_genes = set(df_table3['GeneSymbol'].unique())

    final_records = []

    for row in df_potential.itertuples():
        tp = row.time
        method = row.method
        model = row.model
        clv = row.cluster
        cluster_dir = os.path.join(RESULT_DIR, f"{tp}_{method}")
        marker_file = os.path.join(cluster_dir, f"{tp}_{method}_{model}_marker.csv")
        if not os.path.isfile(marker_file):
            continue
        df_marker = pd.read_csv(marker_file)
        subdf = df_marker[df_marker['cluster']==clv].copy()
        glist = subdf['gene'].unique().tolist()
        overlap = roc_genes.intersection(glist)
        overlap_count = len(overlap)
        ratio = overlap_count / len(glist) if len(glist)>0 else 0.0

        final_records.append({
            'time': tp,
            'method': method,
            'model': model,
            'cluster': clv,
            'marker_count': len(glist),
            'overlap_count': overlap_count,
            'overlap_ratio': ratio
        })

    df_final = pd.DataFrame(final_records)
    df_final_path = os.path.join(RESULT_DIR, "roc_overlap_summary.csv")
    df_final.to_csv(df_final_path, index=False)
    print(f"Overlap summary => {df_final_path}")

    # find best clusters in terms of overlap_count
    best_count_list = []
    for (tp, method, model), group in df_final.groupby(['time','method','model']):
        group2 = group.sort_values('overlap_count', ascending=False)
        if len(group2)>0:
            best_count_list.append(group2.iloc[0].to_dict())
    df_best_count = pd.DataFrame(best_count_list)
    df_best_count.sort_values('overlap_count', ascending=False, inplace=True)
    df_best_count_file = os.path.join(RESULT_DIR,"best_overlap_count.csv")
    df_best_count.to_csv(df_best_count_file, index=False)

    # find best clusters in terms of overlap_ratio
    best_ratio_list = []
    for (tp, method, model), group in df_final.groupby(['time','method','model']):
        group2 = group.sort_values('overlap_ratio', ascending=False)
        if len(group2)>0:
            best_ratio_list.append(group2.iloc[0].to_dict())
    df_best_ratio = pd.DataFrame(best_ratio_list)
    df_best_ratio.sort_values('overlap_ratio', ascending=False, inplace=True)
    df_best_ratio_file = os.path.join(RESULT_DIR,"best_overlap_ratio.csv")
    df_best_ratio.to_csv(df_best_ratio_file, index=False)

    # barplot for overlap_count
    plt.figure(figsize=(10,6))
    sns.barplot(data=df_best_count, x='time', y='overlap_count', hue='method')
    plt.title("Max Overlap Count across (time, method, model) - ignoring model dimension")
    plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    plt.savefig(os.path.join(RESULT_DIR,"barplot_overlap_count.png"), bbox_inches='tight')
    plt.close()

    # barplot for overlap_ratio
    plt.figure(figsize=(10,6))
    sns.barplot(data=df_best_ratio, x='time', y='overlap_ratio', hue='method')
    plt.title("Max Overlap Ratio across (time, method, model) - ignoring model dimension")
    plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
    plt.savefig(os.path.join(RESULT_DIR,"barplot_overlap_ratio.png"), bbox_inches='tight')
    plt.close()

    print("All done. Check your 'results_new' folder.")

if __name__=="__main__":
    main()
