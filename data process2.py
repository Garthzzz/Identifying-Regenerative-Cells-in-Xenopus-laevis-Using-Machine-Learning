import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

################################################################################
# 1. User Parameters
################################################################################
# Root directories that contain your four sets of outputs,
# e.g. "overlap filter ratio=0, resolution=1.0, kmeans n=10", etc.
# Suppose they are stored in a higher-level folder "E:\5243\proj1\Identification-of-Regenerative-Organizing-Cells..."
# and inside we have 4 subfolders. We'll list them here:
ROOT_DIRS = [
    r"E:\5243\proj1\Identification-of-Regenerative-Organizing-Cells-in-Xenopus-laevis-Using-Multi-Cluster-Analysis-and-Machine-Learning\overlap filter ratio=0, resolution=1.0, kmeans n=10",
    r"E:\5243\proj1\Identification-of-Regenerative-Organizing-Cells-in-Xenopus-laevis-Using-Multi-Cluster-Analysis-and-Machine-Learning\overlap filter ratio=0, resolution=2.0, kmeans n=45",
    r"E:\5243\proj1\Identification-of-Regenerative-Organizing-Cells-in-Xenopus-laevis-Using-Multi-Cluster-Analysis-and-Machine-Learning\overlap filter ratio=0.2, resolution=1.0, kmeans n=10",
    r"E:\5243\proj1\Identification-of-Regenerative-Organizing-Cells-in-Xenopus-laevis-Using-Multi-Cluster-Analysis-and-Machine-Learning\overlap filter ratio=0.2, resolution=2.0, kmeans n=45"
]

# The times, methods, models used. Usually 5x3 -> 15 combos, each combo has 3 models => total 45
TIMEPOINTS = ["0","1","2","3","all"]
METHODS = ["leiden","louvain","kmeans"]
MODELS = ["logreg","xgboost","svm"]

# Table3 path:
TABLE3_PATH = r"E:\5243\proj1\aav9996_tables3.xlsx"
TABLE3_SHEET = "ROC markers"

# Output directory to save final combined results:
FINAL_OUTPUT_DIR = r"E:\5243\proj1\FINAL_ANALYSIS"

################################################################################
# 2. Helper function: remove .L/.S suffix
################################################################################
def remove_LS_suffix(gene_name: str) -> str:
    """
    Remove trailing '.L' or '.S' if they exist.
    Example: 'fgf7.L' -> 'fgf7'
    """
    return re.sub(r'\.(L|S)$', '', gene_name)

################################################################################
# 3. Main analysis
################################################################################
def main():
    # Define final_records here to avoid NameError
    final_records = []
    os.makedirs(FINAL_OUTPUT_DIR, exist_ok=True)

    # 3.1 Read Table3 and create a set of ROC genes (all in lower-case or keep as is)
    df_table3 = pd.read_excel(TABLE3_PATH, sheet_name=TABLE3_SHEET)
    # assume the column is "GeneSymbol"
    # It's safer to unify case. But if not, at least we keep the original:
    roc_set = set(df_table3["GeneSymbol"].dropna().unique())

    # 3.2 We'll iterate over each of the 4 root directories:
    #     Each root dir has subfolders: "0_leiden", "0_louvain", "0_kmeans", "1_leiden", etc.
    # Our plan:
    #   For each (tp, method, model), read CSV => group by cluster => find overlap with table3 => track best cluster by count/ratio
    # Then do barplots. Then unify strong rocs genes. Then unify "very strong" across 15 combos, ...
    
    big_records = []  # store (root_name, time, method, model, cluster, overlap_count, overlap_ratio, cluster_markers_noSuffix) for all clusters
    # We'll also store best overlap by count & ratio in a separate structure:
    best_count_info = []
    best_ratio_info = []

    # we also want to keep track of "the cluster that has the max overlap" for each (tp, method, model).
    # Then we can unify them by "tp+method" for 3 models => find strong rocs gene. We'll do that after we gather everything.

    # We'll define an index to identify "which root dir" => e.g. " ratio=0, res=1.0, k=10" => short name
    # You can define a short_name or parse from path.
    def short_dir_name(path):
        return os.path.basename(path)  # or do something more advanced

    for root_dir in ROOT_DIRS:
        root_label = short_dir_name(root_dir)  # e.g. "overlap filter ratio=0, resolution=2.0, kmeans n=45"
        # We'll create a sub-output folder inside FINAL_OUTPUT_DIR to store partial results
        out_subdir = os.path.join(FINAL_OUTPUT_DIR, root_label)
        os.makedirs(out_subdir, exist_ok=True)

        # We'll gather a list for final summarization
        local_records = []

        # parse all timepoint/method combos => e.g. "0_leiden", "0_louvain", "0_kmeans", ...
        for tp in TIMEPOINTS:
            for meth in METHODS:
                # folder example: "{root_dir}/{tp}_{meth}"
                subfolder = os.path.join(root_dir, f"{tp}_{meth}")
                if not os.path.isdir(subfolder):
                    # skip
                    continue
                # Inside subfolder, we expect {tp}_{meth}_{model}_marker.csv
                for model in MODELS:
                    marker_file = os.path.join(subfolder, f"{tp}_{meth}_{model}_marker.csv")
                    if not os.path.isfile(marker_file):
                        continue

                    df_marker = pd.read_csv(marker_file)
                    # df_marker columns: time, method, model, cluster, gene, score, rank_in_cluster
                    # we group by cluster
                    for clust_id, subdf in df_marker.groupby('cluster'):
                        # remove .L/.S from all genes
                        raw_genes = subdf['gene'].unique().tolist()
                        # map remove_LS_suffix
                        processed_genes = [ remove_LS_suffix(g) for g in raw_genes ]
                        # compute overlap
                        overlap = roc_set.intersection(processed_genes)
                        overlap_count = len(overlap)
                        ratio = overlap_count / len(processed_genes) if len(processed_genes)>0 else 0.0

                        rec = {
                            'root_label': root_label,  # e.g. ratio=0, res=2.0, ...
                            'time': tp,
                            'method': meth,
                            'model': model,
                            'cluster': clust_id,
                            'marker_count': len(processed_genes),
                            'overlap_count': overlap_count,
                            'overlap_ratio': ratio,
                            # store processed gene list => we can save if this cluster is best
                            'processed_genes': processed_genes  
                        }
                        local_records.append(rec)

        # convert local_records => df
        df_local = pd.DataFrame(local_records)
        out_csv = os.path.join(out_subdir, "all_clusters_overlap.csv")
        df_local.to_csv(out_csv, index=False)

        # Now find best cluster by overlap_count & overlap_ratio for each (time, method, model)
        best_count_list = []
        best_ratio_list = []
        gcols = ['time','method','model']
        for gkey, gdf in df_local.groupby(gcols):
            gdf2_count = gdf.sort_values('overlap_count', ascending=False)
            if len(gdf2_count)>0:
                best_c = gdf2_count.iloc[0].to_dict()
                best_count_list.append(best_c)
            gdf2_ratio = gdf.sort_values('overlap_ratio', ascending=False)
            if len(gdf2_ratio)>0:
                best_r = gdf2_ratio.iloc[0].to_dict()
                best_ratio_list.append(best_r)

        df_best_count = pd.DataFrame(best_count_list)
        df_best_ratio = pd.DataFrame(best_ratio_list)

        df_best_count.to_csv(os.path.join(out_subdir, "best_overlap_count.csv"), index=False)
        df_best_ratio.to_csv(os.path.join(out_subdir, "best_overlap_ratio.csv"), index=False)

        # barplot for count
        plt.figure(figsize=(10,6))
        # use a stable sort for time => so we might do something like .sort_values
        df_best_count_sorted = df_best_count.sort_values('time')
        sns.barplot(data=df_best_count_sorted, x='time', y='overlap_count', hue='method')
        plt.title(f"Max Overlap Count in {root_label}")
        plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
        plt.savefig(os.path.join(out_subdir,"barplot_best_overlap_count.png"), bbox_inches='tight')
        plt.close()

        # barplot for ratio
        plt.figure(figsize=(10,6))
        df_best_ratio_sorted = df_best_ratio.sort_values('time')
        sns.barplot(data=df_best_ratio_sorted, x='time', y='overlap_ratio', hue='method')
        plt.title(f"Max Overlap Ratio in {root_label}")
        plt.legend(bbox_to_anchor=(1.05,1), loc='upper left')
        plt.savefig(os.path.join(out_subdir,"barplot_best_overlap_ratio.png"), bbox_inches='tight')
        plt.close()

        # We also want to save the "best cluster" genes themselves into a folder
        best_genes_dir = os.path.join(out_subdir, "best_clusters_genes")
        os.makedirs(best_genes_dir, exist_ok=True)

        for idx, row in df_best_count.iterrows():
            # store the cluster that has the max overlap_count
            # we'll keep them in a text/csv file
            t = row['time']
            mth = row['method']
            mod = row['model']
            cl = row['cluster']
            genes = row['processed_genes'] if isinstance(row['processed_genes'], list) else []
            # write them
            fname = f"{t}_{mth}_{mod}_{cl}_bestCount_genes.txt"
            outpath = os.path.join(best_genes_dir, fname)
            with open(outpath,'w',encoding='utf-8') as fw:
                for g in genes:
                    fw.write(g+"\n")

        for idx, row in df_best_ratio.iterrows():
            # store the cluster that has the max overlap_ratio
            t = row['time']
            mth = row['method']
            mod = row['model']
            cl = row['cluster']
            genes = row['processed_genes'] if isinstance(row['processed_genes'], list) else []
            fname = f"{t}_{mth}_{mod}_{cl}_bestRatio_genes.txt"
            outpath = os.path.join(best_genes_dir, fname)
            with open(outpath,'w',encoding='utf-8') as fw:
                for g in genes:
                    fw.write(g+"\n")

        # Now we do the "strong rocs gene" logic for each (time+method).
        # i.e. for each (time, method), we have 3 models => best cluster by overlap_count => gather their genes
        # if a gene appears >= 2 times among the 3 => strong rocs gene
        # We'll call them "best_count" for now, but you can similarly do "best_ratio" if you want.

        # group df_best_count by (time, method)
        strong_records = []
        for (tp, meth), subg in df_best_count.groupby(['time','method']):
            # subg has up to 3 rows (for the 3 models)
            # let's gather gene sets
            all_genes_list = []
            for i2, r2 in subg.iterrows():
                gset = set(r2['processed_genes']) if isinstance(r2['processed_genes'], list) else set()
                all_genes_list.append(gset)
            # We want to see how many times each gene appears across the 3 sets
            from collections import Counter
            c_all = Counter()
            for gset in all_genes_list:
                for gx in gset:
                    c_all[gx]+=1
            # strong rocs gene => appear >=2 times
            strong_genes = [k for k,v in c_all.items() if v>=2]
            # record
            rec_st = {
                'root_label': root_label,
                'time': tp,
                'method': meth,
                'strong_genes': strong_genes,
                'strong_count': len(strong_genes)
            }
            strong_records.append(rec_st)

        df_strong = pd.DataFrame(strong_records)
        df_strong_out = os.path.join(out_subdir, "strong_rocs_genes_per_timeMethod.csv")
        df_strong.to_csv(df_strong_out, index=False)

        # now we store each group's strong genes in a text
        strong_dir = os.path.join(out_subdir, "strong_genes_text")
        os.makedirs(strong_dir, exist_ok=True)
        for idx, row2 in df_strong.iterrows():
            tp = row2['time']
            meth = row2['method']
            sgenes = row2['strong_genes']
            fname = f"{tp}_{meth}_strong_genes.txt"
            pathx = os.path.join(strong_dir, fname)
            with open(pathx,'w',encoding='utf-8') as fw:
                for g in sorted(sgenes):
                    fw.write(g+"\n")

        # we now keep these results for global union
        # so that we can do "very strong rocs gene" across 15 groups
        # each (tp, method) is 1 group, => total 5x3=15
        # but we also want to track root_label => so let's define a unique group ID: (root_label, time, method)
        # each group => we have strong_genes
        for idx, row2 in df_strong.iterrows():
            strong_glist = row2['strong_genes']
            final_records.append({
                'root_label': root_label,
                'time': row2['time'],
                'method': row2['method'],
                'strong_genes': strong_glist
            })

    # 4) Now we gather all (root_label, time, method) => up to 4 (root) x 5 (time) x 3 (method) = 60 groups
    # or maybe you consider each "time+method" ignoring root_label => user wants to treat same (time,method) across different root settings as distinct or the same?
    # The user says "any gene appear in >= 45 groups => very strong rocs gene."
    # but we have possibly 4 * 15 = 60 groups total, or if you treat "time+method" as 15 groups. The user specifically said "any gene appear in 9 or more => very strong." So presumably out of 15 or 60?


    all_strong_map = {}  # key: groupkey -> set of strong_genes
    for rec in final_records:
        groupkey = (rec['root_label'], rec['time'], rec['method'])
        sg = set(rec['strong_genes'])
        all_strong_map[groupkey] = sg

    # Now we want to count how many groups each gene appears in
    from collections import Counter
    global_counter = Counter()
    for groupkey, sgeneset in all_strong_map.items():
        for g in sgeneset:
            global_counter[g]+=1

    # We define "very strong" if a gene appears in >=30 groups
    very_strong = [k for k,v in global_counter.items() if v>=30]

    # Save them
    vs_path = os.path.join(FINAL_OUTPUT_DIR,"very_strong_rocs_genes.txt")
    with open(vs_path,'w',encoding='utf-8') as fw:
        for g in sorted(very_strong):
            fw.write(g+"\n")

    print(f"Very strong rocs genes (>=30 groups) saved => {vs_path}")


if __name__=="__main__":
    main()
