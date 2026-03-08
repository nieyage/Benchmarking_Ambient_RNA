import scanpy as sc  # scanpy 1.9.3; anndata 0.12.6 
import pandas as pd
# pip install openpyxl
import numpy as np
import scipy.sparse as sp
from cellbender.remove_background.downstream import anndata_from_h5
import anndata as ad
import pickle
from sklearn.metrics import r2_score, auc
from scipy.stats import kendalltau, wilcoxon
import glob, os

import plotly.express as px

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
mpl.use('Agg') 

import seaborn as sns
sns.set_theme(style="white")


def load_cellnames(contam_ratio_path):
    df = pd.read_csv(contam_ratio_path, index_col=0)
    non_osn_cells = df.index[df['celltype'] == 'nonOSN'].tolist()
    osn_cells = df.index[df['celltype'] == 'OSN'].tolist()
    return non_osn_cells, osn_cells


def calc_or_fraction_per_non_osn(adata, genes_of_interest, umi_threshold=0):
    """
    This function is used for calculating cell fraction and umi fraction of non-osns. In non-osns, cells are considered to not expressing or-specific genes and regarded as contamination.
    for non-osn cells: 
    1. cell fraction: proportion of cells expressing Omp (and other OR genes)
    2. UMI fraction: (UMI of all Omp and other OR genes expression) / (UMI of all genes)

    adata: denoised count data, raw represents barcodes, column represents genes.

    genes_of_interest: genes remained to calculate cell and UMI fraction.

    umi_threshold: threshold of cell UMI. Cell with a UMI higher than threshold is considered to positively expressing and counted.
    """
    # make all genes unique
    adata.var_names_make_unique()
    valid_genes = [g for g in adata.var_names if g.lower() in [x.lower() for x in genes_of_interest]] # filter genes we are interested 

    # extract count data
    expr = adata[:, valid_genes].X
    if sp.issparse(expr):
        expr = expr.toarray()

    # task 1: cell fraction exprssing Omp(and other OR genes)
    if expr.ndim == 1 or expr.shape[1] == 1:  # single genes
        positive_cells = expr.flatten() > umi_threshold
    else:  # multiple genes
        positive_cells = (expr > umi_threshold).any(axis=1)
    non_osn_cell_fraction = np.sum(positive_cells) / len(positive_cells) * 100

    # task 2: umi fraction exprssing Omp (and other OR genes)
    # raw count data
    total_expr = adata.X
    if sp.issparse(total_expr):
        total_expr = total_expr.toarray()
    non_osn_umi_fraction = np.sum(expr, axis=1) / np.sum(total_expr, axis=1) * 100

    return non_osn_cell_fraction, non_osn_umi_fraction


def calc_or_fraction_per_osn(adata, genes_of_interest, umi_threshold=2):
    """
    Vectorized version for OSN cells:
    1. cell fraction of expressing other non-specific OR genes
    2. UMI fraction: (UMI of other OR) / (UMI of all OR genes)
    """

    print(len(adata.obs_names))

    adata.var_names_make_unique()
    valid_genes = [g for g in genes_of_interest if g in adata.var_names]

    expr = adata[:, valid_genes].X  
    if sp.issparse(expr):
        expr = expr.toarray()  # to dense format

    # cell-specific OR
    specific_or_idx = np.argmax(expr, axis=1) # column index of maximum expression value in each cell
    max_vals = expr[np.arange(expr.shape[0]), specific_or_idx] # the maximum expression value of each cell

    # OR genes expressing in osn more than or equal to umi_threshold=2 are considered as cell-specific OR genes
    valid_cells = max_vals >= umi_threshold   # A boolean array of length equal to the number of cells. [True, False, True, True, False, ...] 

    # print(len(valid_cells))
    print(valid_cells.sum())
    # invalid_barcodes = adata.obs_names[~valid_cells]
    # invalid_max_vals = max_vals[~valid_cells]
    # for bc, mv in zip(invalid_barcodes, invalid_max_vals):
    #     print(bc, mv)

    if valid_cells.sum() == 0:
        return np.nan, np.array([])

    expr_valid = expr[valid_cells, :]
    specific_idx_valid = specific_or_idx[valid_cells]

    other_umi = expr_valid.sum(axis=1) - expr_valid[np.arange(expr_valid.shape[0]), specific_idx_valid]
    total_umi = expr_valid.sum(axis=1) # overall umi of all OR genes
    osn_umi_fraction = other_umi.astype(float) / total_umi * 100

    # print(len(osn_umi_fraction))

    osn_cell_fraction = (osn_umi_fraction > 0).mean() * 100

    return osn_cell_fraction, osn_umi_fraction


def calc_fraction(path, non_osn_cells, osn_cells, method, non_osn_genes_of_interest, osn_genes_of_interest):
    
    if method == "CellClear":
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
    else:
        adata = sc.read_h5ad(path)

    non_osn_cell_fraction, non_osn_umi_fraction = np.nan, np.array([])
    valid_non_osn = list(set(non_osn_cells) & set(adata.obs_names))   # intersect with non-osn of corrected matrix

    if len(valid_non_osn) > 0:
        non_osn_adata = adata[valid_non_osn, :].copy()
        non_osn_cell_fraction, non_osn_umi_fraction = calc_or_fraction_per_non_osn(non_osn_adata, non_osn_genes_of_interest)

    osn_cell_fraction, osn_umi_fraction = np.nan, np.array([])
    valid_osn = list(set(osn_cells) & set(adata.obs_names))    # intersect with osn of corrected matrix

    # print(len(valid_osn))

    if len(valid_osn) > 0:
        osn_adata = adata[valid_osn, :].copy()
        osn_cell_fraction, osn_umi_fraction = calc_or_fraction_per_osn(osn_adata, osn_genes_of_interest)

    return {
        "non_osn_cell_fraction": non_osn_cell_fraction,
        "non_osn_umi_fraction": non_osn_umi_fraction,
        "osn_cell_fraction": osn_cell_fraction,
        "osn_umi_fraction": osn_umi_fraction,
    }


def calc_fractions_multiple_datasets(dataset_configs, methods, non_osn_genes_of_interest, osn_genes_of_interest):
    """"
    methods: path config method
    """
    non_osn_cell_fraction_results = {}
    non_osn_umi_fraction_results = {}
    osn_cell_fraction_results = {}
    osn_umi_fraction_results = {}

    for dataset_name, config in dataset_configs.items():
        print(f"Processing dataset: {dataset_name}")

        # load non-osn and osn according to contamination ratio csv files(row)
        non_osn_cells, osn_cells = load_cellnames(config["contam_ratio"])
        # print(len(osn_cells))

        non_osn_cell_fraction_results[dataset_name] = {}
        non_osn_umi_fraction_results[dataset_name] = {}
        osn_cell_fraction_results[dataset_name] = {}
        osn_umi_fraction_results[dataset_name] = {}

        for method in methods:
            print(f"Processing method: {method}")
            path = config["methods"][method]

            results = calc_fraction(path, non_osn_cells, osn_cells, method, non_osn_genes_of_interest, osn_genes_of_interest)

            non_osn_cell_fraction_results[dataset_name][method] = results["non_osn_cell_fraction"]
            non_osn_umi_fraction_results[dataset_name][method] = results["non_osn_umi_fraction"]
            osn_cell_fraction_results[dataset_name][method] = results["osn_cell_fraction"]
            osn_umi_fraction_results[dataset_name][method] = results["osn_umi_fraction"]

    return (
        non_osn_cell_fraction_results,
        non_osn_umi_fraction_results,
        osn_cell_fraction_results,
        osn_umi_fraction_results,
    )


# matplotlib type 
def plot_cell_fraction_bar(cell_fraction_results, datasets_to_plot, output, ylabel_title, plot_title, figsize=(16,8)):
    """
    Plot bar chart of cell fraction results for selected datasets.
    
    Parameters:
    -----------
    cell_fraction_results : dict
        Nested dict of the form {dataset: {method: fraction, ...}, ...}
    datasets_to_plot : list, optional
        List of datasets to plot. If None, plot all.
    output : str
        Path to save the figure
    ylabel_title : str
        Y-axis label
    plot_title : str
        Figure title
    figsize : tuple
        Figure size
    """
    df = pd.DataFrame(
        [(dataset, method, fraction) 
         for dataset, methods in cell_fraction_results.items() 
         for method, fraction in methods.items()],
        columns=["Dataset", "Method", "Fraction"]
    )
    # Filter datasets if specified
    if datasets_to_plot is not None:
        df = df[df["Dataset"].isin(datasets_to_plot)]

    palette = {
        "rawdata": "#BFBFBF",
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    # Get the list of unique datasets and methods
    datasets = list(df['Dataset'].unique())
    methods = list(df['Method'].unique())
    n_methods = len(methods)

    fig, ax = plt.subplots(figsize=figsize)
    bar_width = 0.1  # width of each bar
    gap = 0.2 # horizontal gap between different datasets
    offsets = [] # store x-axis offsets for all bars
    labels_pos = [] # store label positions for dataset names

    for i, dataset in enumerate(datasets):
        # Select data for this dataset
        sub = df[df['Dataset']==dataset].sort_values('Fraction', ascending=False)
        # base_x = i * (n_methods * (bar_width+0.05) + gap)  # starting position of x for this dataset group (for libraries)
        base_x = i * (n_methods * bar_width + gap)  # starting position of x for this dataset group

        labels_pos.append(base_x + (n_methods - 1) * bar_width / 2)

        for j, (idx, row) in enumerate(sub.iterrows()):
            # xi = base_x + j * (bar_width+0.05) # x-coordinate for this bar in current group
            xi = base_x + j * bar_width
            offsets.append(xi) # store position
            bar = ax.bar(xi, row['Fraction'], color=palette[row['Method']], width=bar_width, label=row['Method'])
            # ax.bar_label(bar, fmt="%.1f", fontsize=6, padding=2, rotation=90)
    
    # Set dataset names as x-axis tick labels
    ax.set_xticks(labels_pos)
    ax.set_xticklabels(datasets, fontsize=8)
    ax.set_ylabel(ylabel_title, fontsize=8)

    legend_order = ["rawdata", "DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]     # Define fixed display order of methods in legend
    handles, labels = ax.get_legend_handles_labels()    # Extract handles and labels from plotted bars
    ordered_handles = [handles[labels.index(l)] for l in legend_order if l in labels] # Reorder handles to match desired order
    ax.legend(ordered_handles, legend_order, title="Method", bbox_to_anchor=(1.05, 1), loc="upper left", prop={'size':8})     # Draw legend on the right side of the figure
    ax.get_legend().get_title().set_fontsize(8)

    ax.set_title(plot_title, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()    # Adjust layout to fit all elements neatly
    plt.savefig(output, dpi=300)
    plt.close()



def plot_decont_efficiency_and_diff(cell_fraction_results,
                                    datasets_to_plot,data_palette,
                                    method_order,
                                    outdir_df_eff,outdir_df_eff_var,
                                    outdir_eff_stripplot,outdir_eff_barplot,outdir_eff_lineplot, outdir_eff_radar,
                                    outdir_diff,outdir_var,
                                    outdir_eff_meanplot
                                    ):
    """
    cell_fraction_results: 
    datasets_to_plot: which datasets need to plot
    data_palette: dataset color
    method_order: the order of methods to plot
    outdir_df_eff: the decontamination efficacy csv (.csv)
    outdir_df_eff_var: the decontamination efficacy variance (.csv)

    outdir_eff_stripplot: output path of decontamination efficacy stripplot. The plot: x-axis is dataset, y-axis is decontamination efficacy, each dataset includes different methods.
    outdir_eff_barplot: output path of decontamination efficacy barplot. The plot: x-axis is method, y-axis is decontamination efficacy, each method includes different datasets.
    outdir_eff_lineplot: output path of decontamination efficacy lineplot. The plot: x-axis is dataset, y-axis is decontamination efficacy, each dataset includes different methods, lines within methods.
    outdir_eff_radar: output path of decontamination efficacy radar plot. 
    outdir_diff: output path of decontamination efficacy difference barplot(only in two datasets). x-axis is method, y-axis is difference of decontamination efficacy between two datasets.
    outdir_var: output path of decontamination efficacy variance barplot. x-axis is method, y-axis is variance of decontamination efficacy across different datasets in ascending order.
    outdir_eff_meanplot: path of mean decontamination efficacy barplot. x-axis is method, y-axis is mean decontamination efficacy across different datasets in decending order. Each black point means dataset.
    """
    method_palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    df = pd.DataFrame(
        [(dataset, method, fraction) 
         for dataset, methods in cell_fraction_results.items() 
         for method, fraction in methods.items()],
        columns=["Dataset", "Method", "Fraction"]
    )
    # print(cell_fraction_results)

    # Compute Decontamination efficacy 
    efficacy_data = []
    for ds in datasets_to_plot:
        sub = df[df["Dataset"] == ds]
        raw_val = sub[sub["Method"] == "rawdata"]["Fraction"].values[0]
        for _, row in sub.iterrows():
            if row["Method"] == "rawdata":
                eff = 0  # baseline
            else:
                eff = (raw_val - row["Fraction"]) / raw_val * 100
            efficacy_data.append({
                "Dataset": ds,
                "Method": row["Method"],
                "Efficacy": eff
            })
    df_eff = pd.DataFrame(efficacy_data) 
    df_eff = df_eff[df_eff["Method"] != "rawdata"] # exclude rawdata
    
    df_eff_raw = df_eff.copy()  # data names order
    df_eff = df_eff.sort_values("Efficacy", ascending=False) #sorted by efficacy

    # print(df_eff)
    df_eff.to_csv(outdir_df_eff)

    # ------------------------------
    # Plot efficacy stripplot 
    def efficacy_stripplot():
        # x is dataset

        plt.figure(figsize=(8,8))
        ax = sns.stripplot(
            data=df_eff,
            x="Dataset",
            y="Efficacy",
            hue="Method",
            palette=method_palette,
            dodge=False,
            jitter=0.2,        
            size=9,            
            order=datasets_to_plot,
            hue_order=method_order,
            linewidth=0.7
        )

        df_mean = df_eff.groupby("Dataset", as_index=False)["Efficacy"].mean()
        for i, dataset in enumerate(datasets_to_plot):
            mean_val = df_mean.loc[df_mean["Dataset"] == dataset, "Efficacy"].values[0]
            ax.scatter(i, mean_val, color="black", s=120, zorder=5, label="_nolegend_")

        # Wilcoxon test
        if len(datasets_to_plot) == 2:
                ds1, ds2 = datasets_to_plot
                merged = (
                    df_eff.pivot(index="Method", columns="Dataset", values="Efficacy")
                    .dropna(subset=[ds1, ds2])
                )
                stat, pval = wilcoxon(merged[ds1], merged[ds2])

                p_str = f"{pval:.2e}"
                y_max = df_eff["Efficacy"].max()
                y_sig = y_max * 1.1
                x1, x2 = 0, 1
                ax.plot([x1, x1, x2, x2], [y_sig-1, y_sig, y_sig, y_sig-1], lw=1.2, c='k')
                ax.text((x1 + x2) / 2, y_sig + 1, f"p = {p_str}", ha='center', va='bottom', fontsize=12, color='black')
                ax.set_ylim(top=y_max * 1.25)

        ax.get_legend().remove()
        plt.xlim(-0.5, len(datasets_to_plot) - 0.5)
        plt.ylabel("Percentage (%)")
        plt.title("Decontamination Efficacy")
        plt.tight_layout()

        if outdir_eff_stripplot is not None:
            plt.savefig(outdir_eff_stripplot)
        plt.close()
    # ------------------------------
    
    # ------------------------------
    # Plot efficacy barplot
    def efficacy_barplot():
        # x is method   
        
        plt.figure(figsize=(8, 8))
        # ax = sns.barplot(
        #     data=df_eff, 
        #     x="Method", 
        #     y="Efficacy", 
        #     hue="Dataset", 
        #     palette=data_palette,
        #     order=method_order,
        #     hue_order=datasets_to_plot
        # )

        # mean values of each method across datasets in descend order
        method_means = (
            df_eff.groupby('Method')['Efficacy']
            .mean()
            .sort_values(ascending=False)
        )   
        method_order_new = method_means.index.tolist()

        x = np.arange(len(method_order_new))
        total_width = 0.8  

        for i, m in enumerate(method_order_new):
            sub = df_eff[df_eff['Method'] == m]
            sub = sub.set_index("Dataset").reindex(datasets_to_plot).dropna(subset=["Efficacy"])
        
            n = len(sub)
            bar_width = total_width / n 
            offsets = np.linspace(-total_width/2 + bar_width/2, total_width/2 - bar_width/2, n)
            for j, (ds, val) in enumerate(zip(sub.index, sub["Efficacy"])):
                plt.bar(
                    x[i] + offsets[j],
                    val,
                    width=bar_width,
                    color=data_palette.get(ds, 'gray'),
                    label=ds if i == 0 else ""
                )
        plt.xticks(x, method_order_new)
        plt.ylabel("Percentage (%)")
        plt.title("Decontamination Efficacy")
        plt.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()

        if outdir_eff_barplot is not None:
            plt.savefig(outdir_eff_barplot)
        plt.close()

    # ------------------------------
    # Plot efficacy lineplot
    def efficacy_lineplot():
        plt.figure(figsize=(8,8))
        ax = sns.lineplot(
            data=df_eff,
            x="Dataset",
            y="Efficacy",
            hue="Method",
            palette=method_palette,
            marker="o",
            linewidth=1.5,
            markersize=8
        )
        ax.set_title("Decontamination efficacy ", fontsize=8)
        ax.set_ylabel("Percentage (%)", fontsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        if outdir_eff_lineplot is not None:
            plt.savefig(outdir_eff_lineplot)
        plt.close()

    # ------------------------------
    # Plot efficacy radar plot

    # print(df_eff_raw)
    def efficacy_radar():
        # it requires that 'Dataset' column of 'df_plot' are sorted by dataset name.
        fig = px.line_polar(
            df_eff_raw,
            r="Efficacy",
            theta="Dataset",
            color="Method",
            line_close=True,
            markers=True,
            color_discrete_map=method_palette
        )

        # remove background color
        fig.update_layout(
            title=f"Radar plot of decontamination efficacy",  # radar plot title
            paper_bgcolor="white",   # whole figure background
            plot_bgcolor="white",    # plot area background
            font=dict(color="black"),  # set the color of axis label, legend and so on as black
            polar=dict(
                bgcolor="white",              # radar plot background
                radialaxis=dict(gridcolor="lightgrey"),  # radial grid circles
                angularaxis=dict(gridcolor="lightgrey")  # angular grid lines
            )
        )
        if outdir_eff_radar is not None:
            fig.write_image(outdir_eff_radar)

    # ------------------------------
    # Plot efficacy difference barplot
    def efficacy_diff_barplot():

        # Compute efficacy difference between two datasets for each method
            
        if len(datasets_to_plot) == 2:
            pivot = df_eff.pivot(index="Method", columns="Dataset", values="Efficacy")
            pivot_diff = (pivot[datasets_to_plot[0]] - pivot[datasets_to_plot[1]]).abs().sort_values(ascending=False)
            
            plt.figure(figsize=(8,8))
            ax = sns.barplot(
                x=pivot_diff.index,
                y=pivot_diff.values,
                palette=[method_palette.get(m, "gray") for m in pivot_diff.index]
            )

            plt.ylabel("Percentage (%)")
            plt.title("The difference of decontamination efficacy across libraries")
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()

            if outdir_diff is not None:
                plt.savefig(outdir_diff)
            plt.close()

    # ------------------------------
    # Compute decontamination efficacy variance arcoss different datasets for each method
    def efficacy_var():
        var_df = (
                df_eff.groupby("Method", as_index=False)["Efficacy"]
                .var()
                .rename(columns={"Efficacy": "Efficacy_var"})
        )
        var_df = var_df.sort_values("Efficacy_var", ascending=True)
        colors = [method_palette[m] for m in var_df["Method"].tolist()]

        var_df.to_csv(outdir_df_eff_var)

        plt.figure(figsize=(8, 8))
        ax = sns.barplot(
            data=var_df,
            x="Method",
            y="Efficacy_var",
            order=var_df["Method"].tolist(),
            palette=colors
        )

        ax.set_title("Variance of decontamination efficacy", fontsize=8)
        ax.set_ylabel("Variance", fontsize=8)

        ax.tick_params(axis='x', rotation=45, labelsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        if outdir_var is not None:
            plt.savefig(outdir_var, dpi=300, bbox_inches="tight")
        plt.close()

    # ------------------------------
    # Compute mean decontamination efficacy across different datasets for each method
    def mean_efficacy_barplot():
        # mean and std
        df_stats = (
            df_eff.groupby("Method", as_index=False)
            .agg({'Efficacy': ['mean', 'std']})
        )
        df_stats.columns = ['Method', 'Mean_Efficacy', 'Std_Efficacy']

        # decrease by mean
        df_stats = df_stats.sort_values("Mean_Efficacy", ascending=False)
        # print(df_stats)

        method_order = df_stats["Method"].tolist()
        df_eff["Method"] = pd.Categorical(df_eff["Method"], categories=method_order, ordered=True)

        plt.figure(figsize=(8, 8))
        x = np.arange(len(method_order))
        bar_width = 0.6

        plt.bar(
            x,
            df_stats["Mean_Efficacy"],
            yerr=df_stats["Std_Efficacy"],
            capsize=5,
            color=[method_palette.get(m, "gray") for m in method_order],
            edgecolor="black",
            width=bar_width
        )

        for i, method in enumerate(method_order):
            sub = df_eff[df_eff["Method"] == method]
            jitter = np.linspace(-0.08, 0.08, len(sub))
            for j, eff in enumerate(sub["Efficacy"]):
                plt.scatter(
                    x[i] + jitter[j],
                    eff,
                    s=80,
                    color="black",
                    edgecolor="white",
                    zorder=10
                )

        plt.xticks(x, method_order, rotation=30, ha="right")
        plt.ylabel("Decontamination Efficacy (%)")
        plt.title("Mean Decontamination Efficacy per Method (Descending Order)")
        plt.tight_layout()

        if outdir_eff_meanplot is not None:
            plt.savefig(outdir_eff_meanplot)
        plt.close()

    efficacy_stripplot()
    efficacy_barplot()
    efficacy_lineplot()
    efficacy_radar()
    efficacy_diff_barplot()
    efficacy_var()
    mean_efficacy_barplot()



def plot_umi_fraction_boxplots(umi_fractions_results, datasets_to_plot, output, ylabel_title, plot_title, figsize=(16,8)):
    palette = {
        "rawdata": "#BFBFBF",
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    if datasets_to_plot is not None:
        datasets = [ds for ds in umi_fractions_results.keys() if ds in datasets_to_plot]
    else:
        datasets = list(umi_fractions_results.keys())

    methods = list(next(iter(umi_fractions_results.values())).keys())
    n_methods = len(methods)
    n_datasets = len(datasets)

    fig, ax = plt.subplots(figsize=figsize)
    box_width = 0.1
    inner_gap = 0.05   # spacing between methods within dataset
    outer_gap = 0.4    # spacing between datasets
    labels_pos = []

    for i, dataset in enumerate(datasets):
        # sort methods by median for that dataset
        method_medians = {m: np.median(umi_fractions_results[dataset][m]) for m in methods}
        methods_sorted = sorted(methods, key=lambda m: method_medians[m], reverse=True)

        group_width = n_methods * (box_width + inner_gap)
        base_x = i * (group_width + outer_gap)
        labels_pos.append(base_x + group_width / 2 - inner_gap / 2)

        # draw boxplots
        for j, method in enumerate(methods_sorted):
            vals = umi_fractions_results[dataset][method]
            pos = base_x + j * (box_width + inner_gap)
            bplot = ax.boxplot(
                vals,
                positions=[pos],
                widths=box_width,
                patch_artist=True,
                showfliers=False
            )
            for patch in bplot['boxes']:
                patch.set_facecolor(palette[method])
                patch.set_linewidth(0.6)
            for median in bplot['medians']:
                median.set(color='k', linewidth=0.8)


    ax.set_xticks(labels_pos)
    ax.set_xticklabels(datasets, fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel_title, fontsize=8)
    ax.set_ylim(bottom=0)
    
    x_min = min(labels_pos) - (n_methods * box_width)
    x_max = max(labels_pos) + (n_methods * box_width)
    ax.set_xlim(x_min, x_max)

    legend_order = ["rawdata", "DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]
    handles = [Patch(facecolor=palette[m], edgecolor='k') for m in legend_order]
    ax.legend(handles, legend_order, title="Method", bbox_to_anchor=(1.05,1), loc="upper left", prop={'size':8})
    ax.get_legend().get_title().set_fontsize(8)
    
    ax.set_title(plot_title, fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(output, dpi=300)
    plt.close()
    

def plot_fraction_auc(cell_fraction_results, dataset_order, dataset_to_ratio, outdir, descend_flag=True):
    """
    This function is used to plot relative decontamination efficacy based on cell_fraction_results and calculate AUC. 
    AUC, the higher, the better.

    'relative decontamination efficacy' is defined as follows. 
    1. Set cell fraction of rawdata as base
    2. Calculate decontamination efficacy of each method relative to rawdata
    3. The maximum value of decontamination efficacy is set to 1
    4. Other decontamination efficacy is normalized based on maximum value to get 'relative decontamination efficacy'
    
    params:
    cell_fraction_results: cell fraction, with columns 'dataset', 'method', 'value'
    dataset_order: order of x-axis to show
    dataset_to_ratio: x-axis coordinate value, in order to calculate auc

    """
    
    palette = {
        "rawdata": "#BFBFBF",
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    methods_to_plot = [m for m in palette.keys() if m != "rawdata"]

    # Build DataFrame from nested dict
    df = pd.DataFrame([
        (dataset, method, value) 
        for dataset, methods in cell_fraction_results.items()
        for method, value in methods.items()
    ], columns=["Dataset", "Method", "Fraction"])

    # Ensure dataset column is categorical with your defined order
    df["Dataset"] = pd.Categorical(df["Dataset"], categories=dataset_order, ordered=True)

    # Compute decontamination efficacy
    efficacy_records = []
    for dataset in dataset_order:
        raw_val = cell_fraction_results[dataset]["rawdata"]
        for method in methods_to_plot:
            val = cell_fraction_results[dataset][method]
            # eff = abs(raw_val - val) / raw_val * 100  
            eff = (raw_val - val) / raw_val * 100
            efficacy_records.append((dataset, method, eff))

    df_eff = pd.DataFrame(efficacy_records, columns=["Dataset", "Method", "Efficacy"])
    # print(df_eff)

    auc_results = {}
    for method in methods_to_plot:
        # Reindex to ensure correct order
        sub = df_eff[df_eff["Method"] == method].set_index("Dataset").reindex(dataset_order)
        x = [dataset_to_ratio[d] for d in dataset_order]  # x-axis: downsample ratio
        y = sub["Efficacy"].values

        # Normalize y to 100%
        y_scaled = np.array(y / y.max() * 100 if y.max() != 0 else y)
        x = np.array(x)

        # avoid minus auc
        order = np.argsort(x)
        x_sorted = x[order]
        y_sorted = y_scaled[order]

        # Compute AUC
        # auc_val = auc(x, y_scaled / 100)
        # auc_results[method] = auc_val
        
        auc_val = auc(x_sorted, y_sorted / 100)
        auc_results[method] = auc_val

        # Plot
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(x, y_scaled, marker="o", color=palette[method], label=method)
        ax.fill_between(x, 0, y_scaled, color=palette[method], alpha=0.2)

        # Annotate AUC
        x_center = (x[0] + x[-1]) / 2
        y_center = y_scaled.max() / 2
        ax.text(x_center, y_center, f"AUC={auc_val:.2f}", fontsize=8, ha='center', va='center', weight='bold')

        # X-axis labels
        ax.set_xticks(x)
        ax.set_xticklabels([int(v*100) for v in x], fontsize=8)

        if descend_flag:
            ax.set_xlim(max(x)+0.05, min(x)-0.05) # descend
        else:
            ax.set_xlim(min(x)-0.05,max(x)+0.05) # ascend

        ax.set_xlabel("Cell downsample (%)", fontsize=9)
        ax.set_ylabel("Relative efficacy (%)", fontsize=9)
        ax.tick_params(axis='both', labelsize=8)
        ax.set_title(f"{method}", fontsize=10, weight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{method}_fraction_decont_eff_auc.pdf")
        plt.savefig(outpath, dpi=300)
        plt.close()
        print(f"Saved {outpath}")

    return auc_results

def plot_auc_from_efficacy(df_eff, dataset_order, dataset_to_ratio, outdir, descend_flag=True, outdir_auc_stemplot=None):
    """
    This function is used to plot 'relative decontamination efficacy' based on pre-calculated decontamination efficacy and calculate AUC. 
    AUC, the higher, the better.

    note: the another format of 'plot_fraction_auc' function

    params:
    df_eff: pre-calculated decontamination efficacy or other metrics, with columns "Dataset", "Method", "Value"
    dataset_order: order of x-axis to show
    dataset_to_ratio: x-axis coordinate value, in order to calculate auc
    outdir: Each method draw an auc plot into this folder
    descend_flag: The order of x-axis. if true, in descend order. Default True.
    outdir_auc_stemplot: All method draw an auc plot by stemplot
    """
    
    palette = {
        "rawdata": "#BFBFBF",
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    methods_to_plot = [m for m in palette.keys() if m != "rawdata"]

    auc_results = {}
    for method in methods_to_plot:
        # Reindex to ensure correct order
        sub = df_eff[df_eff["Method"] == method].set_index("Dataset").reindex(dataset_order)
        x = [dataset_to_ratio[d] for d in dataset_order]  # x-axis: downsample ratio
        y = sub["Value"].values

        # Normalize y to 100%
        y_scaled = np.array(y / y.max() * 100 if y.max() != 0 else y)
        x = np.array(x)

        # avoid minus auc
        order = np.argsort(x)
        x_sorted = x[order]
        y_sorted = y_scaled[order]

        # Compute AUC
        # auc_val = auc(x, y_scaled / 100)
        # auc_results[method] = auc_val
        
        auc_val = auc(x_sorted, y_sorted / 100)
        auc_results[method] = auc_val

        # Plot
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(x, y_scaled, marker="o", color=palette[method], label=method)
        ax.fill_between(x, 0, y_scaled, color=palette[method], alpha=0.2)

        # Annotate AUC
        x_center = (x[0] + x[-1]) / 2
        y_center = y_scaled.max() / 2
        ax.text(x_center, y_center, f"AUC={auc_val:.2f}", fontsize=8, ha='center', va='center', weight='bold')

        # X-axis labels
        ax.set_xticks(x)
        ax.set_xticklabels([int(v*100) for v in x], fontsize=8)

        if descend_flag:
            ax.set_xlim(max(x)+0.05, min(x)-0.05) # descend
        else:
            ax.set_xlim(min(x)-0.05,max(x)+0.05) # ascend

        ax.set_xlabel("Cell downsample (%)", fontsize=9)
        ax.set_ylabel("Relative efficacy (%)", fontsize=9)
        ax.tick_params(axis='both', labelsize=8)
        ax.set_title(f"{method}", fontsize=10, weight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        plt.tight_layout()
        if outdir is not None:
            outpath = os.path.join(outdir, f"{method}_decont_eff_auc.pdf")
            plt.savefig(outpath, dpi=300)
            plt.close()
            print(f"Saved {outpath}")

    ########### auc stemplot
    plt.figure(figsize=(8, 8))
    methods = methods_to_plot
    values = [auc_results[m] for m in methods]

    # plt.hlines(methods, xmin=0, xmax=values, colors="black")
    # plt.scatter(values, methods, color="#FF8C00")
    for m, val in zip(methods, values):
        color = palette.get(m, "black")  
        plt.hlines(m, xmin=0, xmax=val, colors=color, linewidth=2)  
        plt.scatter(val, m, color=color, s=100, zorder=3)  

    plt.axvline(x=0.5, color="red", linestyle="--", linewidth=1)

    plt.xlabel("AUC")
    plt.ylabel("Method")
    plt.title("AUC comparison of decontamination methods")
    plt.xticks(rotation=45, ha="right")
    plt.gca().invert_yaxis() # invert yaxis
    plt.tight_layout()

    if outdir_auc_stemplot is not None:
        plt.savefig(outdir_auc_stemplot, dpi=300)
    plt.close()

    return auc_results

def save_fraction_results_to_csv(results_dict, filename, long_format=True):
    rows = []
    for dataset, methods in results_dict.items():
        for method, values in methods.items():
            if isinstance(values, (list, tuple)):  # if vector
                for v in values:
                    rows.append({"dataset": dataset, "method": method, "value": v})
            else:  # if single value
                rows.append({"dataset": dataset, "method": method, "value": values})
    
    df = pd.DataFrame(rows)
    if not long_format:
        df = df.pivot(index="dataset", columns="method", values="value").reset_index()
    
    df.to_csv(filename, index=False)
    print(f"Saved {filename}")
    return df


