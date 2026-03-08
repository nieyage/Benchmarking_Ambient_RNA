import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, median_absolute_error
from scipy.stats import kendalltau, wilcoxon, mannwhitneyu
import os
import statsmodels.api as sm

import plotly.express as px

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from brokenaxes import brokenaxes
import seaborn as sns
sns.set_theme(style="white")

def load_cellnames(contam_ratio_path):
    df = pd.read_csv(contam_ratio_path, index_col=0)
    non_osn_cells = df.index[df['celltype'] == 'nonOSN'].tolist()
    osn_cells = df.index[df['celltype'] == 'OSN'].tolist()
    return non_osn_cells, osn_cells

def regression_metrics(y_true, y_pred):
    y_true = np.array(y_true).ravel()
    y_pred = np.array(y_pred).ravel()

    # Concordance Correlation Coefficient, CCC, [-1,1]
    # cov = np.mean((y_true - np.mean(y_true)) * (y_pred - np.mean(y_pred)))
    # ccc = (2 * cov) / (np.var(y_true) + np.var(y_pred) + (np.mean(y_true) - np.mean(y_pred)) ** 2)

    mean_x = np.mean(y_true)
    mean_y = np.mean(y_pred)
    var_x = np.var(y_true, ddof=1)  
    var_y = np.var(y_pred, ddof=1)
    # Pearson correlation
    rho = np.corrcoef(y_true, y_pred)[0,1]
    ccc = (2 * rho * np.sqrt(var_x) * np.sqrt(var_y)) / \
                (var_x + var_y + (mean_x - mean_y)**2)

    # Coefficient of Determination (0,1)
    # r2 = r2_score(y_true, y_pred)

    X = np.array(y_pred).ravel()
    X = sm.add_constant(X)  # intercept
    model = sm.OLS(y_true, X).fit()
    r2 = model.rsquared

    # Kendall's tau [-1,1]
    tau, p_value = kendalltau(y_true, y_pred)

    return {"CCC": ccc, "R2": r2, "Kendall_tau": tau}

def error_metrics(y_true, y_pred):
    y_true = np.array(y_true).ravel()
    y_pred = np.array(y_pred).ravel()

    # SMAPE [0, 2]
    denominator = (np.abs(y_true) + np.abs(y_pred)) / 2.0
    mask = denominator != 0
    if np.any(mask):
        smape = np.mean(np.abs(y_true[mask] - y_pred[mask]) / denominator[mask])
    else:
        smape = np.nan

    # RMSLE [0,+∞)
    rmsle = np.sqrt(np.mean((np.log1p(y_pred) - np.log1p(y_true)) ** 2))

    # MedAE [0,+∞)
    medae = median_absolute_error(y_true, y_pred)

    return {"SMAPE": smape, "RMSLE": rmsle, "MedAE": medae}

def merged_contam_ratio_evaluation(csv_files, method_cols, min_samples=5): 
    """
    csv files: all contamination ratio of all methods, the colname of truth is 'contam_OR'
    """
    osn_results = []
    non_osn_results = []
    for f in csv_files:
        df = pd.read_csv(f)
        dataset_name = os.path.splitext(os.path.basename(f))[0]
        rename_dict = {
            "cellbender": "CellBender",
            "fastcar": "FastCAR",
            "scar": "scAR",
            "scCDC": "scCDC",
            "cellclear": "CellClear",
            "decontx": "DecontX",
            "soupx": "SoupX"
        }
        df.rename(columns=rename_dict, inplace=True)

        osn_df = df[df["celltype"] == "OSN"].copy()
        non_osn_df = df[df["celltype"] == "nonOSN"].copy()

        for celltype, df_sub, res_list in [("OSN", osn_df, osn_results), ("nonOSN", non_osn_df, non_osn_results)]:
            print(f"Processing dataset {dataset_name}, celltype {celltype}...")
            y_true_all = df_sub['contam_OR']

            for method in method_cols:
                y_pred_all = df_sub[method]

                # remove NA
                mask = y_true_all.notna() & y_pred_all.notna()
                y_true = y_true_all[mask].values
                y_pred = y_pred_all[mask].values

                if len(y_true) < min_samples:
                    print(f"skip {dataset_name} - {celltype} - {method} (valid sample {len(y_true)})")
                    continue
                try:
                    r_metrics = regression_metrics(y_true, y_pred)
                    e_metrics = error_metrics(y_true, y_pred)

                    res_list.append({
                        "Dataset": dataset_name,
                        "Celltype": celltype,   
                        "Method": method,
                        **r_metrics,
                        **e_metrics
                    })
                except Exception as e:
                    print(f"computation failure {dataset_name} - {celltype} - {method}: {e}")

    return (pd.DataFrame(non_osn_results), pd.DataFrame(osn_results))

def metric_lineplot(df_metrics, metric, output_file, datasets_to_plot=None, figsize=(8,8)):
    """
    df_metrics: DataFrame, columns include Dataset, Method, metric 
    metric: str, 'CCC' / 'R2' / 'Kendall_tau'
    """
    plt.figure(figsize=figsize)
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
        plot_data = df_metrics[df_metrics["Dataset"].isin(datasets_to_plot)].copy()
        plot_data["Dataset"] = pd.Categorical(plot_data["Dataset"], categories=datasets_to_plot, ordered=True)
    else:
        plot_data = df_metrics.copy()

    ax = sns.lineplot(
        data=plot_data,
        x="Dataset",
        y=metric,
        hue="Method",
        palette=palette,
        marker="o",
        markersize=12,
        linewidth=1.2
    )

    ax.set_title(f"Comparison of {metric} across methods", fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis='both', labelsize=8)
    
    ymin, ymax = ax.get_ylim()
    if ymin != 0:
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8) 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.legend(title='Method', bbox_to_anchor=(1.05, 1), loc='upper left',prop = {'size':8})
    plt.gca().get_legend().get_title().set_fontsize(8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


def mean_metric_lineplot_with_error(df_all, metric, output_file, out_csv, datasets_to_plot=None, figsize=(8,8)):
    """
    df_all: DataFrame, columns include Dataset, Method, metric
    metric: str, e.g., 'Kendall_tau'
    output_file: str
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

    if datasets_to_plot is not None:
        df_all = df_all[df_all["Dataset"].isin(datasets_to_plot)].copy()
        df_all["Dataset"] = pd.Categorical(df_all["Dataset"], categories=datasets_to_plot, ordered=True)

    df_summary = df_all.groupby(["Dataset", "Method"], as_index=False)[metric].agg(
        mean_metric="mean", std_metric="std"
    )
    # mean results
    df_summary.to_csv(out_csv)

    plt.figure(figsize=figsize)
    ax = sns.lineplot(
        data=df_summary,
        x="Dataset",
        y="mean_metric",
        hue="Method",
        palette=palette,
        marker="o",
        markersize=10,
        linewidth=1.8,
        markeredgewidth=0,
        legend=False
    )

    for method in df_summary["Method"].unique():
        method_data = df_summary[df_summary["Method"] == method]
        plt.errorbar(
            x=method_data["Dataset"],
            y=method_data["mean_metric"],
            yerr=method_data["std_metric"],
            fmt="none",
            ecolor=palette[method],
            elinewidth=1.8, 
            capsize=5, 
            alpha=0.8
        )

    ax.set_title(f"Comparison of {metric} across methods", fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel(metric)
    ax.tick_params(axis='both', labelsize=9)

    ymin, ymax = ax.get_ylim()

    # 0 axisline
    if ymin > 0:
        ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300)
    plt.close()


def metric_barplot_with_points(df_metrics, datasets_to_plot, metric, output_file, flag=False):
    df_metrics = df_metrics[df_metrics["Dataset"].isin(datasets_to_plot)].copy()

    palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }

    mean_vals = df_metrics.groupby("Method")[metric].mean()
    std_vals = df_metrics.groupby("Method")[metric].std()

    # flag means descending or ascending by values, False = descending
    methods = mean_vals.sort_values(ascending=flag).index.tolist()
    mean_vals = mean_vals.reindex(methods)
    std_vals = std_vals.reindex(methods)
    # print(mean_vals)

    x = np.arange(len(methods))*0.8  # bar gap width
    fig, ax = plt.subplots(figsize=(4,4))
    bar_width = 0.4


    for i, method in enumerate(methods):
        ax.bar(
            x[i],
            mean_vals[method],
            yerr=std_vals[method],
            color=palette.get(method, "lightgray"),
            width=bar_width,
            capsize=4
        )
        subset = df_metrics[df_metrics["Method"] == method]
        y_vals = subset[metric].values
        # add points representing datasets
        ax.scatter(np.full_like(y_vals, x[i]), y_vals,
                   color="black", s=10, zorder=3)
    
    ax.set_xticks(x)
    ax.set_xticklabels(methods, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("")
    ax.tick_params(axis="y", labelsize=8)
    
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    
    ax.set_title(f"{metric}", fontsize=8)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

def metric_flower_plot(df_metrics, metric, datasets_to_plot, outdir):
    """
    Draws a "flower plot" for a specific metric across two datasets. 
    Each dataset is represented as a flower: points for different methods
    are arranged around the mean to form petals. The plot also performs
    a paired Wilcoxon test between the two datasets and annotates the p-value.

    Parameters:
    -----------
    df_metrics : pd.DataFrame
        DataFrame containing metrics with columns: ['Dataset', 'Method', metric, ...]
    metric : str
        The column name of the metric to plot (e.g., 'Kendall_tau')
    datasets_to_plot : list of str
        e.g. List of two dataset names to include in the plot
    outdir : str
        Output directory to save the figure
    """

    # Define color palette for each method  
    method_palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
     # Filter data to include only selected datasets
    plot_data = df_metrics[df_metrics["Dataset"].isin(datasets_to_plot)].copy()
    print(plot_data)

    # Map datasets to x-axis positions
    dataset_to_x = {ds: i for i, ds in enumerate(datasets_to_plot)}
    radius = 0.15  # radius for flower spread (x jitter)

    plt.figure(figsize=(4, 4))
    ax = plt.gca()

    # Draw each dataset as a flower
    for ds in datasets_to_plot:
        sub = plot_data[plot_data["Dataset"] == ds]
        x_center = dataset_to_x[ds]
        n_points = len(sub)

        # Evenly distribute methods in a circle around dataset center — only in x
        angles = np.linspace(0, 2*np.pi, n_points, endpoint=False)
        x = x_center + np.cos(angles) * radius  # jitter x
        y = sub[metric].values                  # keep true metric

        # Plot each point (petal)
        for xi, yi, m in zip(x, y, sub["Method"]):
            ax.scatter(xi, yi, color=method_palette.get(m, "gray"), s=70, alpha=0.8)

        # Plot mean value as black dot
        mean_val = sub[metric].mean()
        ax.scatter(x_center, mean_val, color="black", s=120, zorder=5, label="Mean" if ds == datasets_to_plot[0] else "")

    # Set axis ticks and labels
    ax.set_xticks(list(dataset_to_x.values()))
    ax.set_xticklabels(datasets_to_plot)
    ax.set_ylabel(metric)
    ax.set_ylim(-1,1)
    plt.xlim(-0.5, len(datasets_to_plot) - 0.5)

    # Wilcoxon test
    if len(datasets_to_plot) == 2:
        ds1, ds2 = datasets_to_plot
        vals1 = plot_data[plot_data["Dataset"] == ds1][metric].values
        vals2 = plot_data[plot_data["Dataset"] == ds2][metric].values
        stat, p_value = wilcoxon(vals1, vals2)
        y_max = max(max(vals1), max(vals2)) + 0.05
        ax.plot([dataset_to_x[ds1], dataset_to_x[ds2]], [y_max+0.2, y_max+0.2], color='k', lw=1.2)
        ax.text((dataset_to_x[ds1] + dataset_to_x[ds2]) / 2, y_max + 0.3, f"p = {p_value:.2e}", ha='center')

    sns.despine(ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{metric}_flower.pdf"), dpi=300, bbox_inches="tight")
    plt.close()


def metric_radar_plot(df_metrics, metric, datasets_to_plot, outdir):

    # Define color palette for each method  
    method_palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }

    # filter data
    df_plot = df_metrics[(df_metrics["Dataset"].isin(datasets_to_plot))].copy()

    # it requires that 'Dataset' column of 'df_plot' are sorted by dataset name.
    fig = px.line_polar(df_plot,
                        r=metric,
                        theta="Dataset",
                        color="Method",
                        line_close=True,
                        color_discrete_map=method_palette,
                        markers=True
                        )
    
    # remove background color
    fig.update_layout(
        title=f"Radar plot of {metric} for nonOSN",  # radar plot title
        paper_bgcolor="white",   # whole figure background
        plot_bgcolor="white",    # plot area background
        font=dict(color="black"),  # set the color of axis label, legend and so on as black
        polar=dict(
            bgcolor="white",              # radar plot background
            radialaxis=dict(gridcolor="lightgrey"),  # radial grid circles
            angularaxis=dict(gridcolor="lightgrey")  # angular grid lines
        )
    )
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"non-OSN_{metric}_radar.pdf")

    fig.write_image(outpath)




def metric_diff_bar(df_metrics, metric, datasets_to_plot, outdir):
    palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }
    # Pivot table: rows=Method, columns=Dataset, values=metric
    plot_data = df_metrics[df_metrics["Dataset"].isin(datasets_to_plot)].copy()
    ds1, ds2 = datasets_to_plot

    # Only keep methods that exist in both datasets
    pivot = plot_data.pivot_table(index="Method", columns="Dataset", values=metric)
    pivot["Diff"] = (pivot[ds1] - pivot[ds2]).abs() # Calculate abs Δmetric
    pivot = pivot.sort_values("Diff", ascending=False) 

    plt.figure(figsize=(8,8))
    bars = plt.bar(
        pivot.index,
        pivot["Diff"],
        color = [palette.get(m,"gray") for m in pivot.index]
    )
    
    plt.title(f"The difference of Kendall tau across libraries")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    fname = f"{metric}_diff_{ds1}_vs_{ds2}"
    plt.savefig(os.path.join(outdir, f"{fname}.pdf"))
    plt.close()

def metric_var_bar(df_metrics, metric, outdir_metric_var, output_file, datasets_to_plot=None, figsize=(8,8)):
    palette = {
        "DecontX": "#9ECAE1", 
        "SoupX": "#4293C6",
        "scAR": "#0052A1",
        "CellBender": "#88419D",
        "FastCAR": "#AE017E",
        "scCDC": "#F768A1",
        "CellClear": "#FDCDE5"
    }

    if datasets_to_plot is not None:
        df_plot = df_metrics[df_metrics["Dataset"].isin(datasets_to_plot)].copy()
    else:
        df_plot = df_metrics.copy()
    
    var_data = df_plot.groupby("Method")[metric].var().reset_index()
    var_data = var_data.sort_values(by=metric) 

    # print(var_data)
    var_data.to_csv(outdir_metric_var)

    plt.figure(figsize=figsize)
    # plt.bar(
    #     var_data["Method"], 
    #     var_data[metric], 
    #     color=[palette[m] for m in var_data["Method"]]
    # )

    sns.barplot(
        data=var_data,
        x="Method",
        y=metric,
        palette=palette,    
        order=var_data["Method"]   
    )
    
    plt.xlabel("Method")
    plt.xticks(rotation=45, ha="right")
    plt.title(f"Variance of {metric} across datasets")
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()



def plot_overall_contam_ratio(csv_files, palette, dataset_labels, method_cols, datasets_to_plot, outdir):
    """
    This funtion is used to plot overall contamination ratio across different using various tools
    
    csv_files: overall comtamination ratio table, raw represents barcode, column represents different tools
    palette: color set of bar
    dataset_label: dataset name needed to show
    method_cols: different tools' name
    datasets_to_plot: dataset that needs to plot
    outdir: output file direction
    """
    # rename the csv column names to standard name
    rename_dict = {
        "cellbender": "CellBender",
        "fastcar": "FastCAR",
        "scar": "scAR",
        "scCDC": "scCDC",
        "cellclear": "CellClear",
        "decontx": "DecontX",
        "soupx": "SoupX",
        "contam_OR": "ground_truth"
    }

    all_data = []
    for file in csv_files:
        df = pd.read_csv(file)
        dataset_name = os.path.basename(file).replace(".csv", "")
        df = df.rename(columns=rename_dict)
        methods = [c for c in df.columns if c in method_cols]
        
        df_melt = df.melt(
            id_vars="celltype",
            value_vars=methods,
            var_name="Method",
            value_name="Contamination"
        )
        df_melt["Dataset"] = dataset_name
        all_data.append(df_melt)

    df_all = pd.concat(all_data, ignore_index=True)
    df_all["Method"] = pd.Categorical(df_all["Method"], categories=method_cols, ordered=True)
    df_all["Dataset_label"] = df_all["Dataset"].map(dataset_labels)
    
    # convert ratio decimals to percentage
    df_all["Contamination"] = df_all["Contamination"] * 100

    # df_all.to_csv("/data/R04/zhangchao/joint_nuclei/07_figures/all_results.csv", index=False)


    def plot_all_broken_box(data, celltype):
        plot_data = data[data["celltype"] == celltype]
        low_range = (0, 2)    
        high_range = (2, 100) 

        fig, (ax_high, ax_low) = plt.subplots(
            2, 1, sharex=True, figsize=(8, 4),
            gridspec_kw={'height_ratios':[4,1], 'hspace':0.05}  
        )
        # high
        ax_high = sns.boxplot(
            data=plot_data,
            x="Method", y="Contamination",
            hue="Dataset_label",
            palette=palette,
            showfliers=False,
            order=method_cols,
            width=0.6,
            ax=ax_high
        )
        # low
        ax_low = sns.boxplot(
            data=plot_data,
            x="Method", y="Contamination",
            hue="Dataset_label",
            palette=palette,
            showfliers=False,
            order=method_cols,
            width=0.6,
            ax=ax_low
        )
        ax_high.set_xlabel("")
        ax_high.set_ylim(high_range)
        ax_high.set_ylabel("")
        ax_high.get_legend().remove()

        ax_low.set_xlabel("")
        ax_low.set_yticks(np.arange(low_range[0], low_range[1]+0.01, 0.5))
        ax_low.set_ylim(low_range)
        ax_low.set_ylabel("")
        ax_low.tick_params(top=False)
        ax_low.get_legend().remove()

        handles, labels = ax_low.get_legend_handles_labels() 
        ax_high.legend(handles, labels, title="Dataset", loc='upper left', bbox_to_anchor=(1.02,1), borderaxespad=0,  title_fontsize=8, prop = {'size':8})

        fig.text(0.05, 0.55,"Contamination percentage(%)", va="center", rotation="vertical")
        fig.suptitle(f"Comparision of {celltype} cells' contamination ratio in different datasets across methods")

        # add brokenaxis mark line
        d = .015
        kwargs = dict(transform=ax_low.transAxes, color='k', clip_on=False)
        ax_low.plot((-d, +d),(1-d,1+d), **kwargs)
        ax_low.plot((1-d,1+d),(1-d,1+d), **kwargs)
        kwargs.update(transform=ax_high.transAxes)
        ax_high.plot((-d,+d),(-d,+d), **kwargs)
        ax_high.plot((1-d,1+d),(-d,+d), **kwargs)

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{celltype}_box.pdf")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"saved {outpath} successfully.")

    def plot_all_broken_bar(data, celltype):
        plot_data = data[data["celltype"] == celltype].copy()
        mean_data = plot_data.groupby(["Method", "Dataset_label"], sort=False)["Contamination"].mean().reset_index()

        datasets_order = datasets_to_plot

        low_range = (0, 3)
        high_range = (3, 90)

        method_cols = ["DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
        x = np.arange(len(method_cols))  # bar gap across method

        fig, (ax_high, ax_low) = plt.subplots(
            2, 1, sharex=True, 
            figsize=(16,8),
            gridspec_kw={'height_ratios':[1,1], 'hspace':0.05}
        )
        
        width =  0.8 / len(datasets_order) # bar width
        for i, dataset in enumerate(datasets_order):
            subset = mean_data[mean_data["Dataset_label"] == dataset]
            heights = [subset[subset["Method"]==m]["Contamination"].values[0] for m in method_cols]
            # offsets = x + i*width - width*(len(datasets_order)-1)/2
            offsets = x - 0.8/2 + i*width + width/2
            ax_high.bar(offsets, heights, width=width, color=palette[dataset], label=dataset)
            ax_low.bar(offsets, heights, width=width, color=palette[dataset], label=dataset)
        
        ax_high.set_ylim(high_range)
        ax_high.tick_params(axis='y', labelsize=8)
        ax_high.set_xticks([])  
        ax_high.set_xlabel('')  
        ax_high.spines['bottom'].set_visible(False)

        ax_low.set_xticks(x)
        ax_low.set_xticklabels(method_cols, fontsize=8)
        ax_low.set_yticks(np.arange(low_range[0], low_range[1]+0.01, 0.5))
        ax_low.set_ylim(low_range)
        ax_low.tick_params(axis='y', labelsize=8)

        for ax in [ax_high, ax_low]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
        
        fig.text(0.05, 0.55,"Percentage (%)", va="center", rotation="vertical", fontsize=8)
        handles, labels = ax_low.get_legend_handles_labels()
        fig.legend(
            handles, labels, title="Dataset", 
            loc='upper center', ncol=len(datasets_order),  
            bbox_to_anchor=(0.5, 0.93),  
            title_fontsize=8, prop={'size':8}
        )

        outpath = os.path.join(outdir, f"{celltype}_bar.pdf")
        plt.tight_layout()
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"saved {outpath} successfully.")

    def plot_violin_omp_or(data, celltype):
        plot_data = data[(data["celltype"] == celltype) & (data["Method"] == "ground_truth") & data["Dataset_label"].isin(datasets_to_plot)].copy()
        plot_data = plot_data.dropna(subset=["Contamination"]) # remove na values
        plot_data["Method"] = plot_data["Method"].astype(str)
        # print(plot_data)

        median_order = (
            plot_data.groupby("Dataset_label")["Contamination"]
            .median()
            .sort_values(ascending=True)
            .index.tolist()
        )

        mean_order = (
            plot_data.groupby("Dataset_label")["Contamination"]
            .mean()
            .sort_values(ascending=True)
            .index.tolist()
        )        
        
        fig, ax = plt.subplots(figsize=(8,8))
        sns.violinplot(data=plot_data, 
                       x="Dataset_label", y="Contamination", hue="Dataset_label", 
                       palette=palette, split=False, 
                       order=datasets_to_plot,   # Sort by data input order
                    #    order=median_order,     # Sort by median
                    #    order=mean_order,       # Sort by mean
                       inner_kws=dict(box_width=10, whis_width=1, color=".9"), # color represents white
                       linewidth=1.2, ax=ax)

        ax.set_ylabel("Contamination (%)", fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        if celltype=="nonOSN":
            ax.set_ylim(0,0.2)
        else:
            ax.set_ylim(0,10)
        
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        outpath = os.path.join(outdir, f"{celltype}_omp_or_violin.pdf")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_violin_omp_or_test(data, celltype):
        plot_data = data[
            (data["celltype"] == celltype) &
            (data["Method"] == "ground_truth") &
            data["Dataset_label"].isin(datasets_to_plot)
        ].copy()
        plot_data = plot_data.dropna(subset=["Contamination"])
        
        group1 = ["NC1","NC2","NC3","NC4"]
        group2 = ["SA3","SA4","SA5","SA6"]

        # test
        data1 = plot_data[plot_data["Dataset_label"].isin(group1)]["Contamination"]
        data2 = plot_data[plot_data["Dataset_label"].isin(group2)]["Contamination"]
        # H1: data1 < data2 （NC less than SA）
        stat, pvalue = mannwhitneyu(data1, data2, alternative="less")
        print(f"{celltype} NC vs SA p-value:", pvalue)

        plot_order = group1 + [""] + group2
        plot_data_temp = plot_data.copy()
        plot_data_temp.loc[~plot_data_temp["Dataset_label"].isin(group1+group2), "Dataset_label"] = ""

        fig, ax = plt.subplots(figsize=(8,8))
        sns.violinplot(
            data=plot_data_temp,
            x="Dataset_label", y="Contamination", hue="Dataset_label",
            palette=palette,
            split=False,
            order=plot_order,
            inner="box",         
            linewidth=1.2,
            inner_kws=dict(box_width=10, whis_width=1, color=".9"),
            cut=0, # Trim sharp corners
            ax=ax
        )

        ax.set_ylabel("Contamination (%)", fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        y_max = np.percentile(plot_data["Contamination"], 95)
        if celltype=="nonOSN":
            ax.set_ylim(0, min(y_max, 0.2))
            # ax.set_ylim(0, 0.4)
        else:
            ax.set_ylim(0, min(y_max, 10))

        if pvalue < 0.001:
            star = "***"
        elif pvalue < 0.01:
            star = "**"
        elif pvalue < 0.05:
            star = "*"
        else:
            star = "ns"
        
        y_min, y_max = ax.get_ylim()
        x1 = plot_order.index(group1[0])
        x2 = plot_order.index(group2[-1])
        y = y_max * 1.05
        h = y_max * 0.03
        ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='black')
        ax.text((x1 + x2) / 2, y + h, star, ha='center', va='bottom', fontsize=12)

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{celltype}_omp_or_test_violin.pdf")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()
    
    def plot_box_omp_or_test(data, celltype):
        plot_data = data[
            (data["celltype"] == celltype) &
            (data["Method"] == "ground_truth") &
            data["Dataset_label"].isin(datasets_to_plot)
        ].copy()
        plot_data = plot_data.dropna(subset=["Contamination"])
        
        group1 = ["NC1", "NC2", "NC3", "NC4"]
        group2 = ["SA3", "SA4", "SA5", "SA6"]

        data1 = plot_data[plot_data["Dataset_label"].isin(group1)]["Contamination"]
        data2 = plot_data[plot_data["Dataset_label"].isin(group2)]["Contamination"]
        stat, pvalue = mannwhitneyu(data1, data2, alternative="less")
        print(f"{celltype} NC vs SA p-value:", pvalue)

        plot_order = group1 + [""] + group2
        plot_data_temp = plot_data.copy()
        plot_data_temp.loc[
            ~plot_data_temp["Dataset_label"].isin(group1 + group2), 
            "Dataset_label"
        ] = ""

        fig, ax = plt.subplots(figsize=(8, 8))
        sns.boxplot(
            data=plot_data_temp,
            x="Dataset_label",
            y="Contamination",
            order=plot_order,
            palette=palette,  
            fliersize=3,
            linewidth=1.2,
            showfliers=False,
            ax=ax
        )

        ax.set_ylabel("Contamination (%)", fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.xticks(rotation=45, ha="right")

        if pvalue < 0.001:
            star = "***"
        elif pvalue < 0.01:
            star = "**"
        elif pvalue < 0.05:
            star = "*"
        else:
            star = "ns"

        ymin, ymax = ax.get_ylim()
        y = ymax * 0.92
        h = ymax * 0.04

        x1 = plot_order.index(group1[0])
        x2 = plot_order.index(group2[-1])

        ax.plot([x1, x1, x2, x2],
                [y, y + h, y + h, y],
                lw=1.5, c='black')

        ax.text((x1 + x2) / 2, y + h * 1.2,
                star, ha='center', va='bottom',
                fontsize=12, fontweight='bold')

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{celltype}_omp_or_test_boxplot.pdf")
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close()

    def plot_split_voilin(data, celltype):
        # Filter data for the given celltype, exclude "OMP" method, and only include selected datasets
        plot_data = data[(data["celltype"] == celltype) & (data["Method"] != "ground_truth") & data["Dataset_label"].isin(datasets_to_plot)].copy()
        plot_data = plot_data.dropna(subset=["Contamination"]) # remove na values
        plot_data["Method"] = plot_data["Method"].astype(str)

        method_order = ["DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]

        plt.figure(figsize=(12,8))
        # Draw split violin plot: x-axis = Method, y-axis = Contamination, hue = Dataset_label
        sns.violinplot(data=plot_data, x="Method", y="Contamination", hue="Dataset_label", split=True, 
               density_norm = 'count', 
               order= method_order,
               inner_kws=dict(box_width=10, whis_width=1, color=".9"),
               fill=True, palette=palette)
        
        # Get current axis
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False) 
        ax.set_ylim(bottom=0)

        plt.tight_layout()
        outpath = os.path.join(outdir, f"{celltype}_omp_or_split_violin.pdf")
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close()


    for celltype in ["nonOSN", "OSN"]:
        # plot_all_broken_box(df_all, celltype)
        plot_all_broken_bar(df_all, celltype)
        plot_violin_omp_or(df_all, celltype)
        # plot_violin_omp_or_test(df_all,celltype)
        # plot_box_omp_or_test(df_all,celltype)
        # plot_split_voilin(df_all, celltype)





