import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, median_absolute_error
from scipy.stats import kendalltau, wilcoxon
import os
import statsmodels.api as sm

# 04
from plot_contam_ratio import *

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from brokenaxes import brokenaxes
import seaborn as sns
sns.set_theme(style="white")


def run_all_dataset():
     
    # load 8 new dataset
    csv_files = [
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC1/NC1_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC2/NC2_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC3/NC3_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC4/NC4_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA3/SA3_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA4/SA4_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA6/SA6_merged.csv"
    ]
    method_cols = ["DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
    name_map = {
        "NC1_merged": "NC1",
        "NC2_merged": "NC2",
        "NC3_merged": "NC3",
        "NC4_merged": "NC4",
        "SA3_merged": "SA3",
        "SA4_merged": "SA4",
        "SA5_merged": "SA5",
        "SA6_merged": "SA6",

    }
    data_palette = {
            "NC1": "#4B3C6C",
            "NC2": "#5F5AA5",
            "NC3": "#ADAEEB",
            "NC4": "#C8D2F0",
            "SA3": "#CFDC82",
            "SA4": "#A7C53B",
            "SA5": "#5E891B",
            "SA6": "#283B0A"
    }

    # calculate metric 
    non_osn_df_metrics, osn_df_metrics = merged_contam_ratio_evaluation(csv_files, method_cols)

    non_osn_df_metrics['Dataset'] = non_osn_df_metrics['Dataset'].map(name_map)
    non_osn_df_metrics.to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/all_datasets_non_osn_metrics.csv", index=False)
    osn_df_metrics['Dataset'] = osn_df_metrics['Dataset'].map(name_map)
    osn_df_metrics.to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/all_datasets_osn_metrics.csv", index=False)

    # print(non_osn_df_metrics)

    datasets_to_plot = ["NC1","NC2","NC3","NC4","SA3","SA4","SA5","SA6"]
    # reclaim the method cols
    method_cols = ["ground_truth", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]

    # figure 4C 
    plot_overall_contam_ratio(csv_files, data_palette, name_map, method_cols, datasets_to_plot, outdir="/data/R04/zhangchao/joint_nuclei/figures_new/libraries")
    # figure 4D
    metric_radar_plot(df_metrics=non_osn_df_metrics, metric="RMSLE", datasets_to_plot=datasets_to_plot, outdir="/data/R04/zhangchao/joint_nuclei/figures_new/libraries")
    metric_var_bar(non_osn_df_metrics,
                metric='RMSLE',
                outdir_metric_var='/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_RMSLE_variance.csv',
                output_file='/data/R04/zhangchao/joint_nuclei/figures_new/libraries/non-OSN_metric_RMSLE_var.pdf',
                datasets_to_plot=datasets_to_plot
                )

    # figure 4E
    metric_radar_plot(df_metrics=non_osn_df_metrics, metric="Kendall_tau", datasets_to_plot=datasets_to_plot, outdir="/data/R04/zhangchao/joint_nuclei/figures_new/libraries")
    metric_var_bar(non_osn_df_metrics,
                metric='Kendall_tau',
                outdir_metric_var='/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_kendall_tau_variance.csv',
                output_file='/data/R04/zhangchao/joint_nuclei/figures_new/libraries/non-OSN_metric_Kendall_tau_var.pdf',
                datasets_to_plot=datasets_to_plot
                )


    # figure 5C
    # all 8 datasets are merged for plot
    def plot_lib_merge_osn_violin(csv_path, outpath):
        df = pd.read_csv(csv_path, sep=None)  
        df.columns = df.columns.str.strip()  
        df = df[df['celltype'].str.contains("OSN", case=False)]

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
        # method_order = ["ground_truth", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
        palette = {
            "ground_truth": '#F7FCB9',
            "DecontX": "#9ECAE1", 
            "SoupX": "#4293C6",
            "scAR": "#0052A1",
            "CellBender": "#88419D",
            "FastCAR": "#AE017E",
            "scCDC": "#F768A1",
            "CellClear": "#FDCDE5"
        }
        df_renamed = df.rename(columns=rename_dict)
        df_long = df_renamed.melt(
            id_vars=["celltype", "subtype"],
            value_vars=list(palette.keys()),
            var_name="Method",
            value_name="Fraction"
        )
        df_long["Fraction"] = df_long["Fraction"] * 100
        
        method_medians = (
            df_long.groupby("Method")["Fraction"]
            .median()
            .sort_values(ascending=True)
        )
        method_order = method_medians.index.tolist()

        plt.figure(figsize=(8,8))
        ax = sns.violinplot(
                data=df_long,
                x="Method",
                y="Fraction",
                palette=palette,
                inner="box",
                linewidth=1.2,
                order=method_order
        )
        ax.set_ylabel("Background RNA per cell (%)", fontsize=12)
        plt.xticks(rotation=45, ha="right")
        plt.title("contamination ratio", fontsize=14, weight="bold")
        plt.tight_layout()
        plt.savefig(outpath, dpi=300)
        plt.close()

    plot_lib_merge_osn_violin(csv_path = "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/merged/all_new_dataset_merged.csv", 
                                outpath= "/data/R04/zhangchao/joint_nuclei/figures_new/libraries/eight_lib_OSN_merge_violin.pdf")


    # figure 5D and 5E
    metric_barplot_with_points(osn_df_metrics, datasets_to_plot, 'RMSLE',  '/data/R04/zhangchao/joint_nuclei/figures_new/libraries/OSN_metric_RMSLE.pdf', flag=True)
    metric_barplot_with_points(osn_df_metrics, datasets_to_plot, 'Kendall_tau',  '/data/R04/zhangchao/joint_nuclei/figures_new/libraries/OSN_metric_Kendall_tau.pdf')


run_all_dataset()



#  downsample 
def run_downsample():
    datasets_to_plot=["SA5", "SA5_1200", "SA5_800", "SA5_400"]

    repeats = ["repeat_1", "repeat_2", "repeat_3"]
    cell_num = ["1200", "800", "400"]

    method_cols = ["DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]
    name_map = {
            "SA5_merged": "SA5",
            "SA5_1200_merged": "SA5_1200",
            "SA5_800_merged": "SA5_800",
            "SA5_400_merged": "SA5_400"
    }
    data_palette = {
            "SA5": "#b7282e",
            "SA5_1200": "#c44438",
            "SA5_800": "#d16d5b",
            "SA5_400": "#dc917b"
    }
    root_dir = "/data/R04/zhangchao/joint_nuclei/03_analysis"

    for rep in repeats:
            num_csv_files = [f"{root_dir}/contam_ratio/SA5/SA5_merged.csv"] + [
                f"{root_dir}/SA5_downsampling/SA5_{num}/{rep}/overall_contam/SA5_{num}_merged.csv" for num in cell_num
            ]

            non_osn_df_metrics, osn_df_metrics = merged_contam_ratio_evaluation(num_csv_files, method_cols)

            non_osn_df_metrics['Dataset'] = non_osn_df_metrics['Dataset'].map(name_map)
            osn_df_metrics['Dataset'] = osn_df_metrics['Dataset'].map(name_map)

            output_dir = f"{root_dir}/SA5_downsampling/temp_files/{rep}"
            os.makedirs(output_dir, exist_ok=True)

            non_osn_df_metrics.to_csv(os.path.join(output_dir, "all_datasets_non_osn_metrics.csv"), index=False)
            osn_df_metrics.to_csv(os.path.join(output_dir, "all_datasets_osn_metrics.csv"), index=False)

            print(f"Finished processing {rep}.")

    # nonOSN
            
    all_dfs = []
    for rep in repeats:
        csv_file = os.path.join("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files", rep, "all_datasets_non_osn_metrics.csv")
        df = pd.read_csv(csv_file)
        df["repeat"] = rep  
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, axis=0, ignore_index=True)


    # figure 4H
    mean_metric_lineplot_with_error(df_all, metric="RMSLE", 
                               output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/non_osn_RMSLE.pdf",
                               out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_RMSLE.csv",
                               datasets_to_plot=datasets_to_plot)
    # figure 4I
    mean_metric_lineplot_with_error(df_all, metric="Kendall_tau", 
                               output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/non_osn_Kendall_tau.pdf",
                               out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_Kendall_tau.csv",
                               datasets_to_plot=datasets_to_plot)

    # figure S4B
    repeat_combined_csv_files=[
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
         "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/SA5_1200_combined_rep_contam_ratio.csv",
         "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/SA5_800_combined_rep_contam_ratio.csv",
         "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/SA5_400_combined_rep_contam_ratio.csv"
    ]
    name_map = {
            "SA5_merged": "SA5",
            "SA5_1200_combined_rep_contam_ratio": "SA5_1200",
            "SA5_800_combined_rep_contam_ratio": "SA5_800",
            "SA5_400_combined_rep_contam_ratio": "SA5_400"
    }
    method_cols = ["ground_truth", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
    plot_overall_contam_ratio(repeat_combined_csv_files, data_palette, name_map, method_cols, datasets_to_plot, 
                              outdir="/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/combined_rep")


run_downsample()



# down depth
def run_downdepth():
    datasets_to_plot=["SA5", "SA5_75pct", "SA5_50pct", "SA5_25pct"]

    repeats = ["repeat_1", "repeat_2", "repeat_3"]
    depth_percentages = ["75", "50", "25"]

    method_cols = ["DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]
    name_map = {
            "SA5_merged": "SA5",
            "SA5_75pct_merged": "SA5_75pct",
            "SA5_50pct_merged": "SA5_50pct",
            "SA5_25pct_merged": "SA5_25pct"
    }
    data_palette = {
                "SA5": "#396C9A",
                "SA5_75pct": "#6A92C2",
                "SA5_50pct": "#97B4D9",
                "SA5_25pct": "#C9D8EB"
    }
    root_dir = "/data/R04/zhangchao/joint_nuclei/03_analysis"

    for rep in repeats:
            depth_csv_files = [f"{root_dir}/contam_ratio/SA5/SA5_merged.csv"] + [
                f"{root_dir}/SA5_downdepth/SA5_{pct}_percentage/{rep}/overall_contam/SA5_{pct}pct_merged.csv" for pct in depth_percentages
            ]

            non_osn_df_metrics, osn_df_metrics = merged_contam_ratio_evaluation(depth_csv_files, method_cols)

            non_osn_df_metrics['Dataset'] = non_osn_df_metrics['Dataset'].map(name_map)
            osn_df_metrics['Dataset'] = osn_df_metrics['Dataset'].map(name_map)

            output_dir = f"{root_dir}/SA5_downdepth/temp_files/{rep}"
            os.makedirs(output_dir, exist_ok=True)

            non_osn_df_metrics.to_csv(os.path.join(output_dir, "all_datasets_non_osn_metrics.csv"), index=False)
            osn_df_metrics.to_csv(os.path.join(output_dir, "all_datasets_osn_metrics.csv"), index=False)

            print(f"Finished processing {rep}.")

    # nonOSN
            
    all_dfs = []
    for rep in repeats:
        csv_file = os.path.join("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files", rep, "all_datasets_non_osn_metrics.csv")
        df = pd.read_csv(csv_file)
        df["repeat"] = rep  
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, axis=0, ignore_index=True)

    # figure 4L
    mean_metric_lineplot_with_error(df_all, metric="RMSLE", 
                               output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/non_osn_RMSLE.pdf",
                               out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_RMSLE.csv",
                               datasets_to_plot=datasets_to_plot)
    # figure 4M 
    mean_metric_lineplot_with_error(df_all, metric="Kendall_tau", 
                               output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/non_osn_Kendall_tau.pdf",
                               out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_Kendall_tau.csv",
                               datasets_to_plot=datasets_to_plot)


    # figure S4F
    repeat_combined_csv_files=[
        "/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/SA5_75pct_combined_rep_contam_ratio.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/SA5_50pct_combined_rep_contam_ratio.csv",
        "/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/SA5_25pct_combined_rep_contam_ratio.csv"
    ]
    name_map = {
            "SA5_merged": "SA5",
            "SA5_75pct_combined_rep_contam_ratio": "SA5_75pct",
            "SA5_50pct_combined_rep_contam_ratio": "SA5_50pct",
            "SA5_25pct_combined_rep_contam_ratio": "SA5_25pct"
    }
    method_cols = ["ground_truth", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
    plot_overall_contam_ratio(repeat_combined_csv_files, data_palette, name_map, method_cols, datasets_to_plot, 
                              outdir="/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/combined_rep")
    

run_downdepth()