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

# helper functions
from plot_corrected_mat import *
from plot_contam_ratio import mean_metric_lineplot_with_error

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
mpl.use('Agg') 

import seaborn as sns
sns.set_theme(style="white")

# load OR genes
or_genes = pd.read_csv("/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/groundtruth/merged_OR_genes.csv")
or_genes = list(set(or_genes['V1'].unique()))  # 1140 OR genes
non_osn_genes_of_interest = or_genes + ["Omp"]


# contam_ratio: contam_ratio of all methods in current dataset
# methods: corrected matrix, 'rawdata' means filtered matrix after annotation without any correction

# cellbender_adata = anndata_from_h5("/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/my_cellbender_output_file_filtered.h5")  
# cellbender_adata.write("/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/my_cellbender_output_file_filtered.h5ad")
def run_all_dataset():    
    dataset_configs = {
        "NC1": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC1/NC1_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/NC1_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC1/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/NC1_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/NC1/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/NC1_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_NC1_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/NC1_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/NC1_corrected_mat.h5ad"
            }
        },
        "NC2": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC2/NC2_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/NC2_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC2/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/NC2_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/NC2/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/NC2_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_NC2_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/NC2_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/NC2_corrected_mat.h5ad"
            }
        },
        "NC3": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC3/NC3_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/NC3_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC3/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/NC3_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/NC3/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/NC3_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_NC3_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/NC3_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/NC3_corrected_mat.h5ad"
            }
        },
        "NC4": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/NC4/NC4_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/NC4_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/NC4/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/NC4_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/NC4/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/NC4_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_NC4_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/NC4_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/NC4_corrected_mat.h5ad"
            }
        },
        "SA3": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA3/SA3_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA3_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA3/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA3_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA3/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA3_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA3_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA3_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA3_corrected_mat.h5ad"
            }
        },
        "SA4": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA4/SA4_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA4_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA4/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA4_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA4/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA4_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA4_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA4_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA4_corrected_mat.h5ad"
            }
        },
        "SA5": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA5_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA5/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA5_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA5/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA5_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA5_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA5_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA5_corrected_mat.h5ad"
            }
        },
        "SA6": {
            "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA6/SA6_merged.csv",
            "methods":{
                "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA6_annotated_mat.h5ad",
                "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA6/cellbender_output_file_filtered.h5ad",
                "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA6_corrected_mat.h5ad",
                "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA6/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA6_corrected_mat.h5ad",
                "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA6_matrix/",
                "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA6_corrected_mat.h5ad",
                "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA6_corrected_mat.h5ad"
            }
        }           
    }
    
    # path config method
    methods = ["rawdata", "DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]

    (
        non_osn_cell_fraction_results,
        non_osn_umi_fraction_results,
        osn_cell_fraction_results,
        osn_umi_fraction_results
    ) = calc_fractions_multiple_datasets(dataset_configs, methods, non_osn_genes_of_interest, or_genes)

    # save as pkl
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_cell_fraction.pkl", "wb") as f:
        pickle.dump(non_osn_cell_fraction_results, f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_umi_fraction.pkl", "wb") as f:
        pickle.dump(non_osn_umi_fraction_results, f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_cell_fraction.pkl", "wb") as f:
        pickle.dump(osn_cell_fraction_results, f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_umi_fraction.pkl", "wb") as f:
        pickle.dump(osn_umi_fraction_results, f)
    
    # load
    # with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_cell_fraction.pkl", "rb") as f:
    #     non_osn_cell_fraction = pickle.load(f)

    # save as csv 
    save_fraction_results_to_csv(non_osn_cell_fraction_results, "/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_cell_fraction.csv")
    save_fraction_results_to_csv(non_osn_umi_fraction_results, "/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_umi_fraction.csv")
    save_fraction_results_to_csv(osn_cell_fraction_results, "/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_cell_fraction.csv")
    save_fraction_results_to_csv(osn_umi_fraction_results, "/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_umi_fraction.csv")

# run_all_dataset()
    
def plot_all_dataset():
    # plot 
    # load cell fraction
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_cell_fraction.pkl", "rb") as f:
        non_osn_cell_fraction_results = pickle.load(f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_umi_fraction.pkl", "rb") as f:
        non_osn_umi_fraction_results = pickle.load(f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_cell_fraction.pkl", "rb") as f:
        osn_cell_fraction_results = pickle.load(f)
    with open("/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/osn_umi_fraction.pkl", "rb") as f:
        osn_umi_fraction_results = pickle.load(f)

    # non-osn plot
    datasets_to_plot = ["NC1","NC2","NC3","NC4","SA3","SA4","SA5","SA6"]
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
    method_order = ["DecontX", "SoupX", "scAR", "CellBender", "FastCAR", "scCDC", "CellClear"]
        
    # figure 4F
    plot_decont_efficiency_and_diff(non_osn_cell_fraction_results, 
                                    datasets_to_plot, data_palette=data_palette, 
                                    method_order= method_order,
                                    outdir_df_eff="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_decont_efficacy.csv",
                                    outdir_df_eff_var="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/non_osn_decont_efficacy_variance.csv",
                                    outdir_eff_stripplot=None,
                                    outdir_eff_barplot=None,
                                    outdir_eff_lineplot=None,
                                    outdir_eff_radar="/data/R04/zhangchao/joint_nuclei/figures_new/libraries/non-OSN_decont_efficacy_radar.pdf",
                                    outdir_diff=None,
                                    outdir_var="/data/R04/zhangchao/joint_nuclei/figures_new/libraries/non-OSN_decont_efficacy_variance.pdf",
                                    outdir_eff_meanplot=None
                                    ) 
    # figure 4G  
    dataset_to_ratio = {
        "NC1": 0.125,
        "NC2": 0.25,
        "NC3": 0.375,
        "NC4": 0.5,
        "SA3": 0.625,
        "SA4": 0.75,
        "SA5": 0.875,
        "SA6": 1
    }
    # figure 4G
    plot_fraction_auc(non_osn_cell_fraction_results, datasets_to_plot, dataset_to_ratio, 
                    outdir="/data/R04/zhangchao/joint_nuclei/figures_new/libraries/decont_eff_auc",
                    descend_flag= False
                    )

    # osn plot
    # figure 5F
    plot_decont_efficiency_and_diff(osn_cell_fraction_results,
                                    datasets_to_plot, data_palette=data_palette, 
                                    method_order= method_order,
                                    outdir_df_eff="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/OSN_decont_efficacy.csv",
                                    outdir_df_eff_var="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/OSN_decont_efficacy_variance.csv",
                                    outdir_eff_stripplot=None,
                                    outdir_eff_barplot=None,
                                    outdir_eff_lineplot=None,
                                    outdir_eff_radar=None,
                                    outdir_diff=None,
                                    outdir_var=None,
                                    outdir_eff_meanplot="/data/R04/zhangchao/joint_nuclei/figures_new/libraries/OSN_mean_decont_efficacy_barplot.pdf"
                                    )
        
    # convert umi fraction per cell to mean value in different methods
    osn_umi_fracrion_mean_dict = {
        dataset: {method: np.mean(arr) for method, arr in methods.items()}
        for dataset, methods in osn_umi_fraction_results.items()
    }
    # print(osn_umi_fraction_results['NC1']['rawdata'])
    # print(osn_umi_fracrion_mean_dict)

    # figure 5G
    plot_decont_efficiency_and_diff(osn_umi_fracrion_mean_dict,
                                    datasets_to_plot, data_palette=data_palette, 
                                    method_order= method_order,
                                    outdir_df_eff="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/OSN_umi_decont_efficacy.csv",
                                    outdir_df_eff_var="/data/R04/zhangchao/joint_nuclei/03_analysis/temp_files_new/OSN_umi_decont_efficacy_variance.csv",
                                    outdir_eff_stripplot=None,
                                    outdir_eff_barplot=None,
                                    outdir_eff_lineplot=None,
                                    outdir_eff_radar=None,
                                    outdir_diff=None,
                                    outdir_var=None,
                                    outdir_eff_meanplot="/data/R04/zhangchao/joint_nuclei/figures_new/libraries/OSN_mean_umi_decont_efficacy_barplot.pdf"
                                    )

plot_all_dataset()






###################################################################################################
# SA5 downsampling non-osn into 1200, 800, 400

def run_downsample_dataset():
    repeats = ["repeat_1", "repeat_2", "repeat_3"]

    for rep in repeats:
        print(f"Processing {rep}.")
        
        downsample_dataset_configs = {
            "SA5": {
                "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
                "methods":{
                    "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA5_annotated_mat.h5ad",
                    "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA5/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA5_corrected_mat.h5ad",
                    "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA5/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA5_corrected_mat.h5ad",
                    "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA5_matrix/",
                    "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA5_corrected_mat.h5ad",
                    "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA5_corrected_mat.h5ad"
                }
            },
            "SA5_1200": {
                "contam_ratio": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/overall_contam/SA5_1200_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/scCDC/scCDC_1200_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/decontx/decontx_1200_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_1200/{rep}/soupx/soupx_1200_corrected_mat.h5ad"
                }
            },
            "SA5_800": {
                "contam_ratio": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/overall_contam/SA5_800_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/scCDC/scCDC_800_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/decontx/decontx_800_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_800/{rep}/soupx/soupx_800_corrected_mat.h5ad"
                }
            },
            "SA5_400": {
                "contam_ratio": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/overall_contam/SA5_400_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/scCDC/scCDC_400_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/decontx/decontx_400_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/SA5_400/{rep}/soupx/soupx_400_corrected_mat.h5ad"
                }
            }
        }

        methods = ["rawdata", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
        (
            non_osn_cell_fraction_results,
            non_osn_umi_fraction_results,
            osn_cell_fraction_results,
            osn_umi_fraction_results
        ) = calc_fractions_multiple_datasets(downsample_dataset_configs, methods, non_osn_genes_of_interest, or_genes)

        # save as pickle
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/non_osn_cell_fraction_results.pkl", "wb") as f:
            pickle.dump(non_osn_cell_fraction_results, f)
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/non_osn_umi_fraction_results.pkl", "wb") as f:
            pickle.dump(non_osn_umi_fraction_results, f)

        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/osn_cell_fraction_results.pkl", "wb") as f:
            pickle.dump(osn_cell_fraction_results, f)
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/osn_umi_fraction_results.pkl", "wb") as f:
            pickle.dump(osn_umi_fraction_results, f)

        # load
        # with open("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/osn_cell_fraction_results.pkl", "rb") as f:
        #     test = pickle.load(f)

        # save as csv 
        save_fraction_results_to_csv(non_osn_cell_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/non_osn_cell_fraction.csv")
        save_fraction_results_to_csv(non_osn_umi_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/non_osn_umi_fraction.csv")
        save_fraction_results_to_csv(osn_cell_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/osn_cell_fraction.csv")
        save_fraction_results_to_csv(osn_umi_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/osn_umi_fraction.csv")


        datasets_to_plot = ["SA5","SA5_1200","SA5_800","SA5_400"]
        data_palette = {
                "SA5": "#b7282e",
                "SA5_1200": "#c44438",
                "SA5_800": "#d16d5b",
                "SA5_400": "#dc917b"
        }
        method_order = ["DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]

        # generate decont efficacy for each repeats
        plot_decont_efficiency_and_diff(non_osn_cell_fraction_results,
                                        datasets_to_plot, data_palette=data_palette,
                                        method_order=method_order,
                                        outdir_df_eff= f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/{rep}/non_osn_decont_efficacy.csv",
                                        outdir_df_eff_var=None,
                                        outdir_eff_stripplot=None,
                                        outdir_eff_barplot=None,
                                        outdir_eff_lineplot=None,
                                        outdir_eff_radar=None,
                                        outdir_diff=None,
                                        outdir_var=None,
                                        outdir_eff_meanplot=None
                                        )
        
# run_downsample_dataset()

def plot_downsample_dataset():

    repeats = ["repeat_1", "repeat_2", "repeat_3"]
    datasets_to_plot = ["SA5","SA5_1200","SA5_800","SA5_400"]

    # load all decontamination efficacy and combine
    all_dfs = []
    for rep in repeats:
        csv_file = os.path.join("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files", rep, "non_osn_decont_efficacy.csv")
        df = pd.read_csv(csv_file)
        df["repeat"] = rep  
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, axis=0, ignore_index=True)
    # print(df_all)

    # figure 4J
    mean_metric_lineplot_with_error(df_all, metric="Efficacy", 
                                output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/non_osn_mean_decont_efficacy.pdf",
                                out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/non_osn_mean_decont_efficacy.csv",
                                datasets_to_plot=datasets_to_plot)
        

    dataset_order = ["SA5","SA5_1200","SA5_800","SA5_400"]
    dataset_to_ratio = {
        "SA5": 1,
        "SA5_1200": 0.75,
        "SA5_800": 0.5,
        "SA5_400": 0.25
    }

    # figure 4K and S4E
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/non_osn_mean_decont_efficacy.csv")
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric": "Value"})
    # print(df_mean_eff)
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir="/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/mean_decont_eff_auc", 
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/mean_decont_eff_auc/mean_decont_eff_auc_stemplot.pdf"
                                    )

    # figure S4C
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_RMSLE.csv")
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric": "Value"})
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir=None,
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/mean_RMSLE_auc_stemplot.pdf"
                                    )
    pd.DataFrame.from_dict(auc_results, orient="index").to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_RMSLE_AUC.csv")


    # figure S4D
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_Kendall_tau.csv")
    # min-max normalization to avoid minus AUC
    df_mean_eff["mean_metric_norm"] = df_mean_eff.groupby("Method")["mean_metric"].transform(
        lambda x: (x - x.min()) / (x.max() - x.min()) if (x.max() - x.min()) != 0 else 0
    )
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric_norm": "Value"})
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir=None,
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downsampling/mean_Kendall_tau_auc_stemplot.pdf"
                                    )
    pd.DataFrame.from_dict(auc_results, orient="index").to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downsampling/temp_files/mean_Kendall_tau_AUC.csv")

plot_downsample_dataset()





###################################################################################################
# SA5 down depth into 75%, 50%, 25%

def run_depth_dataset():
    repeats = ["repeat_1", "repeat_2", "repeat_3"]

    for rep in repeats:
        print(f"Processing {rep}.")

        depth_dataset_configs = {
            "SA5": {
                "contam_ratio":"/data/R04/zhangchao/joint_nuclei/03_analysis/contam_ratio/SA5/SA5_merged.csv",
                "methods":{
                    "rawdata": "/data/R04/zhangchao/joint_nuclei/03_analysis/raw/SA5_annotated_mat.h5ad",
                    "CellBender": "/data/R04/zhangchao/joint_nuclei/03_analysis/cellbender/SA5/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   "/data/R04/zhangchao/joint_nuclei/03_analysis/FastCAR/SA5_corrected_mat.h5ad",
                    "scAR":      "/data/R04/zhangchao/joint_nuclei/03_analysis/scar/output/SA5/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC": "/data/R04/zhangchao/joint_nuclei/03_analysis/scCDC/SA5_corrected_mat.h5ad",
                    "CellClear": "/data/R02/tanj93/project/OSN_AmbientRNA_Benchmarking/20_8samples/corrected_rds/cellclear_SA5_matrix/",
                    "DecontX":   "/data/R04/zhangchao/joint_nuclei/03_analysis/decontx/SA5_corrected_mat.h5ad",
                    "SoupX":     "/data/R04/zhangchao/joint_nuclei/03_analysis/soupx/SA5_corrected_mat.h5ad"
                }
            },
            "SA5_75pct": {
                "contam_ratio":f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/overall_contam/SA5_75pct_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/scCDC/scCDC_75pct_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/decontx/decontx_75pct_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_75_percentage/{rep}/soupx/soupx_75pct_corrected_mat.h5ad"
                }
            },
            "SA5_50pct": {
                "contam_ratio":f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/overall_contam/SA5_50pct_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/scCDC/scCDC_50pct_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/decontx/decontx_50pct_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_50_percentage/{rep}/soupx/soupx_50pct_corrected_mat.h5ad"
                }
            },
            "SA5_25pct": {
                "contam_ratio":f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/overall_contam/SA5_25pct_merged.csv",
                "methods":{
                    "rawdata": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/annotated_outs/filtered_feature_bc_matrix.h5ad",
                    "CellBender": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/cellbender/cellbender_output_file_filtered.h5ad",
                    "FastCAR":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/fastcar/fastcar_corrected_mat.h5ad",
                    "scAR":      f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/scar/filtered_feature_bc_matrix_denoised_mRNA.h5ad",
                    "scCDC":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/scCDC/scCDC_25pct_corrected_mat.h5ad",
                    "CellClear": f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/cellclear/matrix/",
                    "DecontX":   f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/decontx/decontx_25pct_corrected_mat.h5ad",
                    "SoupX":     f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/SA5_25_percentage/{rep}/soupx/soupx_25pct_corrected_mat.h5ad"
                }
            }
        }

        methods = ["rawdata", "DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
        (
            non_osn_cell_fraction_results,
            non_osn_umi_fraction_results,
            osn_cell_fraction_results,
            osn_umi_fraction_results
        ) = calc_fractions_multiple_datasets(depth_dataset_configs, methods, non_osn_genes_of_interest, or_genes)

        # save as pickle
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/non_osn_cell_fraction_results.pkl", "wb") as f:
            pickle.dump(non_osn_cell_fraction_results, f)
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/non_osn_umi_fraction_results.pkl", "wb") as f:
            pickle.dump(non_osn_umi_fraction_results, f)

        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/osn_cell_fraction_results.pkl", "wb") as f:
            pickle.dump(osn_cell_fraction_results, f)
        with open(f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/osn_umi_fraction_results.pkl", "wb") as f:
            pickle.dump(osn_umi_fraction_results, f)

        # load
        # with open("/data/R04/zhangchao/joint_nuclei/03_analysis/19d_downdepth/temp_files/non_osn_cell_fraction_results.pkl", "rb") as f:
        #     non_osn_umi_fraction_results_loaded = pickle.load(f)

        # save as csv 
        save_fraction_results_to_csv(non_osn_cell_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/non_osn_cell_fraction.csv")
        save_fraction_results_to_csv(non_osn_umi_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/non_osn_umi_fraction.csv")
        save_fraction_results_to_csv(osn_cell_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/osn_cell_fraction.csv")
        save_fraction_results_to_csv(osn_umi_fraction_results, f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/osn_umi_fraction.csv")

        datasets_to_plot = ["SA5","SA5_75pct","SA5_50pct","SA5_25pct"]
        data_palette = {
                "SA5": "#396C9A",
                "SA5_75pct": "#6A92C2",
                "SA5_50pct": "#97B4D9",
                "SA5_25pct": "#C9D8EB"
        }
        method_order = ["DecontX", "SoupX", "scAR","CellBender","FastCAR","scCDC","CellClear"]
            
        # generate decont efficacy for each repeats
        plot_decont_efficiency_and_diff(non_osn_cell_fraction_results,
                                            datasets_to_plot, data_palette=data_palette,
                                            method_order=method_order,
                                            outdir_df_eff=f"/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/{rep}/non_osn_decont_efficacy.csv",
                                            outdir_df_eff_var=None,
                                            outdir_eff_stripplot=None,
                                            outdir_eff_barplot=None,
                                            outdir_eff_lineplot=None,
                                            outdir_eff_radar=None,
                                            outdir_diff=None,
                                            outdir_var=None,
                                            outdir_eff_meanplot=None
                                        )

# run_depth_dataset()


def plot_depth_dataset():

    repeats = ["repeat_1", "repeat_2", "repeat_3"]
    datasets_to_plot = ["SA5","SA5_75pct","SA5_50pct","SA5_25pct"]

    # load all decontamination efficacy and combine
    all_dfs = []
    for rep in repeats:
        csv_file = os.path.join("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files", rep, "non_osn_decont_efficacy.csv")
        df = pd.read_csv(csv_file)
        df["repeat"] = rep  
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, axis=0, ignore_index=True)
    # print(df_all)

    # figure 4N
    mean_metric_lineplot_with_error(df_all, metric="Efficacy", 
                                output_file="/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/non_osn_mean_decont_efficacy.pdf",
                                out_csv="/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/non_osn_mean_decont_efficacy.csv",
                                datasets_to_plot=datasets_to_plot)
        

    # figure 4O and S4I
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/non_osn_mean_decont_efficacy.csv")
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric": "Value"})
    # print(df_mean_eff)
    dataset_order = ["SA5","SA5_75pct","SA5_50pct","SA5_25pct"]
    dataset_to_ratio = {
            "SA5": 1,
            "SA5_75pct": 0.75,
            "SA5_50pct": 0.5,
            "SA5_25pct": 0.25
    }
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir="/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/mean_decont_eff_auc",
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/mean_decont_eff_auc/mean_decont_eff_auc_stemplot.pdf"
                                    )


    # figure S4G
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_RMSLE.csv")
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric": "Value"})
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir=None,
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/mean_RMSLE_auc_stemplot.pdf"
                                    )
    pd.DataFrame.from_dict(auc_results, orient="index").to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_RMSLE_AUC.csv")

    # figure S4H
    df_mean_eff = pd.read_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_Kendall_tau.csv")
    # min-max normalization to avoid minus AUC
    df_mean_eff["mean_metric_norm"] = df_mean_eff.groupby("Method")["mean_metric"].transform(
        lambda x: (x - x.min()) / (x.max() - x.min()) if (x.max() - x.min()) != 0 else 0
    )
    df_mean_eff = df_mean_eff.rename(columns={"mean_metric_norm": "Value"})
    auc_results = plot_auc_from_efficacy(df_mean_eff, dataset_order, dataset_to_ratio, 
                                    outdir=None,
                                    descend_flag=True,
                                    outdir_auc_stemplot= "/data/R04/zhangchao/joint_nuclei/figures_new/downdepth/mean_Kendall_tau_auc_stemplot.pdf"
                                    )
    pd.DataFrame.from_dict(auc_results, orient="index").to_csv("/data/R04/zhangchao/joint_nuclei/03_analysis/SA5_downdepth/temp_files/mean_Kendall_tau_AUC.csv")

plot_depth_dataset()







