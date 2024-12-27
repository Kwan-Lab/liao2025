import os

analysis_name = "top12_no11"
fig_path = os.path.join("/nfs/turbo/umms-kykwan/projects/alex_kwan/paper_figures/figures", analysis_name)
os.makedirs(fig_path, exist_ok=True)
os.makedirs(os.path.join(fig_path,"vectors"), exist_ok=True)

data_path = "/nfs/turbo/umms-kykwan/projects/alex_kwan/figures/scanpy_cellranger_fixed_sex"
