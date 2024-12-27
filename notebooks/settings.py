import os
fig_path = "../figures/scanpy_cellranger_fixed_sex"
os.makedirs(fig_path, exist_ok=True)

# Use the subset 250K or the entire thing?
use_subset = True

# Drop the "bad" samples?
drop_samples = True


# Use cellranger or kallisto?
cellranger = True