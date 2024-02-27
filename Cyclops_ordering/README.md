Documentation for the CYCLOPS algorithm can be found <https://github.com/ranafi/CYCLOPS-2.0>

In brief, run_CYCLOPS_script.jl contains the parameters and julia code used for this analysis. Excitatory Neurons (subclusters 3 &5) were used to get an ordering of samples. That ordering (a circadian prediction for each individual) is then assumed for other cell types.

Fits, Plots, Models, Parameters are all subdirs that CYCLOPS creates automatically. Other sub directories present are:

-   *downstream_output\_.\** directories contain downstream results from astrocytes, excitatory neurons, etc. They were generated from the notebook ../Analysis_scripts/perform_downstream_analysis.qmd

-   *proteomics* contains cycling and differential cycling results from TMT_proteomics data from ROSMAP

-   *metabolon* contains cycling and differential cycling results from metabolon data from ROSMAP

-   *PSEA_celltype_diffs* contains the output from ../Analysis_scripts/PSEA_across_celltypes.qmd

-   smoothness.txt and StatErrSummaryFreezeCovs_4EG.csv are quality metrics of the CYCLOPS ordering on this data.

Exc3and5_FiltByEdgeRDefault_fixedipBulkChenZhang_condAndBatchCovs_04ContrVar4EGdefault_noTransferFit
