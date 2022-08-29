# mTOR
R code for Pokhilko, et.al. Mathematical modelling of brain mTOR activity identifies selective vulnerability of cell types and signalling pathways. 
bioRxiv 2022

To simulate kinetics of mTOR pathway (Fig. S1, Fig. S2), run
mTOR_ode_kinetics_insulin.R; mTOR_ode_kinetics_PDGF.R; mTOR_ode_kinetics_NRG.R

For Fig. 2, Fig. S4,S5 run 
correlation_plots.R

For Fig. 3B-E,run
S6Ka_EC50_plots.R

For Fig. 3A run
dose_response_plots.R

For Fig. 4 run
S6Ka_EC50_ADplots.R

For Fig.5 run
summary_Fig.R

normalized expression values for mTOR components in brain cell types are stored in 
snRNAseq_grouped_upquant_0422.csv; snRNAseq_grouped_upquant_AD_0422.csv (for AD)
these files were generated using mTOR_expression.R (also used for Fig. S3, S6)

results on S6Ka and EC50 are stored in files
S6Ka_max_EC50_all.csv; S6Ka_max_EC50_all_AD.csv (for AD)
these files were generated using S6Ka_EC50_calculate.R
