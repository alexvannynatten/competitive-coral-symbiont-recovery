# competitive-coral-symbiont-recovery

Data and code for: Marine heatwaves transform coral symbioses with enduring effects

Authors: Van Nynatten, A 1; Cunning, R 2; Tietjen, KL 1; Baum, JK 1

1 Department of Biology, University of Victoria, Victoria, British Columbia, Canada 
2 Conservation Research Department, John G. Shedd Aquarium, 1200 South DuSable Lake Shore Drive, Chicago, IL 60605, USA.

**Directory map**
_analyses:_ contains all scripts used to run analyses and generate the figures and tables \
_data:_ all raw and intermediate (labelled according to script used to generate) datafiles used in the analyses \
_figures_tables:_ find output from analyses

**File details**
01_symportal_results_df.csv - Tidied SymPortal ITS2 profile count table output for analysis
02_coralnet_results_df.csv - Tidied CoralNet visual survey data for analysis
04_symbiont_lineages.csv - Symbiont lineage classifications based on Unifrac distances
20250527T170756_unifrac_profile_distances_C_sqrt.dist - SymPortal output for Cladocopium interprofile distances
20250527T170756.profiles.absolute.abund_and_meta.txt - SymPortal output for sequence counts of ITS2 profiles by sample
coral_size.csv - Maximum width (in cm) of tagged Pocillopora grandis colonies
coralnet_labels.csv - Label set for CoralNet classifications
coralnet_raw_5July2024_2013-2023b.csv - CoralNet output for benthic 
dhw_5km_97a5_e43d_e9ef.csv - SST and DHW data from Coral Reef Watch for Kiritimati
KI_shape_files - Kiritimati map shape files
sample_metadata.csv - Metadata for each of the coral tissue samples included in the ITS2 DNA metabarcoding analysis
survey_metadata.csv - Ecological data for each of the survey locations at Kiritimati
note: raw fastq files for ITS2 DNA metabarcoding not included but accession number for each sample can be found in sample_metadata

01_a_runs_symportal_offline.sh  - Runs SymPortal locally
01_b_cleans_symportal_data.R - Tidies symportal output and incorperates sample and survey metadata
02_a_cleans_coralnet_output.R - Decodes coralNet labels, cleans survey data, and incorporates survey metadata
02_b_tests_coral_cover.R - Tests differences in coral cover across expeditions and plots results
03_a_plots_sampling_sites.R - Plots sites surveyed and human disturbance gradient across Kiritimati
03_b_plots_sst_dhw.R - Plots SST and DHW at KI over the period surveyed
03_c_summarizes_samples_collected.R - Summarizes the samples collected and surveys conducted by expedition
04_a_tests_symbiont_diversity.R - Visualizes UniFrac distances and plots profile incidence by expedition
04_b_tests_durusdinium_changes.R - Compares proportion of colonies hosting Durusdinium by expedition
04_c_tests_C15_occupancy.R - Detection probability of C15 based on pre-heatwave sampling
05_a_tests_postheatwave_recovery.R - Tests the recovery of pre-heatwave symbiont lineages in Pocillopora
05_b_tests_tracked_colonies.R - Tests the recovery of pre-heatwave lineages in tracked colonies
05_c_tests_surviving_populations.R - Compares symbiont assemblages in surviving colonies with others
06_a_tests_disturbance_impact.R - Tests the recovery of pre-heatwave lineages across disturbance gradient
06_b_tests_durusdinium_endurance.R - Tests differences in durusdinium hosted before and after heatwave
