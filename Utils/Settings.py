from pathlib import Path


## DIRECTORIES

root_data = Path("/alzheimer/Roberto/Allen_Institute/")
root_github_repo = r"/home/roberto/Github/Allen_Institute_Neuropixel"
neuropixel_dataset = Path(f'{root_data}/Visual') # directory Allen dataset
neuropixel_dataset_behavior = Path(f'{root_data}/Visual_Behavior') # directory Allen dataset behavior
output_folder_calculations = Path(f"{root_data}/Processed_data_2023")
output_folder_figures_calculations = Path(f"{root_data}/temp")
output_folder_processed_lfps = Path(f"{root_data}/Processed_lfps")
output_folder = Path(f"{root_github_repo}/Output_figures")
output_folder_supplementary = Path(f"{root_github_repo}/Output_figures/Supplementary")
manuscript_folder = Path(f"{root_github_repo}/Manuscript")
utils_folder = Path(f"{root_github_repo}/Utils")


## ANALYSIS PARAMETERS

fs_lfp = 1250

# bandpass filter

lowcut = 120.0
highcut = 250.0

# ripple detection

ripple_dur_lim = [0.015, 0.250]
min_ripple_distance = 0.05
ripple_power_sp_thr = 100
ripple_thr = 5 #SD

# window for calculation ripple strength RIVD for correlation analysis

start_w = 0.1
stop_w = 0.2

minimum_ripples_count_lag_analysis = 500
minimum_ripples_count_spike_analysis = 200
minimum_ripples_count_generated_in_lateral_or_medial_spike_analysis = 100

thr_rip_cluster = 0.06 # belong to same ripple if closer than, in s

clip_ripples_clusters = (-60, 60) # ms

# threshold to discard channels with weak ripple activity

var_thr = 5

# windows for spike analysis

window_spike_hist = (0.25, 0.25)
window_spike_per_ripple = (0.12, 0.12) #pre post

# inh-exc waveform duration threshold

waveform_dur_thr = 0.4

# parameter for good units selection

waveform_PT_ratio_thr = 5
isi_violations_thr = .5
amplitude_cutoff_thr = .1
presence_ratio_thr = .1

# Nature style figure letters

Adapt_for_Nature_style = False

