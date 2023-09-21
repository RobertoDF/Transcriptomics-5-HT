import math
import random
from scipy.stats import pearsonr as pearsonr

from IPython.display import Markdown
import itertools
from joblib import Parallel, delayed
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import SelectKBest, mutual_info_classif
from sklearn.model_selection import train_test_split
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
from sklearn.model_selection import cross_val_score

from sklearn.preprocessing import StandardScaler

import random
from sklearn.pipeline import Pipeline
from matplotlib.patches import Rectangle, Patch
from allensdk.brain_observatory.behavior.behavior_project_cache.behavior_neuropixels_project_cache import VisualBehaviorNeuropixelsProjectCache
from Utils.Settings import output_folder_calculations, neuropixel_dataset_behavior
from collections import namedtuple
import pandas as pd
import numpy as np

from tqdm.notebook import tqdm
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from scipy.stats import zscore

import pickle
from multiprocessing import Pool, Manager
pd.set_option('display.max_columns', None)

from itertools import compress
from itertools import groupby
from operator import itemgetter

from Utils.Settings import neuropixel_dataset, utils_folder
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import scipy

from collections import namedtuple
from collections import ChainMap
from Utils.Settings import waveform_dur_thr, window_spike_per_ripple, var_thr, output_folder_calculations, ripple_dur_lim, min_ripple_distance, ripple_power_sp_thr, ripple_thr, minimum_ripples_count_lag_analysis, thr_rip_cluster
from matplotlib.lines import Line2D
import os

from sklearn.feature_selection import SelectKBest, mutual_info_classif
from scipy.stats import zscore
from multiprocessing import Pool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib import animation
from matplotlib.colors import rgb2hex
from matplotlib.patches import Ellipse
from numpy import sum, count_nonzero
from scipy import signal
from scipy.signal import butter, sosfiltfilt
from scipy.signal import hilbert, medfilt, detrend
from scipy.stats import skew, pearsonr
from tqdm import tqdm
from tqdm.notebook import tqdm as tqdm_notebook
from Utils.Settings import lowcut, highcut, ripple_dur_lim

try:

    import sklearn
except:
    pass

try:
    from allensdk.brain_observatory.ecephys.ecephys_project_cache import EcephysProjectCache
    from allensdk.api.queries.ontologies_api import OntologiesApi
    import nrrd
    from allensdk.core.structure_tree import StructureTree
    from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
    from allensdk.core.reference_space_cache import ReferenceSpaceCache, ReferenceSpace
    from allensdk.api.queries.mouse_connectivity_api import MouseConnectivityApi
    from allensdk.config.manifest import Manifest

    # define acronym map
    oapi = OntologiesApi()
    structure_graph = oapi.get_structures_with_sets([1])  # 1 is the id of the adult mouse structure graph
    # This removes some unused fields returned by the query
    structure_graph = StructureTree.clean_structures(structure_graph)
    tree = StructureTree(structure_graph)
    mcc = MouseConnectivityCache()
    structure_tree = mcc.get_structure_tree()
    # get the ids of all the structure sets in the tree
    structure_set_ids = structure_tree.get_structure_sets()

    summary_structures_all = pd.DataFrame(tree.nodes())
    summary_structures = pd.DataFrame(structure_tree.get_structures_by_set_id(
        [2]))  # Structures representing the major divisions of the mouse brain
    summary_structures = pd.concat([summary_structures, pd.DataFrame.from_dict(structure_tree.get_structures_by_acronym(["VS"]))], axis=0)

    summary_structures_finer = pd.DataFrame(structure_tree.get_structures_by_set_id(
        [3]))
    summary_structures_HP = pd.DataFrame(structure_tree.get_structures_by_set_id(
        [688152359]))


    # useful
    # pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
    # pd.DataFrame(tree.nodes())

    acronym_color_map = tree.value_map(lambda x: x['acronym'], lambda y: y['rgb_triplet'])
    acronym_structure_path_map = tree.value_map(lambda x: x['acronym'], lambda y: y['structure_id_path'])
    acronym_structure_id_map = tree.value_map(lambda x: x['acronym'], lambda y: y['id'])
    acronym_structure_graph_order_map = tree.value_map(lambda x: x['acronym'], lambda y: y['graph_order'])
    id_acronym_map = tree.value_map(lambda x: x['id'], lambda y: y['acronym'])

    with open(f"{utils_folder}/summary_structures.pkl", "wb") as fp:
        pickle.dump(summary_structures, fp)

except:
    with open(f"{utils_folder}/acronym_color_map.pkl", "rb") as fp:  # Unpickling
        acronym_color_map = pickle.load(fp)
    with open(f"{utils_folder}/acronym_structure_path_map.pkl", "rb") as fp:  # Unpickling
        acronym_structure_path_map = pickle.load(fp)
    with open(f"{utils_folder}/acronym_structure_id_map.pkl", "rb") as fp:  # Unpickling
        acronym_structure_id_map = pickle.load(fp)
    with open(f"{utils_folder}/acronym_structure_graph_order_map.pkl", "rb") as fp:  # Unpickling
        acronym_structure_graph_order_map = pickle.load(fp)
    with open(f"{utils_folder}/id_acronym_map.pkl", "rb") as fp:  # Unpickling
        id_acronym_map = pickle.load(fp)
    with open(f"{utils_folder}/summary_structures.pkl", "rb") as fp:  # Unpickling
        summary_structures = pickle.load(fp)


def pick_acronym(name, pos):
    # pos is zero, 1 if 2-APN
    s = name.split("-", 2)
    if len(s) == 3:
        name = s[1]
    elif len(s) == 2:
        name = s[pos]
    else:
        name = s[0]
    return name


def color_to_labels(axs, which_axes, minor_or_major, pos=1, *args, **kwargs):
    if which_axes == 'y':

        if minor_or_major == 'minor':
            for ytick in axs.get_yticklabels(minor='True'):
                if ytick.get_text() != 'nan' and ytick.get_text() != []:
                    name = ytick.get_text()
                    name = pick_acronym(name, pos)
                    color = rgb2hex([x / 255 for x in acronym_color_map.get(name)])
                else:
                    color = [0, 0, 0]
                ytick.set_color(color)
        else:
            for ytick in axs.get_yticklabels():
                if ytick.get_text() != 'nan' and ytick.get_text() != []:
                    name = ytick.get_text()
                    name = pick_acronym(name, pos)
                    color = rgb2hex([x / 255 for x in acronym_color_map.get(name)])
                else:
                    color = [0, 0, 0]
                ytick.set_color(color)
    else:
        if minor_or_major == 'minor':
            for xtick in axs.get_xticklabels(minor='True'):
                if xtick.get_text() != 'nan' and xtick.get_text() != []:
                    name = xtick.get_text()
                    name = pick_acronym(name,pos)
                    color = rgb2hex([x / 255 for x in acronym_color_map.get(name)])
                else:
                    color = [0, 0, 0]
                xtick.set_color(color)
        else:
            for xtick in axs.get_xticklabels():
                if xtick.get_text() != 'nan' and xtick.get_text() != []:
                    name = xtick.get_text()
                    name = pick_acronym(name,pos)
                    color = rgb2hex([x / 255 for x in acronym_color_map.get(name)])
                else:
                    color = [0, 0, 0]
                xtick.set_color(color)



def acronym_to_graph_order(area):
    # return (index according to summary_structures(major divisions of mouse brain), and acronym of parent structure)

    if area == "grey":
        graph_order = 2000
    else:
        structure_path = acronym_structure_path_map.get(area)
        parent = set(structure_path).intersection(summary_structures['id']).pop()
        graph_order = summary_structures[summary_structures['id'] == parent]['graph_order'].values[0]

    return graph_order


def acronym_to_main_area(area):
    # return (index according to summary_structures(major divisions of mouse brain), and acronym of parent structure)

    if area == "grey":
        acronym = "grey"
    elif area == "root":
        acronym = "root"
    else:
        structure_path = acronym_structure_path_map.get(area)
        parent = set(structure_path).intersection(summary_structures['id']).pop()
        acronym = summary_structures[summary_structures['id'] == parent]['acronym'].values[0]

    return acronym

def acronym_to_main_area_finer(area):
    # return (index according to summary_structures(major divisions of mouse brain), and acronym of parent structure)

    structure_path = acronym_structure_path_map.get(area)

    if area == "grey":
        acronym = "grey"
    elif len(set(structure_path).intersection(summary_structures_finer['id'])) == 0:
        #print(f"{area} already too generic")
        acronym = area

    else:
        parent = set(structure_path).intersection(summary_structures_finer['id']).pop()
        acronym = summary_structures_finer[summary_structures_finer['id'] == parent]['acronym'].values[0]

    return acronym

def find_ripples_clusters_new(ripples, source_area):

    out = []

    ripples_source_area = ripples[ripples["Probe number-area"] == source_area]

    for ripple_source in tqdm(ripples_source_area.iterrows(), total=ripples_source_area.shape[0]):

        near_ripples = ripples.loc[np.abs(ripples["Start (s)"] - ripple_source[1]["Start (s)"]) < thr_rip_cluster].copy()
        near_ripples["Start (s)"] = near_ripples["Start (s)"] - ripple_source[1]["Start (s)"]

        # if detected twice on source probe we keep the source (0.00 s)
        if near_ripples[near_ripples["Probe number-area"] == source_area].shape[0] > 1:
            to_drop = near_ripples[
                (near_ripples["Probe number-area"] == source_area) & (near_ripples["Start (s)"] != 0)].index
            near_ripples = near_ripples.drop(to_drop)

        # if detected twice or more on one probe we keep the first one temporally
        if np.any(near_ripples["Probe number-area"].duplicated(keep=False)):
            to_drop = near_ripples.index[near_ripples["Probe number-area"].duplicated(keep=False)]
            to_drop = to_drop.drop(
                near_ripples[near_ripples["Probe number-area"].duplicated(keep=False)]["Start (s)"].idxmin())
            near_ripples = near_ripples.drop(to_drop)

        # print(near_ripples)

        positions = near_ripples[["D-V (µm)", "A-P (µm)", "M-L (µm)"]].copy()
        distances = positions.apply(lambda row: calculate_distance(
            near_ripples[near_ripples["Start (s)"] == 0][["D-V (µm)", "A-P (µm)", "M-L (µm)"]], row[["D-V (µm)", "A-P (µm)", "M-L (µm)"]]), axis=1)

        velocities = distances / (near_ripples["Start (s)"].copy().abs() * 1000)

        strengths = near_ripples["Z-scored ∫Ripple"].copy()

        velocities.index = near_ripples["Probe number-area"].copy() + " speed (µm/ms)"
        distances.index = near_ripples["Probe number-area"].copy() + " distance (µm)"
        strengths.index = near_ripples["Probe number-area"].copy() + " Z-scored ∫Ripple"

        ripple_numbers = near_ripples["Ripple number"].copy()
        ripple_numbers.index = near_ripples["Probe number-area"] + " ripple number"

        _ = pd.Series(near_ripples["Start (s)"]).copy()
        _.name = ripple_source[1]["Ripple number"]

        _.index = near_ripples["Probe number-area"] + " lag (s)"
        _["∫Ripple"] = ripple_source[1]["∫Ripple"]
        _["Duration (s)"] = ripple_source[1]["Duration (s)"]
        _["Start (s)"] = ripple_source[1]["Start (s)"]
        _["Strongest ripple M-L (µm)"] = near_ripples.loc[near_ripples["Z-scored ∫Ripple"].idxmax()]["M-L (µm)"]
        if ripples["Probe number-area"].unique().shape[0]>1:
            _["Spatial engagement"] = (near_ripples.shape[0]-1) / (ripples["Probe number-area"].unique().shape[0]-1)
        else:
            _["Spatial engagement"] = np.nan

        _["Source M-L (µm)"] = near_ripples.loc[near_ripples["Start (s)"].idxmin()]["M-L (µm)"]

        #if near_ripples.shape[0] > 1:
        # # linear regression to determine direction
        # y = near_ripples["M-L (µm)"]
        # x = near_ripples["Start (s)"]
        # fit = np.polyfit(x, y, deg=1)
        # # predict = np.poly1d(fit)
        # # plt.plot(x, predict(x), color="k", alpha=0.3)
        # # plt.scatter(x, y)
        #
        # if fit[0] > 0:
        #     _["Direction"] = "M\u2192L"
        # elif fit[0] < 0:
        #     _["Direction"] = "L\u2192M"

        # strenght across M-L axis
        if near_ripples.loc[near_ripples["Start (s)"].idxmin()]["Probe number-area"] == source_area:
            _["Global strength"] = near_ripples["Local strong"].sum()/ripples["Probe number-area"].unique().shape[0]
        else:
            _["Global strength"] = np.nan

        if near_ripples.loc[near_ripples["Start (s)"].idxmin()]["Probe number-area"] == source_area:
            _["Source"] = True
        else:
            _["Source"] = False
        _2 = pd.concat([_, velocities, strengths, distances, ripple_numbers])

        out.append(_2)
    _ = pd.concat(out, axis=1).T
    _.index.name = "Ripple number"
    _.columns.name = ""

    return _.reindex(sorted(_.columns), axis=1)

def calculate_distance(x, y):
    return np.linalg.norm(x-y)

def batch_trajectories(ripples_calcs, kind, func, sessions):
    """
    spatial_info[0]= relative to source/seed
    spatial_info[1]= relative to all others CA1 locations
    """

    print(f"var_thr: {var_thr}")
    input_rip = []
    for session_id in ripples_calcs.keys():
        ripples = ripples_calcs[session_id][3].copy()
        ripples = ripples[ripples["Area"] == "CA1"]
        ripples = ripples.groupby("Probe number-area").filter(lambda group: group["∫Ripple"].var() > var_thr)
        input_rip.append(ripples.groupby("Probe number-area")["M-L (µm)"].mean())

    ml_space = pd.concat(input_rip)

    medial_lim = ml_space.quantile(.33333)
    lateral_lim = ml_space.quantile(.666666)
    center = ml_space.median()

    input_rip = []
    spatial_infos = []

    for session_id in tqdm(ripples_calcs.keys()):
        ripples = ripples_calcs[session_id][3].copy()

        sel_probe = ripples_calcs[session_id][5]

        if ripples[ripples['Probe number'] == sel_probe].shape[0] < 1000:
            continue

        ripples = ripples.sort_values(by="Start (s)").reset_index(drop=True)

        ripples = ripples[ripples["Area"] == "CA1"]

        ripples = ripples.groupby("Probe number-area").filter(lambda group: group["∫Ripple"].var() > var_thr)

        if ripples["Probe number-area"].unique().shape[0] < 2:
            continue

        ripples["Local strong"] = ripples.groupby("Probe number").apply(
            lambda x: x["∫Ripple"] > x["∫Ripple"].quantile(.9)).sort_index(level=1).values

        ripples["Session"] = session_id

        ripples = ripples.reset_index().rename(columns={'index': 'Ripple number'})

        ripples["Z-scored ∫Ripple"] = ripples.groupby("Probe number-area").apply(
            lambda group: zscore(group["∫Ripple"], ddof=1)).droplevel(0)

        pos_matrix = ripples.groupby("Probe number-area")[["D-V (µm)", "A-P (µm)", "M-L (µm)"]].mean().T

        distance_matrix = pos_matrix.corr(method=calculate_distance)
        distance_tab = distance_matrix.where(np.triu(np.ones(distance_matrix.shape)).astype(bool))
        distance_tab.columns.name = None
        distance_tab.index.name = None
        distance_tab = distance_tab.stack().reset_index()
        distance_tab.columns = ['Reference area', 'Secondary area', 'Distance (µm)']

        full_genotype = sessions.loc[session_id]["full_genotype"]

        if kind == "medial":
            position = "Medial"
            source_area = str(ripples["Probe number"].min()) + "-CA1"
            sub_distance_tab = distance_tab[distance_tab["Reference area"]==source_area]
            if (np.any(ripples.groupby("Probe number-area")["M-L (µm)"].mean() < medial_lim)) & (np.any(ripples.groupby("Probe number-area")["M-L (µm)"].mean() > medial_lim)): # check that there are recording sites in other hippocampal sections
                spatial_infos.append([ripples.groupby("Probe number-area").mean().loc[source_area, :], pd.DataFrame(
                    ripples.groupby("Probe number-area").mean().loc[
                    ripples.groupby("Probe number-area").mean().index != source_area, :])])
                input_rip.append((session_id, ripples, sub_distance_tab, source_area, full_genotype, position))
        elif kind == "lateral":
            position = "Lateral"
            source_area = str(ripples["Probe number"].max()) + "-CA1"
            sub_distance_tab = distance_tab[(distance_tab["Reference area"] == source_area) | (distance_tab["Secondary area"] == source_area)]
            sub_distance_tab.columns = ["Secondary area", "Reference area", "Distance (µm)"]
            if (np.any(ripples.groupby("Probe number-area")["M-L (µm)"].mean() < lateral_lim)) & (np.any(ripples.groupby("Probe number-area")["M-L (µm)"].mean() > lateral_lim)):
                spatial_infos.append([ripples.groupby("Probe number-area").mean().loc[source_area, :], pd.DataFrame(
                    ripples.groupby("Probe number-area").mean().loc[
                    ripples.groupby("Probe number-area").mean().index != source_area, :])])
                input_rip.append((session_id, ripples, sub_distance_tab, source_area, full_genotype, position))
        elif kind == "center":
            position = "Center"
            if (np.any(ripples.groupby("Probe number-area")["M-L (µm)"].mean().between(medial_lim, lateral_lim))) & (np.any(~ripples.groupby("Probe number-area")["M-L (µm)"].mean().between(medial_lim, lateral_lim))):
                source_area = ripples.groupby("Probe number-area")["M-L (µm)"].mean().sub(center).abs().idxmin()
                sub_distance_tab = distance_tab[distance_tab["Reference area"]==source_area]
                spatial_infos.append([ripples.groupby("Probe number-area").mean().loc[source_area, :], pd.DataFrame(
                    ripples.groupby("Probe number-area").mean().loc[
                    ripples.groupby("Probe number-area").mean().index != source_area, :])])
                input_rip.append((session_id, ripples, sub_distance_tab, source_area, full_genotype, position))

    with Pool(processes=len(input_rip)) as pool:
        r = pool.starmap_async(func, input_rip)
        list_propagation_chances = r.get()
        pool.terminate()
    trajectories = pd.concat(list_propagation_chances).reset_index(drop=True)

    return trajectories, spatial_infos

def add_distance_col(row, sub_distance_tab):
    return sub_distance_tab[sub_distance_tab["Secondary area"]==row.index]["Distance (µm)"]

def get_propagation_percentage_by_distance(session_id, ripples, sub_distance_tab, source_area, full_genotype, position):

    real_ripple_summary = find_ripples_clusters_new(ripples, source_area)
    propagation_table = real_ripple_summary.loc[:, real_ripple_summary.columns.str.contains('lag')].count()/real_ripple_summary.loc[:, real_ripple_summary.columns.str.contains('lag')].count().max()*100

    propagation_table.index = [w[:5] for w in propagation_table.index]
    propagation_table = pd.DataFrame(propagation_table)
    propagation_table.columns = ["Ripples propagated to (%)"]
    print(session_id)
    propagation_table["Distance (µm)"] = propagation_table.apply(add_distance_col, sub_distance_tab=sub_distance_tab).values
    propagation_table["Session id"] = session_id
    propagation_table["Genotype"] = full_genotype
    propagation_table["Reference location"] = position

    return propagation_table

def batch_process_spike_clusters_per_ripple(func, real_ripple_summary, units, spike_times, target_area,
                                            field_to_use_to_compare, n_cpus, window):

    if units[units[field_to_use_to_compare] == target_area].shape[0] > 0:

        space_sub_spike_times = dict(zip(units[units[field_to_use_to_compare] == target_area].index,
                                         itemgetter(*units[units[field_to_use_to_compare] == target_area].index)(
                                             spike_times)))

        input_process_spike_per_ripple = []
        ripples_features = []
        for index, row in real_ripple_summary.iterrows():
            time_center = row["Start (s)"] + row[row.index.str.contains('lag')].min() # either sum zero or sum a negative value
            ripples_features.append(row[["Start (s)", "Duration (s)", "∫Ripple", "Location seed", "Global strength", "Spatial engagement", "Source M-L (µm)", "Source"]])
            input_process_spike_per_ripple.append((time_center, space_sub_spike_times, window))

        with Pool(processes=n_cpus) as pool:
            r = pool.starmap_async(func, input_process_spike_per_ripple, chunksize=250)
            time_space_sub_spike_times = r.get()
            pool.close()

    return time_space_sub_spike_times, ripples_features, units

def process_spikes_per_ripple(time_center, space_sub_spike_times, window):

    time_space_sub_spike_times = {
        cluster_id: spikes[(spikes > time_center - window[0]) & (spikes < time_center + window[1])] - time_center for
        cluster_id, spikes in space_sub_spike_times.items()}


    return time_space_sub_spike_times


def ripple_finder(sig, fs, threshold_ripples, probe_n, area):

    filtered = butter_bandpass_filter(np.nan_to_num(sig.values), lowcut, highcut, fs, order=6)
    analytic_signal = hilbert(filtered)
    amplitude_envelope = np.abs(analytic_signal)
    amplitude_envelope = pd.Series(amplitude_envelope, index=sig.time)

    if acronym_to_graph_order(area) == 454:  # if area belongs to HPF
        print("HPF specific threshold")
        threshold = np.std(amplitude_envelope) * ripple_thr
    else:
        print("generic threshold")
        threshold = threshold_ripples

    unwrapped = np.unwrap(np.angle(analytic_signal))
    instantaneous_frequency = np.diff(medfilt(unwrapped, int(fs / 1000 * 12)))  # according to tingsley buzsaky 2021
    instantaneous_frequency = instantaneous_frequency / (2 * np.pi) * fs

    peaks, peak_properties = signal.find_peaks(amplitude_envelope, height=threshold, distance=fs/1000*50)

    peak_heights = peak_properties["peak_heights"]

    prominences, left_bases, right_bases = signal.peak_prominences(amplitude_envelope, peaks)
    results_half = signal.peak_widths(amplitude_envelope.rolling(5, center=True).mean(), peaks, rel_height=1, prominence_data=(
        peak_heights - np.std(amplitude_envelope) * 2, left_bases, right_bases))  # width evaluated at: np.std(amplitude_envelope) * 2
    peaks_width_sec = results_half[0] / fs
    mask = (peaks_width_sec > ripple_dur_lim[0]) & (peaks_width_sec < ripple_dur_lim[1])  # (peaks_width_sec > 0.015)   & (peaks_width_sec < 0.250) Tingley
    peak_heights = peak_heights[mask]

    peaks_sec = sig.time.values[peaks][mask]
    peaks_start_sec = sig.time.values[np.rint(results_half[2]).astype(int)][mask]
    peaks_stop_sec = sig.time.values[np.rint(results_half[3]).astype(int)][mask]

    array_diff = np.diff(peaks_start_sec) > min_ripple_distance # if nearer than x keep first one, we check on the start times because it assures peaks are more distant, opposite is not true!
    if array_diff.size != 0:
        array_diff[-1] = True
        start_mask = clean_start_time_ripples(array_diff)
        peak_mask = clean_peak_time_ripples(array_diff, peak_heights)
        peaks_start_sec = peaks_start_sec[start_mask]
        stop_mask = list(array_diff)
        stop_mask.append(True)
        peaks_stop_sec = peaks_stop_sec[stop_mask]
        temp_peak1 = peaks_sec[stop_mask]
        temp_peak2 = peaks_sec[peak_mask]
        peaks_sec = np.maximum(temp_peak1, temp_peak2)
        temp_height1 = peak_heights[stop_mask]
        temp_height2 = peak_heights[peak_mask]
        peak_heights = np.maximum(temp_height1, temp_height2)

        ripples = pd.DataFrame([peaks_sec, peaks_start_sec, peaks_stop_sec,  peak_heights],
                               index=["Peak (s)", "Start (s)", "Stop (s)", "Amplitude (mV)"]).T

        ripples["Probe number"] = probe_n
        ripples["Area"] = area
        ripples["Duration (s)"] = ripples["Stop (s)"] - ripples["Start (s)"]

        instantaneous_frequency_pd = pd.DataFrame(instantaneous_frequency, index=sig.time[:-1].values, columns=["freq"])

        inst_freq = []
        for index, row in ripples.iterrows():
            try:
                inst_freq.append(np.mean(instantaneous_frequency_pd[row["Start (s)"]: row["Stop (s)"]].values))
            except:
                inst_freq.append(np.nan)

        ripples["Instantaneous Frequency (Hz)"] = inst_freq
        ripples["Probe number-area"] = ripples[["Probe number", "Area"]].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)

        amplitude_envelope = pd.Series(amplitude_envelope, index=sig.time)
        #filtered = pd.Series( filtered, index=sig.time)
        sig = sig.to_series().fillna(value=0)

        high_freq_integral = []
        power_peak = []
        freq_peak = []

        ripples_start = ripples["Start (s)"]
        ripples_stop = ripples["Stop (s)"]

        for start, stop in zip(ripples_start, ripples_stop):

            f, pxx = signal.periodogram(sig[(sig.index > start) & (sig.index < stop)], fs=fs, detrend="constant", scaling="spectrum")
            peaks, properties = signal.find_peaks(pxx, prominence=[0.001, None])
            if len(peaks) > 1:
                peaks = peaks.max()
            elif len(peaks) == 1:
                peaks = peaks[0]

            freq_peak.append(f[peaks])
            power_peak.append(pxx[peaks]*1000)
            high_freq_integral.append(np.trapz(amplitude_envelope[(amplitude_envelope.index > start) & (amplitude_envelope.index < stop)]))

        freq_peak = [el if isinstance(el, (int, float)) else np.nan for el in
                     freq_peak]  # don't want empty lists in the series
        power_peak = [el if isinstance(el, (int, float)) else np.nan for el in power_peak]

        ripples["∫Ripple"] = high_freq_integral
        ripples["Peak frequency (Hz)"] = freq_peak
        ripples["Peak power"] = power_peak

        ripples = ripples[~ripples["Start (s)"].duplicated(keep="last")]
        ripples = ripples[~ripples["Stop (s)"].duplicated()]

        # filter by peak freq
        print(probe_n, area, "Ripples retained by peak freq: ", ripples[ripples["Peak frequency (Hz)"] > 100].shape[0],
              ", total: ", ripples.shape[0])
        ripples = ripples[ripples["Peak frequency (Hz)"] > ripple_power_sp_thr]

        # ripples = ripples.convert_dtypes() # this messes up things somehow

        print("Duplicated starts: " + str(
            np.sum(ripples["Start (s)"].duplicated(keep="last"))) + ", Duplicated stops: " + str(
            np.sum(ripples["Stop (s)"].duplicated())))

        print(probe_n, area, "Ripples discarded: ", sum(~mask))# not counting joined ones

        print(probe_n, area, "Ripples detected: ", len(peaks_sec))
    else:
        ripples = pd.DataFrame()
        print(probe_n, area, "no ripples detected!")

    return ripples


def clean_start_time_ripples(array_diff):
    # if closer than x keep first starting time
    start_mask = []
    flag_next = False
    for key, group in groupby(array_diff):
        # print(key, list(group))
        g = list(group)
        if key == False:
            g[0] = True
            flag_next = True
            start_mask.extend(g)
        elif flag_next == True:
            g[0] = False
            start_mask.extend(g)
            flag_next = False
        else:
            start_mask.extend(g)

    start_mask.append(True)

    return start_mask


def clean_peak_time_ripples(array_diff, peak_heights):
    #pick highest peak if ripples closer than x
    flag_next = False
    peak_mask = []
    for key, group in groupby(zip(array_diff, peak_heights), itemgetter(0)):
        g = list(list(zip(*group))[1])
        gg = ([key] * len(g))
        if key == False:
            gg[np.argmax(g)] = True
            flag_next = True
            peak_mask.extend(gg)
        elif flag_next == True:
            gg[0] = False
            peak_mask.extend(gg)
            flag_next = False
        else:
            peak_mask.extend(gg)

    peak_mask.append(True)
    return peak_mask


def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], btype='band', output="sos")
    return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sosfiltfilt(sos, data)
    return y

def select_quiet_part(lfp, start_quiet, stop_quiet):
    out = []
    for start, stop in zip(start_quiet, stop_quiet):
        out.append(lfp.sel(time=slice(start, stop)))
    out = xr.concat(out, dim="time")
    return out


def calculations_per_ripple(start, stop, sliced_lfp_per_probe, lowcut, highcut, fs_lfp, length):

    high_freq_area = []
    areas = []
    probe_n_area = []
    AUC_pos = []
    AUC_neg = []

    for probe_n, lfp in enumerate(sliced_lfp_per_probe):

        for sig in lfp.T:
            if sig.shape[0] == length:
                sig.values = np.nan_to_num(sig.values)
                filtered = butter_bandpass_filter(sig.values, lowcut, highcut, fs_lfp, order=6)
                analytic_signal = hilbert(filtered)
                amplitude_envelope = np.abs(analytic_signal)
                amplitude_envelope = pd.Series(amplitude_envelope, index=sig.time)
                high_freq_area.append(np.trapz(amplitude_envelope
                                               [(amplitude_envelope.index > start - 0.02) & (amplitude_envelope.index < stop + 0.02)]))
                # for AUC
                sig = sig.sel(time=slice(start, start + 0.2))
                AUC_pos.append(np.trapz(sig.values[sig.values > 0]))
                AUC_neg.append(np.trapz(sig.values[sig.values < 0]))
            else:
                #print(f"LFP from {sig.area.values} at {timestamp} is {sig.shape[0]} instead of {length} samples long.")
                high_freq_area.append(np.nan)
                AUC_pos.append(np.nan)
                AUC_neg.append(np.nan)

        areas.append(lfp.area.values)

        for area in lfp.area.values:
            probe_n_area.append((probe_n, area))

    res = pd.DataFrame([high_freq_area, AUC_pos, AUC_neg],
                       columns=pd.MultiIndex.from_tuples(probe_n_area, names=["probe_n", "area"]),
                       index=["Ripple area (mV*s)", "Positive area (mV*s)", "Negative area (mV*s)"])
    sorting_df = pd.DataFrame(
        [np.concatenate(areas), [acronym_to_graph_order(area.split("-")[0]) for area in np.concatenate(areas)]],
        index=["area", "graph_id"]).T
    idx = sorting_df.sort_values(["graph_id", "area"]).index
    res = res.iloc[:, idx]

    return res, start


def clean_ripples_calculations(ripples_calcs):
    """
    clean from mistakes: one case with negative dv,ap and lr coords, "grey" (unassigned) areas and ripples during
    run.
    """

    for session_id, session_data in tqdm(ripples_calcs.items()):
        print(f"Session {session_id}")

        pos = session_data[1]
        sel_probe = session_data[5]

        out = []
        for t in session_data[0][0].columns:
            out.append(pos[(pos["Probe number"] == t[0]) & (pos["Area"] == t[1])])

        pos = pd.concat(out).reset_index(drop=True)

        mask = (pos["D-V (µm)"].values > 0) & (session_data[0][0].columns.get_level_values(1) != "grey")
        ripples_calcs[session_id][0][0] = session_data[0][0].loc[:, session_data[0][0].columns[mask]]
        ripples_calcs[session_id][0][1] = session_data[0][1].loc[:, session_data[0][1].columns[mask]]
        ripples_calcs[session_id][0][2] = session_data[0][2].loc[:, session_data[0][2].columns[mask]]
        # create M-L axis, subtract midpoint
        pos["M-L (µm)"] = pos["L-R (µm)"] - 5691.510009765625
        ripples_calcs[session_id][1] = pos.iloc[pos.index[mask]]

        # clean from running epochs

        start_running = ripples_calcs[session_id][4][1]
        stop_running = ripples_calcs[session_id][4][2]
        behavior = ripples_calcs[session_id][4][0]

        ripple_power = ripples_calcs[session_id][0][0]
        pos_area = ripples_calcs[session_id][0][1]
        neg_area = ripples_calcs[session_id][0][2]

        #  clean ripples if behavior data not available
        for q in range(len(start_running)):
            idxs = ripple_power.loc[
                   ripple_power.index[(ripple_power.index > start_running[q]) & (ripple_power.index < stop_running[q])],
                   :].index
            if idxs.shape[0] > 0:
                pos_area.drop(idxs, inplace=True)
                neg_area.drop(idxs, inplace=True)
                ripple_power.drop(idxs, inplace=True)

        idxs = ripple_power.loc[
               ripple_power.index[(ripple_power.index > behavior.iloc[-1].name) | (ripple_power.index < behavior.iloc[0].name)],
               :].index

        if idxs.shape[0] > 0:
            pos_area.drop(idxs, inplace=True)
            neg_area.drop(idxs, inplace=True)
            ripple_power.drop(idxs, inplace=True)

        ripples = ripples_calcs[session_id][3]

        #ripples = ripples[ripples["Peak frequency (Hz)"]> 120]

        # create M-L axis, subtract midpoint
        ripples["M-L (µm)"] = ripples["L-R (µm)"].copy() - 5691.510009765625

        ripples_calcs[session_id][3] = clean_ripples_from_running(ripples, start_running, stop_running, behavior)
        print(f"{session_id}: ripples number on best probe:{ripples[ripples['Probe number']==sel_probe].shape[0]}")

    return ripples_calcs


def clean_ripples_from_running(ripples, start_running, stop_running, behavior):
    for q in range(len(start_running)):
        idxs = ripples.loc[
               ripples.index[(ripples["Start (s)"] > start_running[q]) & (ripples["Start (s)"] < stop_running[q])],
               :].index
        if idxs.shape[0] > 0:
            ripples.drop(idxs, inplace=True)

    #  clean ripples if behavior data not available
    idxs = ripples.loc[
           ripples.index[(ripples["Start (s)"] > behavior.iloc[-1].name) | (ripples["Start (s)"] < behavior.iloc[0].name)],
           :].index
    if idxs.shape[0] > 0:
        ripples.drop(idxs, inplace=True)

    return ripples.reset_index(drop=True)

def process_spiking_rate_per_ripple(spikes_per_ripple, Sessions, area, sessions):

    out= []
    for genotype in sessions["full_genotype"].unique():
        print(genotype)
        for session_id in tqdm(sessions[sessions["full_genotype"]==genotype].index):
            try:
                units = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].units
                spikes_per_ripple_session = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].time_space_sub_spike_times
                spikes_per_ripple_session_table = pd.DataFrame({k:[len(i[k][i[k]<0])/(window_spike_per_ripple[0]*10)for i in spikes_per_ripple_session]
                                                                for k in dict(ChainMap(*spikes_per_ripple_session)).keys()})
                ##only active in at least 10% of ripples
                #spikes_per_ripple_session_table =  spikes_per_ripple_session_table.loc[:,spikes_per_ripple_session_table.astype(bool).sum()/len(spikes_per_ripple_session)>0.1]

                idx_inh = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]<waveform_dur_thr].index]
                idx_exc = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]>=waveform_dur_thr].index]
                out.append([genotype, session_id, spikes_per_ripple_session_table.loc[:, idx_exc].mean().mean(), spikes_per_ripple_session_table.loc[:, idx_inh].mean().mean()])
            except Exception as e:

                print(e)
                print(session_id, " not analyzed")

    t = pd.DataFrame(out, columns=["Genotype", "Session id", "Mean spiking rate exc", "Mean spiking rate inh"])
    t["Type"] = "Pre-ripple start"

    out= []
    for genotype in sessions["full_genotype"].unique():
        print(genotype)
        for session_id in tqdm(sessions[sessions["full_genotype"]==genotype].index):
            try:
                units = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].units
                spikes_per_ripple_session = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].time_space_sub_spike_times
                spikes_per_ripple_session_table = pd.DataFrame({k:[len(i[k][i[k]>0])/(window_spike_per_ripple[0]*10)for i in spikes_per_ripple_session]
                                                                for k in dict(ChainMap(*spikes_per_ripple_session)).keys()})
                ##only active in at least 10% of ripples
                #spikes_per_ripple_session_table =  spikes_per_ripple_session_table.loc[:,spikes_per_ripple_session_table.astype(bool).sum()/len(spikes_per_ripple_session)>0.1]

                idx_inh = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]<waveform_dur_thr].index]
                idx_exc = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]>=waveform_dur_thr].index]
                out.append([genotype, session_id, spikes_per_ripple_session_table.loc[:, idx_exc].mean().mean(), spikes_per_ripple_session_table.loc[:, idx_inh].mean().mean()])
            except Exception as e:

                print(e)
                print(session_id, " not analyzed")

    t2 = pd.DataFrame(out, columns=["Genotype", "Session id", "Mean spiking rate exc", "Mean spiking rate inh"])
    t2["Type"] = "Post-ripple start"

    mean_spiking_rate = pd.concat([t, t2])

    mean_spiking_rate["Area"] = area

    return mean_spiking_rate


def process_participating_neurons_per_ripple(spikes_per_ripple, Sessions, area, sessions):

    New_Sessions_key = namedtuple("Sessions", ["session_id", "area", "genotype"])
    Data = namedtuple("Data", ["Number_neurons_exc", "Number_neurons_inh"])

    out = {}
    for genotype in tqdm(sessions["full_genotype"].unique()):
        print(genotype)
        for session_id in sessions[sessions["full_genotype"]==genotype].index:
            try:
                units = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].units
                spikes_per_ripple_session = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].time_space_sub_spike_times
                spikes_per_ripple_session_table = pd.DataFrame({k:[len(i[k][i[k]>0]) for i in spikes_per_ripple_session]
                                                                for k in dict(ChainMap(*spikes_per_ripple_session)).keys()}).astype(bool)
                ##only active in at least 10% of ripples
                #spikes_per_ripple_session_table =  spikes_per_ripple_session_table.loc[:,spikes_per_ripple_session_table.astype(bool).sum()/len(spikes_per_ripple_session)>0.1]

                idx_inh = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]<waveform_dur_thr].index]
                idx_exc = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]>=waveform_dur_thr].index]
                s = New_Sessions_key(session_id=session_id, area=area, genotype=genotype)

                d = Data(Number_neurons_exc=spikes_per_ripple_session_table.loc[:, idx_exc].sum(axis=1),
                         Number_neurons_inh=spikes_per_ripple_session_table.loc[:, idx_inh].sum(axis=1))

                out[s] = d

            except Exception as e:

                print(e)
                print(session_id, " not analyzed")

    return out


def process_spiking_hists_per_ripple(spikes_per_ripple, Sessions, area, sessions):

    New_Sessions_key = namedtuple("Sessions", ["session_id", "area", "genotype"])
    Data = namedtuple("Data", ["Histogram_tot", "Histogram_top", "Histogram_bottom", "Histogram_exc", "Histogram_inh", "Bins"])

    spike_hists = {}
    for genotype in tqdm(sessions["full_genotype"].unique()):
        print(genotype)
        for session_id in sessions[sessions["full_genotype"]==genotype].index:
            try:
                units = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].units
                ripples = pd.DataFrame(spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].ripples_features)

                spikes_per_ripple_session = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].time_space_sub_spike_times
                spikes_per_ripple_session_table = pd.DataFrame({k:[len(i[k][i[k]>0])/(window_spike_per_ripple[0]*10)for i in spikes_per_ripple_session]
                                                                for k in dict(ChainMap(*spikes_per_ripple_session)).keys()})
                _ = {k: [d.get(k) for d in spikes_per_ripple_session if len(d.get(k)) > 0] for k in
                     set().union(*spikes_per_ripple_session)}
                spiking_times_ripple_start_centered = {k: np.concatenate(v) if len(v) > 0 else [] for k, v in _.items()}

                spikes_per_ripple_session_top = list(compress(spikes_per_ripple_session, ripples["∫Ripple"] > ripples["∫Ripple"].quantile(.9)))
                _ = {k: [d.get(k) for d in spikes_per_ripple_session_top if len(d.get(k)) > 0] for k in
                     set().union(*spikes_per_ripple_session_top)}
                spiking_times_ripple_start_centered_top = {k: np.concatenate(v) if len(v) > 0 else [] for k, v in _.items()}

                spikes_per_ripple_session_bottom = list(compress(spikes_per_ripple_session, ripples["∫Ripple"] < ripples["∫Ripple"].quantile(.9)))
                _ = {k: [d.get(k) for d in spikes_per_ripple_session_bottom if len(d.get(k)) > 0] for k in
                     set().union(*spikes_per_ripple_session_bottom)}
                spiking_times_ripple_start_centered_bottom = {k: np.concatenate(v) if len(v) > 0 else [] for k, v in
                                                           _.items()}

                idx_inh = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]<waveform_dur_thr].index]
                idx_exc = [x for x in spikes_per_ripple_session_table.columns if x in units[units["waveform_duration"]>=waveform_dur_thr].index]

                spiking_times_ripple_start_centered_exc = {k: spiking_times_ripple_start_centered[k] for k in idx_exc}
                spiking_times_ripple_start_centered_inh = {k: spiking_times_ripple_start_centered[k] for k in idx_inh}
                print("Check consistency: ", len(spiking_times_ripple_start_centered_exc) + len(spiking_times_ripple_start_centered_inh) == len(
                    spiking_times_ripple_start_centered))

                bins = np.linspace(-window_spike_per_ripple[0], window_spike_per_ripple[1], num=100)

                counts_tot, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered.values())),
                                                bins=bins)
                counts_tot = counts_tot / len(spiking_times_ripple_start_centered)

                counts_top, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered_top.values())),
                                                bins=bins)

                counts_top = counts_top / len(spiking_times_ripple_start_centered_top)

                counts_bottom, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered_bottom.values())),
                                                bins=bins)

                counts_bottom = counts_bottom / len(spiking_times_ripple_start_centered_bottom)

                counts_exc, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered_exc.values())),
                                                bins=bins)

                counts_exc = counts_exc / len(spiking_times_ripple_start_centered_exc)

                counts_inh, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered_inh.values())),
                                                bins=bins)
                counts_inh = counts_inh / len(spiking_times_ripple_start_centered_inh)

                s = New_Sessions_key(session_id=session_id, area=area, genotype=genotype)

                d = Data(Histogram_tot=counts_tot/len(spikes_per_ripple_session), Histogram_top=counts_top/len(spikes_per_ripple_session_top),
                         Histogram_bottom=counts_bottom/len(spikes_per_ripple_session_bottom),
                         Histogram_exc=counts_exc/len(spikes_per_ripple_session), Histogram_inh=counts_inh/len(spikes_per_ripple_session), Bins=bins)

                spike_hists[s] = d

            except Exception as e:

                print(e)
                print(session_id, " not analyzed")

    return spike_hists


def process_spiking_hists_per_ripple_per_area(spikes_per_ripple, Sessions, area, sessions):

    New_Sessions_key = namedtuple("Sessions", ["session_id", "area", "genotype"])
    Data = namedtuple("Data", ["Histogram", "Bins"])

    spike_hists = {}
    for genotype in tqdm(sessions["full_genotype"].unique()):
        print(genotype)
        for session_id in sessions[sessions["full_genotype"]==genotype].index:
            if session_id in [key.session_id for key in spikes_per_ripple.keys() if key.area==area]:

                units = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].units

                spikes_per_ripple_session = spikes_per_ripple[Sessions(session_id=session_id, area=area, source_position='medial')].time_space_sub_spike_times
                spikes_per_ripple_session_table = pd.DataFrame({k:[len(i[k][i[k]>0])/(window_spike_per_ripple[0]*10)for i in spikes_per_ripple_session]
                                                                for k in dict(ChainMap(*spikes_per_ripple_session)).keys()})
                _ = {k: [d.get(k) for d in spikes_per_ripple_session if len(d.get(k)) > 0] for k in
                     set().union(*spikes_per_ripple_session)}
                spiking_times_ripple_start_centered_parent_area = {k: np.concatenate(v) if len(v) > 0 else [] for k, v in _.items()}

                bins = np.linspace(-window_spike_per_ripple[0], window_spike_per_ripple[1], num=100)

                print(units[units["parent area"]==area]["ecephys_structure_acronym"].unique())
                for target_area in units[units["parent area"]==area]["ecephys_structure_acronym"].unique():

                        idx_area = [x for x in spikes_per_ripple_session_table.columns if x in units[units["ecephys_structure_acronym"]==target_area].index]

                        spiking_times_ripple_start_centered = {k: spiking_times_ripple_start_centered_parent_area[k] for k in idx_area}

                        print("Area processed: ", target_area, ", number units: ", len(idx_area), ", session: ", session_id)

                        counts_tot, bins = np.histogram(np.concatenate(list(spiking_times_ripple_start_centered.values())),
                                                        bins=bins)
                        counts_tot = counts_tot / len(spiking_times_ripple_start_centered)

                        s = New_Sessions_key(session_id=session_id, area=target_area, genotype=genotype)

                        d = Data(Histogram=counts_tot/len(spikes_per_ripple_session), Bins=bins)

                        spike_hists[s] = d


    return spike_hists


def process_spike_per_trial(ecephys_session_id, window_size, del_window_size, spikes_field, target_field, cache,
                            predict_field=None):

    session = cache.get_ecephys_session(ecephys_session_id=ecephys_session_id)

    metadata = session.metadata
    trials = session.trials

    stimulus_presentations = session.stimulus_presentations
    stimulus_presentations["codes_images"] = pd.Categorical(stimulus_presentations["image_name"]).codes

    # sub_stimulus_presentations = stimulus_presentations[stimulus_presentations["start_frame"].isin(trials["change_frame"])]
    # sub_trials = trials[trials["change_frame"].isin(sub_stimulus_presentations["start_frame"])]

    # change_times = sub_trials.apply(get_change_time_from_stim_table, sub_stimulus_presentations=sub_stimulus_presentations, axis=1)
    # trials['change_time'] =  change_time
    trials['change_time'] = trials.apply(get_change_time_from_stim_table, stimulus_presentations=stimulus_presentations,
                                         axis=1)

    units = session.get_units()
    units = units[units["quality"] != "noise"]

    channels = session.get_channels()
    unit_channels = units.merge(channels, left_on='peak_channel_id', right_index=True)
    unit_channels["parent_area"] = unit_channels["structure_acronym"].apply(lambda area: acronym_to_main_area(area))

    unit_channels = unit_channels[unit_channels["structure_acronym"] != "root"]

    spike_times = session.spike_times
    spike_times = {k: spike_times[k] for k in (unit_channels.index)}

    type_func = lambda row: 'miss' if row['miss'] else (
        'hit' if row['hit'] else ('false_alarm' if row['false_alarm'] else 'aborted'))

    sub_stim_pres = stimulus_presentations[
        (stimulus_presentations["omitted"] == False) & (stimulus_presentations["stimulus_block"] == 0)]

    sub_stim_pres = sub_stim_pres.merge(trials[["change_time", "miss", "hit", "false_alarm", "aborted"]],
                                        left_on="start_time", right_on="change_time", how="left")
    sub_stim_pres.drop_duplicates(inplace=True)
    sub_stim_pres.reset_index(drop=True, inplace=True)

    sub_stim_pres['type'] = sub_stim_pres.apply(type_func, axis=1)

    pres_times = sub_stim_pres["start_time"]

    # Define a lambda function to add a random float to each element
    add_random_float = lambda x: x + random.uniform(60.0, 500.9)

    # Apply the lambda function to each element of the series
    # pres_times_scrambled = pres_times.apply(add_random_float)

    spike_counts = {}
    # spike_counts_scrambled = {}

    for neuron_index, neuron_spikes in spike_times.items():
        spike_counts[neuron_index] = np.searchsorted(neuron_spikes,
                                                     pres_times + window_size + del_window_size) - np.searchsorted(
            neuron_spikes, pres_times + del_window_size)
        # spike_counts_scrambled[neuron_index] = np.searchsorted(neuron_spikes, pres_times_scrambled + window_size) - np.searchsorted(neuron_spikes, pres_times_scrambled)

    if predict_field is not None:

        sub_stim_second_pres = stimulus_presentations[
            (stimulus_presentations["omitted"] == False) & (stimulus_presentations["stimulus_block"] == 5)]
        second_pres_times = sub_stim_second_pres["start_time"]

        spike_counts_second_pres = {}

        for neuron_index, neuron_spikes in spike_times.items():
            spike_counts_second_pres[neuron_index] = np.searchsorted(neuron_spikes,
                                                                     second_pres_times + window_size + del_window_size) - np.searchsorted(
                neuron_spikes, second_pres_times + del_window_size)

        spike_counts_df = pd.DataFrame(spike_counts)
        # spike_counts_scrambled_df = pd.DataFrame(spike_counts_scrambled)
        spike_counts_second_pres_df = pd.DataFrame(spike_counts_second_pres)
        spikes_trials = xr.Dataset(data_vars={"spikes": (["stimulus_presentations", "neurons"], spike_counts_df),
                                              "second_pres_spikes": (
                                              ["stimulus_presentations", "neurons"], spike_counts_second_pres_df)},
                                   coords={"stimulus_presentations": (
                                   "stimulus_presentations", sub_stim_pres["image_name"]),
                                           "neurons": ("neurons", list(spike_times.keys()))},
                                   attrs=dict(window_size=window_size, del_window_size=del_window_size))

        spikes_trials = spikes_trials.assign_coords(
            {"parent_area": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["parent_area"]),
             "waveform_duration": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["waveform_duration"]),
             "area": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["structure_acronym"]),
             "image": ("stimulus_presentations", sub_stim_pres["image_name"]),
             "type": ("stimulus_presentations", sub_stim_pres["type"]),
             "hit": ("stimulus_presentations", np.nan_to_num(sub_stim_pres["hit"].astype(float), False).astype(bool)),
             "miss": ("stimulus_presentations", np.nan_to_num(sub_stim_pres["miss"].astype(float), False).astype(bool)),
             "is_change": (
             "stimulus_presentations", np.nan_to_num(sub_stim_pres["is_change"].astype(float), False).astype(bool))})

        zero_var_cols = np.where(spikes_trials[spikes_field].var(dim='stimulus_presentations').to_numpy() == 0)[0]
        spikes_trials = spikes_trials.sel(
            neurons=spikes_trials.neurons[~np.in1d(spikes_trials.neurons, spikes_trials.neurons[zero_var_cols])])

        X = spikes_trials[spikes_field].values
        y = spikes_trials[target_field].values

        mi = mutual_info_classif(X, y)
        spikes_trials = spikes_trials.assign_coords({"mutual_information_" + target_field: ("neurons", mi)})

        X = spikes_trials[predict_field].values
        y = spikes_trials[target_field].values

        mi = mutual_info_classif(X, y)
        spikes_trials = spikes_trials.assign_coords(
            {"mutual_information_" + target_field + "_" + predict_field: ("neurons", mi)})
    else:

        spike_counts_df = pd.DataFrame(spike_counts)

        spikes_trials = xr.Dataset(data_vars={"spikes": (["stimulus_presentations", "neurons"], spike_counts_df)},
                                   coords={"stimulus_presentations": (
                                   "stimulus_presentations", sub_stim_pres["image_name"]),
                                           "neurons": ("neurons", list(spike_times.keys()))},
                                   attrs=dict(window_size=window_size, del_window_size=del_window_size))

        spikes_trials = spikes_trials.assign_coords(
            {"parent_area": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["parent_area"]),
             "waveform_duration": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["waveform_duration"]),
             "area": ("neurons", unit_channels.loc[spikes_trials["neurons"]]["structure_acronym"]),
             "image": ("stimulus_presentations", sub_stim_pres["image_name"]),
             "type": ("stimulus_presentations", sub_stim_pres["type"]),
             "hit": ("stimulus_presentations", np.nan_to_num(sub_stim_pres["hit"].astype(float), False).astype(bool)),
             "miss": ("stimulus_presentations", np.nan_to_num(sub_stim_pres["miss"].astype(float), False).astype(bool)),
             "is_change": (
             "stimulus_presentations", np.nan_to_num(sub_stim_pres["is_change"].astype(float), False).astype(bool))})

        zero_var_cols = np.where(spikes_trials[spikes_field].var(dim='stimulus_presentations').to_numpy() == 0)[0]
        spikes_trials = spikes_trials.sel(
            neurons=spikes_trials.neurons[~np.in1d(spikes_trials.neurons, spikes_trials.neurons[zero_var_cols])])

        X = spikes_trials[spikes_field].values
        y = spikes_trials[target_field].values

        mi = mutual_info_classif(X, y)
        spikes_trials = spikes_trials.assign_coords({"mutual_information_" + target_field: ("neurons", mi)})

    return spikes_trials


def find_closest_smaller_value(s_arr, target):
    closest_value = -np.inf
    for i in range(len(s_arr)):
        if s_arr[i] < target:
            closest_value = max(closest_value, s_arr[i])
    return closest_value


def get_change_time_from_stim_table(row, stimulus_presentations):
    '''
    Given a particular row in the trials table,
    find the corresponding change time in the
    stimulus presentations table
    '''

    change_frame = row['change_frame']

    if np.isnan(change_frame):
        if row["lick_times"].shape[0] > 0:
            closest_value = find_closest_smaller_value(
                stimulus_presentations[stimulus_presentations["omitted"] == False]["start_time"].values,
                row["lick_times"][0])
            return closest_value
        else:
            return np.nan
    else:

        change_time = stimulus_presentations[stimulus_presentations["start_frame"] == change_frame] \
            ['start_time'].values[0]

        return change_time



def count_spikes(spike_times, pres_times, window_size):
    spike_counts = np.zeros(len(pres_times))
    for i in range(len(pres_times)):
        spike_counts[i] = np.searchsorted(spike_times, pres_times[i] + window_size) - np.searchsorted(spike_times, pres_times[i])
    return spike_counts

def decoding_image_identity(spikes_trials, spikes_field ,target_field):
    if spikes_field is not None:
        X = spikes_trials[spikes_field].values
        y = spikes_trials[target_field].values

        model = LinearDiscriminantAnalysis()
        n_folds = 5

        # Define the train-test split for each fold
        test_size = 0.2

        # Use cross-validation to get the mean score
        scores_ = cross_val_score(model, X, y, cv=n_folds, scoring='accuracy', n_jobs=n_folds)
    else:
        scores_ = []
    return scores_


def decoding_image_identity_across_presentations(spikes_trials, spikes_field, target_field, predict_field=None):
    if predict_field is not None:
        X = spikes_trials[spikes_field].values
        y = spikes_trials[target_field].values

        model = LinearDiscriminantAnalysis()

        model.fit(X, y)

        y_pred = model.predict(spikes_trials[predict_field].values)

        score_second_pres = accuracy_score(y, y_pred)
    else:
        score_second_pres = []

    return score_second_pres


def decoding_image_identity_brain_region(spikes_trials, spikes_field, field_to_compare, target_field,
                                         predict_field=None):
    model = LinearDiscriminantAnalysis()
    scores_first = {}
    scores_second = {}
    scores_second_fit_first = {}
    for area in np.unique(spikes_trials[field_to_compare]):
        data = spikes_trials.where(spikes_trials[field_to_compare] == area, drop=True)

        n_folds = 5

        # Define the train-test split for each fold
        test_size = 0.2

        X = data[spikes_field].values
        y = data[target_field]

        # Use cross-validation to get the mean score
        scores_first[area] = cross_val_score(model, X, y, cv=n_folds, scoring='accuracy', n_jobs=n_folds)

        if predict_field is not None:

            model.fit(X, y)

            y_pred = model.predict(data[predict_field].values)

            scores_second_fit_first[area] = accuracy_score(y, y_pred)

            X = data[predict_field].values
            y = data[target_field]

            # Use cross-validation to get the mean score
            scores_second[area] = cross_val_score(model, X, y, cv=n_folds, scoring='accuracy', n_jobs=n_folds)
        else:
            scores_second = []
            scores_second_fit_first = []

    return scores_first, scores_second, scores_second_fit_first

def sub_func(spikes_trials, target_field, X_top, k, model, n_folds):

    # Use cross-validation to get the mean score
    scores_ = cross_val_score(model, X_top, spikes_trials[target_field].values, cv=n_folds, scoring='accuracy', n_jobs=n_folds)
    return k, scores_


def check_decoding_after_dataset_pruning(spikes_trials, spikes_field, target_field):
    # Define the train-test split for each fold
    test_size = 0.2
    n_folds = 5
    model = LinearDiscriminantAnalysis()
    spike_trials = spikes_trials.sortby("mutual_information_" + target_field, ascending=False)

    input_multiprocessing = []
    for k in itertools.chain(range(1, 21), range(30, 110, 10), [spikes_trials.neurons.shape[0]]):
        X_top = spike_trials.isel(neurons=slice(None, k))[spikes_field].values
        input_multiprocessing.append((spikes_trials, target_field, X_top, k, model, n_folds))

    n_jobs = 10  # or set to -1 to use all available CPUs

    # Call the function on each dataset in parallel
    results = Parallel(n_jobs=n_jobs)(delayed(sub_func)(*arg) for arg in input_multiprocessing)

    scores_keep_best = dict(results)

    input_multiprocessing = []
    for k in itertools.chain([1], range(10, 110, 10), range(200, 1000, 100)):
        X_top = spikes_trials.isel(neurons=slice(None, spikes_trials.neurons.shape[0] - k))[spikes_field].values
        input_multiprocessing.append((spikes_trials, target_field, X_top, k, model, n_folds))

    # Define the number of CPUs to use for parallelization
    n_jobs = 10  # len(input_multiprocessing)  # or set to -1 to use all available CPUs

    # Call the function on each dataset in parallel
    results = Parallel(n_jobs=n_jobs)(delayed(sub_func)(*arg) for arg in input_multiprocessing)
    scores_remove_best = dict(results)

    return scores_keep_best, scores_remove_best
def process_decoding_(ecephys_session_id, window_size, del_window_size, spikes_field, target_field,Session, Data , cache, predict_field=None ):
    text = f"Processing session {ecephys_session_id}"
    """ styled_text = f"<span style='color:red'>{text}</span>"
    display(Markdown(styled_text))"""

    print(text)

    spikes_trials = process_spike_per_trial(ecephys_session_id, window_size, del_window_size, spikes_field,
                                            target_field, cache, predict_field)

    field_to_compare = "area"

    scores_first_per_area, scores_second_per_area, scores_second_fit_first_per_area = decoding_image_identity_brain_region(
        spikes_trials, spikes_field, field_to_compare, target_field, predict_field)

    field_to_compare = "parent_area"

    scores_first_per_parent_area, scores_second_per_parent_area, scores_second_fit_first_per_parent_area = decoding_image_identity_brain_region(
        spikes_trials, spikes_field, field_to_compare, target_field, predict_field)

    scores_first_pres = decoding_image_identity(spikes_trials, spikes_field, target_field)
    scores_second_pres = decoding_image_identity(spikes_trials, predict_field, target_field)

    scores_second_pres_fit_first = decoding_image_identity_across_presentations(spikes_trials, spikes_field,
                                                                                target_field, predict_field)

    scores_keep_best, scores_remove_best = check_decoding_after_dataset_pruning(spikes_trials, spikes_field,
                                                                                target_field)

    s = Session(ecephys_session_id=ecephys_session_id)

    parameters = {"window_size": window_size, "del_window_size": del_window_size, "field_to_compare": field_to_compare,
                  "target_field": target_field}

    d = Data(neurons=spikes_trials.neurons, scores_second_pres_fit_first=scores_second_pres_fit_first,
             scores_first_per_area=scores_first_per_area, scores_second_per_area=scores_second_per_area,
             scores_second_fit_first_per_area=scores_second_fit_first_per_area,
             scores_first_per_parent_area=scores_first_per_parent_area,
             scores_second_per_parent_area=scores_second_per_parent_area,
             scores_second_fit_first_per_parent_area=scores_second_fit_first_per_parent_area,
             scores_first_pres=scores_first_pres,
             scores_second_pres=scores_second_pres, scores_keep_best=scores_keep_best,
             scores_remove_best=scores_remove_best,
             parameters=parameters)

    # scores_summary_dict[s] = d

    text = f"Session  {ecephys_session_id} completed"
    """styled_text = f"<span style='color:green'>{text}</span>"
    display(Markdown(styled_text))"""
    print(text)

    return s, d





