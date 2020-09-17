#!/usr/bin/env python
# coding: utf-8

# Builtins
import traceback
from pprint import pprint, pformat
from itertools import chain
from pathlib import Path
# Who needs these!?
import warnings
warnings.filterwarnings('ignore')

# Installed
import colorama
import pandas
import seaborn
from scipy.stats import pearsonr

# 2lazy 2change names throughout notebooks
import aswanalysis.moretools as aswa_tools
import aswanalysis.workflow as worflow

# TODO use a helper function to
# export common packages to top
# namespace and keep rest of aswa
# bound to parent module name
import aswanalysis as aswa
matplotlib = aswa.plots.matplotlib
cm         = aswa.plots.matplotlib.cm
plt        = aswa.plots.plt
mdtraj     = aswa.mdtraj
pyemma     = aswa.pyemma
coor       = aswa.coor
np         = aswa.np
msm        = aswa.msm
subprocess = aswa.subprocess
shlex      = aswa.shlex

# THE WORKFLOW DATA
from datasetup import *


#==============================================================================#
#==============================================================================#
#==============================================================================#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#==============================================================================#
#------------------------------------------------------------------------------#
print(colorama.Back.CYAN+\
    "Notebook-level printing will follow a scheme", flush=True)
print(
    "with " + \
    colorama.Back.GREEN   + "messages"  + colorama.Back.RESET + ", " + \
    colorama.Back.MAGENTA + "data info" + colorama.Back.RESET + ", " + \
    colorama.Back.RED     + "errors"    + colorama.Back.RESET + ", " + \
    colorama.Back.YELLOW  + "warnings"  + colorama.Back.RESET
, flush=True)
print(colorama.Fore.BLUE + \
    "MAYBE other colorings will be used by notebook also", flush=True)
print("but plain text messages should arise from packages the", flush=True)
print("notebook is running. Some notebook messages will guide you", flush=True)
print("to understanding if problems are happening that violate", flush=True)
print("expectations of the datasets but aren't actually errors." + \
    colorama.Fore.RESET, flush=True)


#------------------------------------------------------------------------------#
analysis_directory = Path("analyses")

_model_name = lambda feat, kwargs: "__".join([feat] + [
    "%s-%s"%(str(k),str(v)) for k,v in kwargs.items()])

print(colorama.Back.GREEN+"Processing these workflow datasets", flush=True)
print(colorama.Back.GREEN+"that were defined in file: datasetup.py", flush=True)
print(colorama.Back.MAGENTA+\
    ", ".join(["%s"%wnm for wnm in workflows]), flush=True)

aswa.create_trajlist(workflows, filterkey=lambda fnm: "stride" not in fnm)
aswa.assign_stride(workflows)
aswa.determine_epochsize(workflows)
aswa.determine_datashape(workflows, atomselection=heavies)

#------------------------------------------------------------------------------#
feat_CB  = "invca_ba"
feat_Ca  = "invca"
features = feat_Ca, feat_CB

feature_selections = {
    feat_Ca : {"add_inverse_distances": {"select_Ca": None}},
    feat_CB : [
        {"add_inverse_distances": {"select_Ca": None}},
        {"add_sidechain_torsions": None,
         "kwargs": {"which": ["chi1"]}},
        {"add_backbone_torsions": None},
]}


#------------------------------------------------------------------------------#
all_models = dict({kk:
    dict({
        k : aswa_tools.copy_include(v)
        for k,v in workflows.items()
    }) for kk in feature_selections
})

print(colorama.Back.MAGENTA + \
    "Here is your data:" + colorama.Back.RESET, flush=True)
print(colorama.Fore.BLUE + \
    "Different Feature Selections to apply to trajectories", flush=True)
print(colorama.Fore.MAGENTA + \
    pformat(feature_selections), flush=True)
print(colorama.Fore.BLUE + \
    "Pre-processed dataset descriptions" + colorama.Fore.RESET, flush=True)

#------------------------------------------------------------------------------#
for feat,wnm,w in aswa_tools.iter_models(all_models):
    da = analysis_directory / wnm
    if not da.is_dir():
        print(colorama.Back.GREEN+"Creating analysis directory", flush=True)
        da.mkdir()

for nm, dataset in all_models[list(feature_selections)[-1]].items():
    print(colorama.Back.BLUE + nm + colorama.Back.RESET, flush=True)
    print(colorama.Fore.MAGENTA + \
        pformat(dataset) + colorama.Fore.RESET, flush=True)

for feat, featz in feature_selections.items():
    for nm, dataset in all_models[feat].items():
        dataset["data_reader"],dataset["data_order"] = aswa.prepare_tica_inputs(
            {nm:dataset}, topologies[heavies], features=featz)


#------------------------------------------------------------------------------#
print(colorama.Back.CYAN + \
    "This parameter block was set by running through", flush=True)
print("the notebook calculations, going back and adjusting", flush=True)
print("and iterating further to test around good values", flush=True)
print(flush=True)
print("the workflow names are: {}".format(list(workflows)), flush=True)

tica_lags  = [1, 3, 10, 25, 50, 100, 250, 500, 1000, 2500] # 0
n_clusters = [25, 50, 100, 200, 300]      # for k-means cluste
d_clusters = [1.3, 2.0, 2.75, 3.5, 5]     # for regspace clust
n_macrostates = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
n_tica_dim = 5

print(colorama.Back.RESET + colorama.Fore.LIGHTBLUE_EX + \
    "tica_lags    : " + colorama.Back.LIGHTBLUE_EX + \
    colorama.Fore.WHITE + pformat(tica_lags), flush=True)
print(colorama.Back.RESET + colorama.Fore.LIGHTBLUE_EX + \
    "n_clusters   : " + colorama.Back.LIGHTBLUE_EX + \
    colorama.Fore.WHITE + pformat(n_clusters), flush=True)
print(colorama.Back.RESET + colorama.Fore.LIGHTBLUE_EX + \
    "d_clusters   : " + colorama.Back.LIGHTBLUE_EX + \
    colorama.Fore.WHITE + pformat(d_clusters), flush=True)
print(colorama.Back.RESET + colorama.Fore.LIGHTBLUE_EX + \
    "n_macrostates: " + colorama.Back.LIGHTBLUE_EX + \
    colorama.Fore.WHITE + pformat(n_macrostates), flush=True)

included   = lambda nm: nm in workflows # ALL

#------------------------------------------------------------------------------#
chosen_lag = 1000 # from tica, collective coords for clustering
keep_dims  = 3   # from tica cumvar, n dims of collective coords
chosen_tau = 1000 # from its, lag for MSM
chosen_k   = 50 # from k-means, n states for MSM
chosen_d   = 2.75 # from reg-space, n states for MSM
chosen_clustering = "cluster_regspace"
    
#------------------------------------------------------------------------------#
calc_name = "tica"
parameter = dict(lag=tica_lags)
wkf_pars  = lambda nm: {"reversible": False} if nm else dict()
fixed     = dict(dim=15) # more than will be used to get good CumVar plot
input_source = "data_reader"
input_key = None
input_val = None
aswa_tools.scan_single_parameter(
    calc_name, all_models, parameter, fixed, wkf_pars,
    included, input_source, input_key, input_val)

#------------------------------------------------------------------------------#
calc_name = "vamp"
parameter = dict(lag=tica_lags)
wkf_pars  = lambda x: dict()
fixed     = dict(dim=15) # more than will be used to get good CumVar plot
input_source = "data_reader"
input_key = None
input_val = None
aswa_tools.scan_single_parameter(
    calc_name, all_models, parameter, fixed, wkf_pars,
    included, input_source, input_key, input_val)

#------------------------------------------------------------------------------#
calc_name = "cluster_kmeans"
parameter = dict(k=n_clusters)
fixed     = dict(max_iter=100)
wkf_pars  = lambda nm: {"stride": 5} if nm else dict()
input_source = "tica"
input_key = "lag"
input_val = chosen_lag
aswa_tools.scan_single_parameter(
    calc_name, all_models, parameter, fixed, wkf_pars,
    included, input_source, input_key, input_val)

#------------------------------------------------------------------------------#
calc_name = "cluster_regspace"
# This is VERY dependent on data,
# calculated right before clustering
parameter = dict(dmin=d_clusters)
fixed     = dict(max_iter=100, max_centers=1000, min_counts=[3, 10])
wkf_pars  = lambda nm: {"stride": 5} if nm else dict()
input_source = "tica"
input_key = "lag"
input_val = chosen_lag
aswa_tools.scan_single_parameter(
    calc_name, all_models, parameter, fixed, wkf_pars,
    included, input_source, input_key, input_val)

#------------------------------------------------------------------------------#
rmsd_analysis_directory = analysis_directory / "rmsd-simple"
if not analysis_directory.is_dir():
    analysis_directory.mkdir()

if not rmsd_analysis_directory.is_dir():
    rmsd_analysis_directory.mkdir()

label_heavies   = "Heavy RMSD [Å]"
label_alphies   = "Cα RMSD [Å]"
label_timesteps = "Time [μs]"
label_longtraj  = "Single Trajectory"

#------------------------------------------------------------------------------#
# Human-read from OpenMM `StateDataReporter` output
md_state_properties  = [
    "Step", "Potential Energy (kJ/mole)",
    "Kinetic Energy (kJ/mole)", "Total Energy (kJ/mole)",
    "Temperature (K)", "Box Volume (nm^3)",
    "Density (g/mL)", "Speed (ns/day)"
]
timestep, e_potential, e_kinetic, e_total, temperature, \
volume, density, md_speed = md_state_properties

md_state_types = _type  = {
    k:v for k,v in zip(md_state_properties,
        [int, float, float, float,
        float, float, float, float]
)}

#------------------------------------------------------------------------------#
step_per_ns = {label_longtraj:int(1/ns_per_step)}
us_per_step = ns_per_step / 1000
top_file    = reference[heavies]
n_features  = len(features)

_steps_to_ns = lambda s: s * ns_per_step
_index_to_time = lambda idx: idx/step_per_ns[label_longtraj]/1000

match_keys = lambda d,s: [
    key for key in d if key.find(s) > -1]


#------------------------------------------------------------------------------#
for feat,nm,dataset in aswa_tools.iter_models(all_models):
    dataset["rmsd_trajs"] = rmsd_trajs = list()
    for ftj in map(lambda f: Path(f), dataset["tica_inputs"]):
        if not ftj.is_file():
            print(colorama.Back.RED + \
                "Traj data missing: %s"%dwtj, flush=True)
            continue

        ftj_rmsd = Path("%s-rmsd-%s.npz" % (
            str(ftj).rstrip(".dcd"), feat))

        if ftj_rmsd.is_file():
            print(colorama.Fore.MAGENTA + \
                "Found RMSDs File: %s"%ftj_rmsd, flush=True)
            
            _rtj = np.squeeze(np.load(str(ftj_rmsd))["rmsd"])
            
        else:
            print(colorama.Fore.GREEN + \
                "Calculating and Saving to File: {}".format(
                ftj_rmsd), flush=True)
            
            _rtj = aswa.calc_rmsd_trajs(
                str(ftj),
                topfile=str(reference[heavies]),
            )[0]
            
            np.savez(str(ftj_rmsd), rmsd=_rtj)

        rmsd_trajs.append(10*_rtj)


#------------------------------------------------------------------------------#
traj_heavies   = mdtraj.load(
    all_models['invca']['ltj']['tica_inputs'],
    top=top_file)

c_alphas       = traj_heavies.topology.select(alphies)
folded_heavies = mdtraj.load(top_file)
folded_alphies = folded_heavies.atom_slice(c_alphas)

traj_alphies = traj_heavies.atom_slice(c_alphas)
traj_rmsd_df = pandas.DataFrame({
    label_heavies   : mdtraj.rmsd(
        traj_heavies, folded_heavies) * 10, #nm to angstrom
    label_alphies   : mdtraj.rmsd(
        traj_alphies, folded_alphies) * 10, #nm to angstrom
    label_timesteps : [
        i * us_per_step for i in range(len(traj_alphies))]
})


#------------------------------------------------------------------------------#
print(colorama.Back.LIGHTBLACK_EX + \
    colorama.Fore.BLUE + "  DONE LOADING DATA  ", flush=True)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
