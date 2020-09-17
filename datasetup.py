#!/usr/bin/env python

__all__ = [
    "datasetup", "namekey", "topologies",
    "workflows", "reference", "heavies",
    "protein", "alphies", "ns_per_step"
]

# When new (smaller) traj files for aswanalysis are made
# they will use this identifier in filename
namekey = '-stride-'

heavies = "mass > 5"
alphies = "name CA"
protein = "protein"

top_protein = "data/topo/chignolin-protein.pdb"
top_heavies = "data/topo/chignolin-heavy.pdb"

ns_per_step = 0.020

topologies = {
    heavies: top_heavies,
    protein: top_protein,
    alphies: None,
}

reference = {heavies: "data/topo/chignolin-folded.pdb"}

workflows = dict(
    ltj = {
      "timestep"  : 0.020,# of saved frames, unit: nanosecond
      "topfile"   : top_protein,
      "selection" : protein,
      "directory" : "./",
      "filename"  : "data/traj/longtraj/protein.*-*.dcd",
      "n_trajs"   : 1,# how many per workload?
    },

    xma_Ca_3 = {
    #xma_Ca = {
      "timestep"  : 0.020,# of saved frames, unit: nanosecond
      "topfile"   : "../../../%s"%top_protein,
      "filename"  : "trajs/*/protein.dcd",
      "selection" : protein,
      "database"  : 'mongo',

      #"directory" : "./data/traj/cgn-Ca/",
      #"n_trajs"   : 24,# how many per workload?
      #"projectname": 'cgn-Ca',
      #"last"      : 120,# last  traj idx for these workloads
      "directory" : "./data/traj/cgn-3/",
      "n_trajs"   : 12,
      "projectname": 'cgn-as2',
      "last"      : 60,# last  traj idx for these workloads
    },

    #umi_Ca_2 = {
    umi_Ca_3 = {
    #umi_Ca = {
      "timestep"  : 0.020,# of saved frames, unit: nanosecond
      "topfile"   : "../../../%s"%top_protein,
      "filename"  : "trajs/*/protein.dcd",
      "selection" : protein,
      "database"  : 'mongo',
        
      #"directory" : "./data/traj/cgn-Ca/",
      #"n_trajs"   : 12,# how many per workload?
      #"projectname": 'cgn-Ca',
      #"first"      : 120,# first traj idx for these workloads
      "directory" : "./data/traj/cgn-3",
      #"directory" : "./data/traj/cgn-2/",
      "n_trajs"   : 24,
      "projectname": 'cgn-as2',
      "first"     : 60,# first traj idx for these workloads
    },
)

datasetup = {
    # This is the 'reference state'
    "reference"    : reference,
    "topologies"   : topologies,
    # Next 2 are single long trajs
    "workflows"    : workflows,
}

