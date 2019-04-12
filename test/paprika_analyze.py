from paprika.analysis import fe_calc
import paprika.restraints_json as rjson

import numpy as np
import json
import parmed as pmd
import pytraj as pt
import os
import sys

import logging
logging.basicConfig(level=logging.INFO)
logging.getLogger().setLevel(logging.DEBUG)
logging.info('Started logging...')

sim_type = sys.argv[1]

window_dir = "windows_apr"

topology = "cb6-but/cb6-but-dum.prmtop"
coordinates = "cb6-but/cb6-but-dum.rst7"

structure = pt.load(coordinates, topology)
rests = rjson.load_restraints()

phases = ["attach", "pull", "release"]

guest_restraints = rests[-3:]

analyze = fe_calc()
analyze.prmtop = structure.topology
analyze.trajectory = sim_type + "_prod.nc"
analyze.path = window_dir

analyze.restraint_list = guest_restraints
analyze.collect_data()
method = ["ti-block", "mbar-block"]
m = method[1]
analyze.methods = [m]
if m == "ti-block":
    analyze.quicker_ti_matrix = True
analyze.bootcycles = 2000
analyze.compute_free_energy(phases)
analyze.compute_ref_state_work(
    [guest_restraints[0], guest_restraints[1], None, None, guest_restraints[2], None]
)

results = analyze.results

fraction_fe_matrices = {
    "attach": "",
    "pull": "",
    "release": ""
}


fraction_sem_matrices = {
    "attach": "",
    "pull": "",
    "release": ""
}

if 'attach' in phases:
    fraction_fe_matrices['attach'] = results['attach'][m]['fraction_fe_matrix']
    fraction_sem_matrices['attach'] = results['attach'][m]['fraction_sem_matrix']

if 'pull' in phases:
    fraction_fe_matrices['pull'] = results['pull'][m]['fraction_fe_matrix']
    fraction_sem_matrices['pull'] = results['pull'][m]['fraction_sem_matrix']

if 'release' in phases:
    fraction_fe_matrices['release'] = results['release'][m]['fraction_fe_matrix']
    fraction_sem_matrices['release'] = results['release'][m]['fraction_sem_matrix']

ref_state_work = results['ref_state_work']

matrices = {
    'fraction_fe_matrices': fraction_fe_matrices,
    'fraction_sem_matrices': fraction_sem_matrices,
    'ref_state_work': ref_state_work
}

def npArraysToLists(d):
    for item in d:
        if isinstance(d[item],dict):
            npArraysToLists(d[item])
        if isinstance(d[item],np.ndarray):
            d[item] = d[item].tolist()

npArraysToLists(matrices)

with open("analysis_{0}.json".format(sim_type), 'w') as f:
    json.dump(matrices, f)
