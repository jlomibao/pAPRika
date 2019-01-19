import parmed as pmd
from paprika import tleap
from paprika import restraints
from paprika import amber
import pytraj as pt
import numpy as np
import os
import shutil
import re

def setup_restraints():
    struct_dir = "cb6-but"
    topology_file = "cb6-but-dum.prmtop"
    coordinate_file = "cb6-but-dum.rst7"
    
    structure = pmd.load_file(struct_dir+'/'+topology_file,
                              struct_dir+'/'+coordinate_file,
                              structure=True)
    host = ":CB6"
    guest = ":BUT"

    H = [host+"@C7", host+"@C31", host+"@C19"]
    G = [guest+"@C", guest+"@C3"]
    D = [":DM1", ":DM2", ":DM3"]

    # Get distance of dummy atom to guest with pytraj
    traj = pt.load(struct_dir+'/'+coordinate_file,struct_dir+'/'+topology_file)
    guest_init_dist = pt.distance(traj,D[0]+' '+G[0])[0]

    attach_fractions = [1.0]
    pull_distances = [guest_init_dist,5.0+guest_init_dist,18.0+guest_init_dist]
    windows = [len(attach_fractions), len(pull_distances)]

    # Guest Restraints
    guest_distance_rest = [D[0], G[0]]
    guest_restraint_target = guest_init_dist
    guest_restraint_target_final = 18.0+guest_init_dist
    guest_restraint_distance_fc = 5.0

    rest = restraints.DAT_restraint()
    rest.continuous_apr = True
    rest.amber_index = True
    rest.topology = structure
    rest.mask1 = guest_distance_rest[0]
    rest.mask2 = guest_distance_rest[1]
    rest.attach['fc_final'] = guest_restraint_distance_fc
    rest.attach['target'] = guest_restraint_target
    rest.pull['target_final'] = guest_restraint_target_final
    rest.pull['num_windows'] = windows[1]
    rest.initialize()
    return rest
