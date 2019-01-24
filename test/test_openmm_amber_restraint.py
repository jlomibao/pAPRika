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

    attach_fractions = [0.0,1.0]
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
    rest.attach['fraction_list'] = attach_fractions
    rest.pull['target_initial'] = guest_restraint_target
    rest.pull['target_final'] = guest_restraint_target_final
    rest.pull['num_windows'] = windows[1]
    rest.initialize()
    return rest

def create_directory():
    work_dir = 'amber_test'
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    path = './'+work_dir+'/'

    rest = setup_restraints()
    window_list = restraints.create_window_list([rest])
    for window in window_list:
        os.makedirs(path+window)
        
    struct_dir = "cb6-but"
    topology_file = "cb6-but-dum.prmtop"
    coordinate_file = "cb6-but-dum.rst7"

    for window in window_list:
        shutil.copyfile(struct_dir+'/'+topology_file,
                        path+window+'/'+topology_file)
        if window[0] == 'a':
            shutil.copyfile(struct_dir+'/'+coordinate_file,
                            path+window+'/'+coordinate_file)
        elif window[0] == 'p':
            structure = pmd.load_file(struct_dir+'/'+topology_file,
                                      struct_dir+'/'+coordinate_file)
            target_difference = rest.phase['pull']['targets'][int(window[1:])]-rest.pull['target_initial']
            for atom in structure.atoms:
                if atom.residue.name == 'BUT':
                    atom.xz += target_difference
            structure.save(path+window+'/'+coordinate_file)
        print(window)
        with open(path+window+'/restraints.in', 'w') as f:
            #problem using amber_restraint_line
            string = restraints.amber_restraint_line(rest,window)
            if string is not None:
                f.write(string)

create_directory()
