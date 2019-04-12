import os
import re
import shutil

import pytraj as pt
import numpy as np
import parmed as pmd

from paprika import restraints
from paprika.restraints import DAT_restraint
from paprika.restraints import static_DAT_restraint
from paprika.restraints import create_window_list
from paprika.restraints import amber_restraint_line
import paprika.restraints_json as rjson

host = ":CB6"
guest = ":BUT"

H = [host+"@C7", host+"@C31", host+"@C19"]
G = [guest+"@C", guest+"@C3"]
D = [":DM1", ":DM2", ":DM3"]

struct_dir = "cb6-but"
topology_file = "cb6-but-dum.prmtop"
coordinate_file = "cb6-but-dum.rst7"

structure = pmd.load_file(struct_dir+'/'+topology_file,
                          struct_dir+'/'+coordinate_file,
                          structure=True)

traj = pt.load(struct_dir+'/'+coordinate_file,
               struct_dir+'/'+topology_file)
guest_init_dist = pt.distance(traj,D[0]+' '+G[0])[0]

attach_string = "0.00 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"

attach_fractions = [float(i) / 100 for i in attach_string.split()]

pull_distances= np.arange(0.0 + guest_init_dist, 18.0 + guest_init_dist, 1.0)

release_fractions = attach_fractions[::-1]

windows = [len(attach_fractions), len(pull_distances), len(release_fractions)]

# Static Restraints
static_distance_rest = [D[0], H[0]]
static_angle_rest1 = [D[1], D[0], H[0]]
static_dihedral_rest1 = [D[2], D[1], D[0], H[0]]
static_angle_rest2 = [D[0], H[0], H[1]]
static_dihedral_rest2 = [D[1], D[0], H[0], H[1]]
static_dihedral_rest3 = [D[0], H[0], H[1], H[2]]
static_restraint_atoms = [
    static_distance_rest,
    static_angle_rest1,
    static_dihedral_rest1,
    static_angle_rest2,
    static_dihedral_rest2,
    static_dihedral_rest3
]
static_restraint_distance_fc = 5.0
static_restraint_angle_fc = 100.0

# Guest Restraints
guest_distance_rest = [D[0], G[0]]
guest_angle_rest1 = [D[1], D[0], G[0]]
guest_angle_rest2 = [D[0], G[0], G[1]]
guest_restraint_atoms = [
    guest_distance_rest,
    guest_angle_rest1,
    guest_angle_rest2,
]
guest_restraint_targets = [guest_init_dist, 180.0, 180.0]
guest_restraint_target_final = [18.0+guest_init_dist, 180.0, 180.0]
guest_restraint_distance_fc = 5.0
guest_restraint_angle_fc = 100.0

# Create restraint objects
static_restraints = []
for index, atoms in enumerate(static_restraint_atoms):
    this = static_DAT_restraint(restraint_mask_list=atoms,
                               num_window_list=windows,
                               ref_structure=structure,
                               force_constant=static_restraint_angle_fc if len(atoms)>2 else static_restraint_distance_fc,
                               amber_index=True)
    static_restraints.append(this)

guest_restraints = []
for index, atoms in enumerate(guest_restraint_atoms):
    if len(atoms) > 2:
        angle = True
    else:
        angle = False
    this = DAT_restraint()
    this.auto_apr = True
    this.amber_index = True
    this.topology = structure
    this.mask1 = atoms[0]
    this.mask2 = atoms[1]
    if angle:
        this.mask3 = atoms[2]
        this.attach['fc_final'] = guest_restraint_angle_fc
        this.release['fc_final'] = guest_restraint_angle_fc
    else:
        this.attach['fc_final'] = guest_restraint_distance_fc
        this.release['fc_final'] = guest_restraint_distance_fc
    this.attach['target'] = guest_restraint_targets[index]
    this.attach['fraction_list'] = attach_fractions

    this.pull['target_final'] = guest_restraint_target_final[index]
    this.pull['num_windows'] = windows[1]

    this.release['target'] = guest_restraint_targets[index]
    this.release['fraction_list'] = release_fractions

    this.initialize()
    guest_restraints.append(this)

rjson.save_restraints(restraint_list=static_restraints+guest_restraints)

# Create directory structure
work_dir = "windows_apr" 
if os.path.exists(work_dir):
    shutil.rmtree(work_dir)
os.makedirs(work_dir)
path = './'+work_dir+'/'

window_list = create_window_list(guest_restraints)
for window in window_list:
    os.makedirs(path+window)

for window in window_list:
    shutil.copyfile(struct_dir+'/'+topology_file, path+window+'/'+topology_file)
    if window[0] == 'a':
        shutil.copyfile(struct_dir+'/'+coordinate_file,
                        path+window+'/'+coordinate_file)
        phase = "attach"
    elif window[0] == 'p':
        structure = pmd.load_file(struct_dir+'/'+topology_file,
                                  struct_dir+'/'+coordinate_file)
        target_difference = guest_restraints[0].phase["pull"]["targets"][int(window[1:])] - guest_restraints[0].pull["target_initial"]

        for atom in structure.atoms:
            if atom.residue.name == "BUT":
                atom.xz += target_difference
        structure.save(path+window+'/'+coordinate_file)
        phase = "pull"
    elif window[0] == 'r':
        structure.save(path+window+'/'+coordinate_file)
        phase = "release"
    with open(path+window+"/rest.in", 'w') as f:
        for rest in static_restraints + guest_restraints:
            string = amber_restraint_line(rest,window)
            f.write(string)
