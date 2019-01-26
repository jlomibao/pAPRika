import parmed as pmd
from paprika import tleap
from paprika import restraints
from paprika import amber
from paprika.openmm import *
from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm.openmm as mm
from simtk.unit import *
from sys import stdout
from mdtraj.reporters import NetCDFReporter
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
    #guest_angle_rest = [D[1], D[0], G[0]]
    guest_restraint_atoms = [
        guest_distance_rest,
        #guest_angle_rest
    ]
    guest_restraint_targets = [guest_init_dist]
    guest_restraint_target_final = [18.0+guest_init_dist]
    guest_restraint_distance_fc = 5.0
    guest_restraint_angle_fc = 100.0

    guest_restraints = []
    for index, atoms in enumerate(guest_restraint_atoms):
        if len(atoms) > 2:
            angle = True
        else:
            angle = False
        this = restraints.DAT_restraint()
        this.auto_apr = True
        this.amber_index = True
        this.topology = structure
        this.mask1 = atoms[0]
        this.mask2 = atoms[1]
        if angle:
            this.mask3 = atoms[2]
            this.attach['fc_final'] = guest_restraint_angle_fc
        else:
            this.attach['fc_final'] = guest_restraint_distance_fc
        this.attach['target'] = guest_restraint_targets[index]
        this.attach['fraction_list'] = attach_fractions
  
        this.pull['target_final'] = guest_restraint_target_final[index]
        this.pull['num_windows'] = windows[1] 

        this.initialize()
        guest_restraints.append(this)

    return guest_restraints

def create_directory():
    work_dir = 'amber_test'
    if os.path.exists(work_dir):
        shutil.rmtree(work_dir)
    os.makedirs(work_dir)
    path = './'+work_dir+'/'

    guest_restraints = setup_restraints()
    window_list = restraints.create_window_list(guest_restraints)
    for window in window_list:
        os.makedirs(path+window)
        
    struct_dir = "cb6-but"
    topology_file = "cb6-but-dum.prmtop"
    coordinate_file = "cb6-but-dum.rst7"
    
    pull_distances = [0.0,5.0,18.0]
    for index, window in enumerate(window_list):
        shutil.copyfile(struct_dir+'/'+topology_file,
                        path+window+'/'+topology_file)
        if window[0] == 'a':
            shutil.copyfile(struct_dir+'/'+coordinate_file,
                            path+window+'/'+coordinate_file)
        elif window[0] == 'p':
            structure = pmd.load_file(struct_dir+'/'+topology_file,
                                      struct_dir+'/'+coordinate_file)
            target_difference = pull_distances[index]
            for atom in structure.atoms:
                if atom.residue.name == 'BUT':
                    atom.xz += target_difference
            structure.save(path+window+'/'+coordinate_file)
        print(window)
        with open(path+window+'/restraints.in', 'w') as f:
            #problem using amber_restraint_line
            for rest in guest_restraints:
                string = restraints.amber_restraint_line(rest,window)
                if string is not None:
                    f.write(string)

def amber_minimize_with_restraint():
    work_dir = 'amber_test'
    file_prefix = 'cb6-but-dum'
    path = './'+work_dir+'/'
    window_list = os.listdir(path)
    
    for window in window_list:
        sim = amber.Simulation()
        sim.path = path+'/'+window
        sim.prefix = 'restraint_min'
        sim.executable = 'sander'
        sim.inpcrd = file_prefix+'.rst7'
        sim.ref = file_prefix+'.rst7'
        sim.topology = file_prefix+'.prmtop'
        sim.restraint_file = 'restraints.in'
        sim.config_gb_min()
        sim.cntrl['ntr'] = 1
        sim.cntrl['maxcyc'] = 5000
        sim.cntrl['ncyc'] = 1000
        sim.cntrl['restraint_wt'] = 50.0
        sim.cntrl['restraintmask'] = "'@DUM'"
        sim.run(fail_ok=False)

def openmm_minimize_with_restraint():
    work_dir = 'amber_test'
    file_prefix = 'cb6-but-dum'
    path = './'+work_dir+'/'
    topology_file = file_prefix+'.prmtop'
    coordinate_file = file_prefix+'.rst7'
    window_list = os.listdir(path)
    
    for window in window_list: 
        structure = pmd.load_file(work_dir+'/'+window+'/'+topology_file,
                                  work_dir+'/'+window+'/'+coordinate_file)

        guest_restraints = setup_restraints()
        
        for i, atom in enumerate(structure.atoms):
            if atom.name == 'DUM':
                atom.mass = 0.0

        prmtop = AmberPrmtopFile(work_dir+'/'+window+'/'+topology_file)
        inpcrd = AmberInpcrdFile(work_dir+'/'+window+'/'+coordinate_file)
        system = prmtop.createSystem(
            nonbondedMethod=NoCutoff,
            constraints=HBonds,
            implicitSolventSaltConc=0.0*moles/liter
        )
        integrator = LangevinIntegrator(
            300*kelvin,
            1/picosecond,
            0.002*picoseconds
        )
        sim = OpenMM_GB_simulation()
        OpenMM_GB_simulation.add_openmm_restraints(sim,system,guest_restraints,'pull',int(window[1:]))
        simulation = Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('Reference'))
        simulation.context.setPositions(inpcrd.positions)
        simulation.minimizeEnergy()

#create_directory()
#amber_minimize_with_restraint()
#openmm_minimize_with_restraint()
