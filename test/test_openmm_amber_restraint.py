import parmed as pmd
from paprika import tleap
from paprika import restraints
from paprika import amber
from paprika.openmm import *
from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm.openmm as mm
from simtk.unit import *
from mdtraj.reporters import NetCDFReporter
from sys import stdout
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
        open(path+window+'/no_restraints.in', 'a').close()
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

        with open(path+window+'/restraints.in', 'w') as f:
            #problem using amber_restraint_line
            for rest in guest_restraints:
                string = restraints.amber_restraint_line(rest,window)
                if string is not None:
                    f.write(string)

def amber_minimize(hasRestraint=False):
    work_dir = 'amber_test'
    file_prefix = 'cb6-but-dum'
    path = './'+work_dir+'/'
    window_list = os.listdir(path)
    
    for window in window_list:
        sim = amber.Simulation()
        sim.path = path+window
        if hasRestraint:
            sim.prefix = 'amber_restraint_min'
            sim.restraint_file = 'restraints.in'
        else:
            sim.prefix = 'amber_min'
            sim.restraint_file = 'no_restraints.in'
        sim.executable = 'sander'
        sim.inpcrd = file_prefix+'.rst7'
        sim.ref = file_prefix+'.rst7'
        sim.topology = file_prefix+'.prmtop'
        sim.config_gb_min()
        sim.cntrl['ntr'] = 1
        sim.cntrl['maxcyc'] = 5000
        sim.cntrl['ncyc'] = 5000
        sim.cntrl['restraint_wt'] = 50.0
        sim.cntrl['restraintmask'] = "'@DUM'"
        sim.run(fail_ok=False)

def openmm_minimize(hasRestraints=False):
    work_dir = 'amber_test'
    file_prefix = 'cb6-but-dum'
    path = './'+work_dir+'/'
    topology_file = file_prefix+'.prmtop'
    coordinate_file = file_prefix+'.rst7'
    window_list = os.listdir(path)

    energy_results = {}    
    for window in window_list: 
        structure = pmd.load_file(work_dir+'/'+window+'/'+topology_file,
                                  work_dir+'/'+window+'/'+coordinate_file)

        # Create Positional Restraints
        pos_rest = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        pos_rest.addGlobalParameter("k",50.0*kilocalories_per_mole/angstroms**2)
        pos_rest.addPerParticleParameter("x0")
        pos_rest.addPerParticleParameter("y0")
        pos_rest.addPerParticleParameter("z0")
        for i, atom in enumerate(structure.positions):
            if structure.atoms[i].name == 'DUM':
                pos_rest.addParticle(i, atom.value_in_unit(nanometers))

        for i, atom in enumerate(structure.atoms):
            if atom.name == 'DUM':
                atom.mass = 0.0

        structure.save(work_dir+'/'+window+'/'+file_prefix+'-0m.prmtop',overwrite=True)
        sim = OpenMM_GB_simulation()
        settings = {
            "platform": "CPU",
            "coordinates": work_dir+'/'+window+'/'+coordinate_file,
            "temperature": 300*kelvin,
            "friction": 1/picosecond,
            "timestep": 0.002*picoseconds,
            "nonbonded_method": NoCutoff,
            "solvent": None,
            "salt": 0.0*moles/liter,
            "constraints": HBonds,
        }
        sim.path = path+window
        sim.coordinates = work_dir+'/'+window+'/'+coordinate_file
        sim.topology = work_dir+'/'+window+'/'+topology_file
        sim.phase = 'pull'
        sim.window = int(window[1:])
        system = sim.setup_system(settings, seed=None)

        if hasRestraints:
            guest_restraints = setup_restraints()
            for rest in guest_restraints:
                print(rest.phase['pull']['targets'][int(window[1:])] * angstroms)
            system = sim.add_openmm_restraints(
                system,guest_restraints,sim.phase,sim.window
            )
        system.addForce(pos_rest)
        simulation = sim.setup_simulation(system,settings)
        sim.minimize(simulation,save=False)
        state = simulation.context.getState(getEnergy=True,
                                            getPositions=True,
                                            getVelocities=True)
        restrt = pmd.openmm.RestartReporter(path+window+'/'+str(hasRestraints)+'_openmm.rst7', 500)
        restrt.report(simulation,state)
        totalEnergy = state.getPotentialEnergy() + state.getKineticEnergy()
        totalEnergy = totalEnergy.in_units_of(kilocalories_per_mole)
        #print(window,"restraints =",hasRestraints,totalEnergy)
        energy_results[window] = totalEnergy
    return energy_results

def get_amber_results(hasRestraints=False):
    energy_results = {}
    work_dir = 'amber_test'
    path = './'+work_dir+'/'
    window_list = os.listdir(path)
    for window in window_list:
        if hasRestraints:
            output_file = "/amber_restraint_min.mdinfo"
        else:
            output_file = "/amber_min.mdinfo"
        with open(path+window+output_file, 'r') as f:
            filelines = f.readlines()
            for i,line in enumerate(filelines):
                if re.search('EAMBER', line):
                    cols = line.split()
                    totalEnergy = float(cols[2])*kilocalories_per_mole
                    energy_results[window] = totalEnergy
                    #print(window,"restraints =",hasRestraints,totalEnergy)
    return energy_results

def test_amber_vs_openmm(hasRestraints=False):
    work_dir = 'amber_test'
    path = './'+work_dir+'/'
    window_list = os.listdir(path)
    amber_energy = get_amber_results(hasRestraints)
    openmm_energy = openmm_minimize(hasRestraints)
    for window in window_list:
        print(window,"restraints =",hasRestraints)
        print("AMBER:",amber_energy[window])
        print("OPENMM:",openmm_energy[window])
        print()
   
#create_directory()
#amber_minimize(hasRestraint=True)
#amber_minimize(hasRestraint=False)

test_amber_vs_openmm(hasRestraints=False)
test_amber_vs_openmm(hasRestraints=True)
