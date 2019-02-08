from paprika import restraints
from paprika.openmm import *
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout
import simtk.openmm.openmm as mm
import parmed as pmd
import pytraj as pt
import numpy as np
import os

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

def openmm_minimize(window):
    work_dir = 'amber_test/'+window
    file_prefix = 'cb6-but-dum'
    path = './'+work_dir+'/'
    topology_file = file_prefix+'.prmtop'
    coordinate_file = file_prefix+'.rst7'
    
    structure = pmd.load_file(path+topology_file, path+coordinate_file)
    for i, atom in enumerate(structure.atoms):
        if atom.name == 'DUM':
            atom.mass = 0.0
    structure.save(path+file_prefix+'-0m.prmtop', overwrite=True)
    structure.save(path+file_prefix+'-0m.rst7', overwrite=True)

    prmtop = AmberPrmtopFile(path+topology_file)
    inpcrd = AmberInpcrdFile(path+coordinate_file)
    
    system = prmtop.createSystem(
        nonbondedMethod = NoCutoff,
        constraints = HBonds,
        implicitSolventSaltConc = 0.1*moles/liter
    )
    integrator = LangevinIntegrator(
        300*kelvin,
        1/picosecond,
        0.002*picoseconds
    )
    platform = mm.Platform.getPlatformByName('CPU')

    pos_rest = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    pos_rest.addGlobalParameter("k", 50.0*kilocalories_per_mole/angstroms**2)
    pos_rest.addPerParticleParameter("x0")
    pos_rest.addPerParticleParameter("y0")
    pos_rest.addPerParticleParameter("z0")
    for i, atom in enumerate(structure.positions):
        if structure.atoms[i].name in ('DUM'):
            pos_rest.addParticle(i, atom.value_in_unit(nanometers))

    host = ':CB6'
    guest = ':BUT'
    H = [host+"@C7", host+"@C31", host+"@C19"]
    G = [guest+"@C", guest+"@C3"]
    D = [":DM1", ":DM2", ":DM3"]

    guest_restraints = setup_restraints()

    traj = pt.load(path+coordinate_file, path+topology_file)
    static_distance_rest = [D[0], H[0]]
    static_init_dist = pt.distance(traj,D[0]+' '+H[0])[0]
    
    dist_restraint = mm.CustomBondForce('k * (r - r_0)^2')
    dist_restraint.addPerBondParameter('k')
    dist_restraint.addPerBondParameter('r_0')

    r_0 = static_init_dist * angstroms
    k = 5 * kilojoules_per_mole / angstroms**2

    dist_restraint.addBond(guest_restraints[0].index1[0], guest_restraints[0].index2[0], [k, r_0])

    system.addForce(pos_rest)
    system.addForce(dist_restraint)

    simulation = Simulation(
        prmtop.topology,
        system,
        integrator,
        platform,
        None,
    )
    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy()
    state = simulation.context.getState(
        getEnergy=True,
        getPositions=True,
        getVelocities=True,
    )
    restrt = pmd.openmm.RestartReporter(path+'pos_rest_openmm.rst7', 500)
    restrt.report(simulation,state)
    totalEnergy = state.getPotentialEnergy() + state.getKineticEnergy()
    totalEnergy = totalEnergy.in_units_of(kilocalories_per_mole)
    return totalEnergy
window = 'p000'
print(openmm_minimize(window))
