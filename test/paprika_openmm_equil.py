import os
import parmed as pmd
import pytraj as pt

from paprika import utils
from paprika import restraints
import paprika.restraints_json as rjson

from simtk.openmm import *
import simtk.openmm.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app
from parmed.openmm.reporters import NetCDFReporter
from parmed.openmm.reporters import RestartReporter

path = "./windows_apr/"
topology = "/cb6-but-dum.prmtop"
coordinates = "/cb6-but-dum.rst7"

settings = {
    "nonbonded_method": app.NoCutoff,
    "constraints": app.HBonds,
    "temperature": 298.15*unit.kelvin,
    "friction": 1/unit.picosecond,
    "timestep": 0.002*unit.picosecond,
    "implicit_solvent": app.HCT,
    "num_steps": 500000,
    "write_freq": 250,
    "traj_file": "openmm_equil.nc"
}

window_list = os.listdir(path)

rests = rjson.load_restraints()

for window in window_list:
    if window[0] == 'a':
        phase = "attach"
    elif window[0] == 'p':
        phase = "pull"
    elif window[0] == 'r':
        phase = "release"
    window_num = int(window[1:])
    
    if os.path.exists(path+window+"/openmm_equil.rst7"):
        continue
    
    integrator = LangevinIntegrator(
        settings["temperature"],
        settings["friction"],
        settings["timestep"],
    )

    barostat = mm.MonteCarloBarostat(1.0*unit.bar, settings["temperature"], 25)

    structure = pmd.load_file(path+window+topology,
                              path+window+coordinates,
                              structure=True)
    # Set mass of Dummy atoms to 0 so they are non-interacting
    for atom in structure.atoms:
        if atom.name == "DUM":
            atom.mass = 0

    topology_0m = "/cb6-but-dum-0m.prmtop"
    coordinates_0m = "/cb6-but-dum-0m.rst7"

    structure.save(path+window+topology_0m, overwrite=True)
    structure.save(path+window+coordinates_0m, overwrite=True)

    prmtop = app.AmberPrmtopFile(path+window+topology_0m)
    inpcrd = app.AmberInpcrdFile(path+window+coordinates_0m)

    system = prmtop.createSystem(
        nonbondedMethod = settings["nonbonded_method"],
        constraints = settings["constraints"],
        implicitSolvent = settings["implicit_solvent"],
    )
    static_restraints = []

    dist_restraint = mm.CustomBondForce("k*(r-r0)^2")
    dist_restraint.addPerBondParameter("k")
    dist_restraint.addPerBondParameter("r0")
    
    rest = rests[0] # first static restraint
    static_init_dist = rest.phase[phase]["targets"][window_num]
    bond_fc = rest.phase[phase]["force_constants"][window_num]
    r0 = static_init_dist * unit.angstroms
    k = bond_fc * unit.kilocalories_per_mole / unit.angstroms**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    dist_restraint.addBond(a1[0], a2[0], [k, r0])
    static_restraints.append(dist_restraint)

    angle_restraint_1 = mm.CustomAngleForce("0.5*k*(theta-theta0)^2")
    angle_restraint_1.addPerAngleParameter("k")
    angle_restraint_1.addPerAngleParameter("theta0")

    rest = rests[1] # second static restraint
    static_init_angle_1 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = static_init_angle_1 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    angle_restraint_1.addAngle(a1[0], a2[0], a3[0], [k, theta0])
    static_restraints.append(angle_restraint_1)

    dihedral_restraint_1 = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
    dihedral_restraint_1.addPerTorsionParameter("k")
    dihedral_restraint_1.addPerTorsionParameter("theta0")

    rest = rests[2] # third static restraint
    static_init_dihedral_1 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = static_init_dihedral_1 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    a4 = utils.index_from_mask(structure, mask=rest.mask4, amber_index=False)
    dihedral_restraint_1.addTorsion(a1[0], a2[0], a3[0], a4[0], [k, theta0])
    static_restraints.append(dihedral_restraint_1)

    angle_restraint_2 = mm.CustomAngleForce("0.5*k*(theta-theta0)^2")
    angle_restraint_2.addPerAngleParameter("k")
    angle_restraint_2.addPerAngleParameter("theta0")

    rest = rests[3] # fourth static restraint
    static_init_angle_2 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = static_init_angle_2 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    angle_restraint_2.addAngle(a1[0], a2[0], a3[0], [k, theta0])
    static_restraints.append(angle_restraint_2)

    dihedral_restraint_2 = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
    dihedral_restraint_2.addPerTorsionParameter("k")
    dihedral_restraint_2.addPerTorsionParameter("theta0")

    rest = rests[4] # fifth static restraint
    static_init_dihedral_2 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = static_init_dihedral_2 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    a4 = utils.index_from_mask(structure, mask=rest.mask4, amber_index=False)
    dihedral_restraint_2.addTorsion(a1[0], a2[0], a3[0], a4[0], [k, theta0])
    static_restraints.append(dihedral_restraint_2)

    dihedral_restraint_3 = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
    dihedral_restraint_3.addPerTorsionParameter("k")
    dihedral_restraint_3.addPerTorsionParameter("theta0")

    rest = rests[5] # sixth static restraint
    static_init_dihedral_3 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = static_init_dihedral_3 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    a4 = utils.index_from_mask(structure, mask=rest.mask4, amber_index=False)
    dihedral_restraint_3.addTorsion(a1[0], a2[0], a3[0], a4[0], [k, theta0])
    static_restraints.append(dihedral_restraint_3)

    guest_restraints = []

    dist_restraint_g1 = mm.CustomBondForce("k*(r-r0)^2")
    dist_restraint_g1.addPerBondParameter("k")
    dist_restraint_g1.addPerBondParameter("r0")
    
    rest = rests[6] # first guest restraint
    guest_dist = rest.phase[phase]["targets"][window_num]
    bond_fc = rest.phase[phase]["force_constants"][window_num]
    r0 = guest_dist * unit.angstroms
    k = bond_fc * unit.kilocalories_per_mole / unit.angstroms**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    dist_restraint_g1.addBond(a1[0], a2[0], [k, r0])
    guest_restraints.append(dist_restraint_g1)

    angle_restraint_g1 = mm.CustomAngleForce("0.5*k*(theta-theta0)^2")
    angle_restraint_g1.addPerAngleParameter("k")
    angle_restraint_g1.addPerAngleParameter("theta0")

    rest = rests[7] # second guest restraint
    guest_angle1 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = guest_angle1 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    angle_restraint_g1.addAngle(a1[0], a2[0], a3[0], [k, theta0])
    static_restraints.append(angle_restraint_g1)
 
    angle_restraint_g2 = mm.CustomAngleForce("0.5*k*(theta-theta0)^2")
    angle_restraint_g2.addPerAngleParameter("k")
    angle_restraint_g2.addPerAngleParameter("theta0")

    rest = rests[8] # third guest restraint
    guest_angle2 = rest.phase[phase]["targets"][window_num]
    angle_fc = rest.phase[phase]["force_constants"][window_num]
    theta0 = guest_angle2 * unit.degrees
    k = angle_fc * unit.kilocalories_per_mole / unit.radians**2

    a1 = utils.index_from_mask(structure, mask=rest.mask1, amber_index=False)
    a2 = utils.index_from_mask(structure, mask=rest.mask2, amber_index=False)
    a3 = utils.index_from_mask(structure, mask=rest.mask3, amber_index=False)
    angle_restraint_g2.addAngle(a1[0], a2[0], a3[0], [k, theta0])
    static_restraints.append(angle_restraint_g2)

    for rest in static_restraints+guest_restraints:
        system.addForce(rest)
    #system.addForce(barostat)

    #platforms = [ mm.Platform.getPlatform(index).getName() for index in range(mm.Platform.getNumPlatforms()) ]
    simulation = app.Simulation(prmtop.topology, system, integrator,
                                mm.Platform.getPlatformByName('OpenCL'))

    simulation.context.setPositions(inpcrd.positions)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(settings["temperature"])

    simulation.reporters.append(
        NetCDFReporter(path+window+'/'+settings["traj_file"],
                       settings["write_freq"])
    )

    simulation.reporters.append(
        app.StateDataReporter(
            path + window + "/openmm_equil.log",
            settings["write_freq"],
            step=True,
            potentialEnergy=True,
            temperature=True,
            speed=True
        )
    )

    simulation.reporters.append(
        RestartReporter(
            path + window + "/openmm_equil.rst7",
            settings["num_steps"]
        )
    )

    simulation.step(settings["num_steps"])
