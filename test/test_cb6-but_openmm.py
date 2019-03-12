from paprika import utils 
from simtk.openmm import *
import simtk.openmm.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app
from mdtraj.reporters import NetCDFReporter
import parmed as pmd
import pytraj as pt

def test_openmm_cb6but_sim(num_rests=0):
    path = './cb6-but_test/guest_inside/'
    topology = 'cb6-but-dum.prmtop'
    coordinates = 'cb6-but-dum.rst7'
    
    if num_rests > 0:
        md_out = 'cb6_but_openmm_rest_{:02d}.csv'.format(num_rests)
        traj_out = 'cb6_but_openmm_rest_{:02d}.nc'.format(num_rests)
    else:
        md_out = 'cb6_but_openmm.csv'
        traj_out = 'cb6_but_openmm.nc'

    structure = pmd.load_file(path+topology, path+coordinates, structure=True)
    traj = pt.load(path+coordinates, path+topology)

    host = ":CB6"
    guest = ":BUT"

    H = [host+"@C7", host+"@C31", host+"@C19"]
    G = [guest+"@C", guest+"@C3"]
    D = [":DM1", ":DM2", ":DM3"] 

    H_i = [0,0,0]
    G_i = [0,0]
    D_i = [0,0,0]

    # Get indices for atom masks
    for i, mask in enumerate(H):
        H_i[i] = utils.index_from_mask(structure, mask, amber_index=False)[0]
    for i, mask in enumerate(G):
        G_i[i] = utils.index_from_mask(structure, mask, amber_index=False)[0]
    for i, mask in enumerate(D):
        D_i[i] = utils.index_from_mask(structure, mask, amber_index=False)[0]

    # Set mass of Dummy atoms to 0 so they are non-interacting
    for i, atom in enumerate(structure.atoms):
        if atom.name == 'DUM':
            atom.mass = 0.0

    topology_0m = 'cb6-but-dum-0m.prmtop'
    coordinates_0m = 'cb6-but-dum-0m.rst7'

    structure.save(path+topology_0m, overwrite=True)
    structure.save(path+coordinates_0m, overwrite=True)

    prmtop = app.AmberPrmtopFile(path+topology_0m)
    inpcrd = app.AmberInpcrdFile(path+coordinates_0m)

    settings = {
        'nonbonded_method': app.NoCutoff,
        'temperature': 298.15*unit.kelvin,
        'friction': 1/unit.picosecond,
        'timestep': 0.002*unit.picosecond,
        'implicit_solvent': app.HCT,
        'numsteps': 500000
    }

    system = prmtop.createSystem(
        nonbondedMethod = settings['nonbonded_method'],
        implicitSolvent = settings['implicit_solvent'],
        removeCMMotion=False,
    )

    integrator = LangevinIntegrator(
        settings['temperature'],
        settings['friction'],
        settings['timestep']
    )

    # Create Positional Restraints for Dummy atoms
    pos_restraint = mm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    pos_restraint.addGlobalParameter('k', 50.0*unit.kilocalories_per_mole/unit.angstroms**2)

    pos_restraint.addPerParticleParameter('x0')
    pos_restraint.addPerParticleParameter('y0')
    pos_restraint.addPerParticleParameter('z0')

    for i, atom in enumerate(structure.positions):
        if structure.atoms[i].name == 'DUM':
            pos_restraint.addParticle(i,atom.value_in_unit(unit.nanometers))

    static_restraints = []
    
    # Create Distance Restraint
    static_distance_rest = [D[0], H[0]]
    static_init_dist = pt.distance(traj, D[0]+' '+H[0])[0]

    dist_restraint = mm.CustomBondForce('k*(r-r0)^2')
    dist_restraint.addPerBondParameter('k')
    dist_restraint.addPerBondParameter('r0')

    r0 = static_init_dist * unit.angstroms
    k = 100 * unit.kilocalories_per_mole / unit.angstroms**2

    dist_restraint.addBond(D_i[0], H_i[0], [k, r0])
    static_restraints.append(dist_restraint)

    # Create Angle Restraint 1
    static_angle_rest_1 = [D[1], D[0], H[0]]
    static_init_angle_1 = pt.angle(traj, D[1]+' '+D[0]+' '+H[0])[0]

    angle_restraint_1 = mm.CustomAngleForce('0.5*k*(theta-theta0)^2')
    angle_restraint_1.addPerAngleParameter('k')
    angle_restraint_1.addPerAngleParameter('theta0')

    theta0 = static_init_angle_1 * unit.degrees
    k = 100 * unit.kilocalories_per_mole / unit.radians**2

    angle_restraint_1.addAngle(D_i[1], D_i[0], H_i[0], [k, theta0])
    static_restraints.append(angle_restraint_1)

    # Create Dihedral Restraint 1
    static_dihedral_rest_1 = [D[2], D[1], D[0], H[0]]
    static_init_dihedral_1 = pt.dihedral(traj, D[2]+' '+D[1]+' '+D[0]+' '+H[0])[0]

    dihedral_restraint_1 = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
    dihedral_restraint_1.addPerTorsionParameter('k')
    dihedral_restraint_1.addPerTorsionParameter('theta0')

    theta0 = static_init_dihedral_1 * unit.degrees
    k = 100 * unit.kilocalories_per_mole / unit.radians**2

    dihedral_restraint_1.addTorsion(D_i[2], D_i[1], D_i[0], H_i[0], [k, theta0])
    static_restraints.append(dihedral_restraint_1)

    # angle/dihedral restraints in unit.kilocalories+per_mole / unit.radians**2
    #system.addForce(pos_restraint)

    if num_rests > 0:
        for rest in static_restraints[0:num_rests]:
            system.addForce(rest)

    simulation = app.Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('CPU'))
    simulation.context.setPositions(inpcrd.positions)

    simulation.reporters.append(NetCDFReporter(path+traj_out, 250))
    simulation.reporters.append(
        app.StateDataReporter(
            path+md_out, 250,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True
        )
    )

    simulation.step(settings['numsteps'])

test_openmm_cb6but_sim(num_rests=0)
