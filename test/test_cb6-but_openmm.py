from simtk.openmm import *
import simtk.openmm.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app
from mdtraj.reporters import NetCDFReporter

# Need to change simulation settings
def test_openmm_cb6but_sim():
    path = './cb6-but_test/guest_inside/'
    topology = 'cb6-but-dum.prmtop'
    coordinates = 'cb6-but-dum.rst7'
    md_out = 'cb6-but_openmm.csv'

    prmtop = app.AmberPrmtopFile(path+topology)
    inpcrd = app.AmberInpcrdFile(path+coordinates)

    settings = {
        'nonbonded_method': app.PME,
        'nonbonded_cutoff': 9*unit.angstrom,
        'constraints': app.HBonds,
        'temperature': 298.15*unit.kelvin,
        'friction': 1/unit.picosecond,
        'timestep': 0.002*unit.picosecond
    }

    system = prmtop.createSystem(
        nonbondedMethod = settings['nonbonded_method'],
        nonbondedCutoff = settings['nonbonded_cutoff'],
        constraints = settings['constraints']
    )
    barostat = mm.MonteCarloBarostat(1.0*unit.bar,298.15*unit.kelvin,25)
    system.addForce(barostat)

    integrator = LangevinIntegrator(
        settings['temperature'],
        settings['friction'],
        settings['timestep']
    )

    simulation = app.Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('CPU'))
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.reporters.append(NetCDFReporter(path+'cb6_but_openmm.nc', 250))
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

    simulation.step(1000)

