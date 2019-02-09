import parmed as pmd
from paprika.openmm import *
from simtk.openmm import *
import simtk.openmm.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app
from mdtraj.reporters import NetCDFReporter
from sys import stdout

def test_openmm_etoh_sim():
    path = './etoh_test/sim_openmm/'
    topology = 'etoh.prmtop'
    coordinates = 'etoh.rst7'

    prmtop = app.AmberPrmtopFile(path+topology)
    inpcrd = app.AmberInpcrdFile(path+coordinates)

    settings = {
        'nonbonded_method': app.PME,
        'constraints': app.HBonds,
        'temperature': 298.15*unit.kelvin,
        'friction': 1/unit.picosecond,
        'timestep': 0.002*unit.picosecond
    }

    system = prmtop.createSystem(
        nonbondedMethod = settings['nonbonded_method'],
        constraints = settings['constraints']
    )

    integrator = LangevinIntegrator(
        settings['temperature'],
        settings['friction'],
        settings['timestep']
    )

    simulation = app.Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('Reference'))
    simulation.context.setPositions(inpcrd.positions)
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    simulation.reporters.append(NetCDFReporter('etoh_openmm.nc', 500))
    simulation.reporters.append(
        app.StateDataReporter(
            stdout, 500,
            step=True,
            potentialEnergy=True,
            temperature=True,
            volume=True,
            density=True
        )
    )
    simulation.step(5000)

test_openmm_etoh_sim()
