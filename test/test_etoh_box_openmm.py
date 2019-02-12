import parmed as pmd
from simtk.openmm import *
import simtk.openmm.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app
from mdtraj.reporters import NetCDFReporter

def test_openmm_etoh_sim():
    path = './etoh_test/sim_openmm/'
    topology = 'etoh.prmtop'
    coordinates = 'etoh.rst7'
    md_out = 'etoh_openmm.csv'

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

    simulation.reporters.append(NetCDFReporter(path+'etoh_openmm.nc', 250))
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
    
    simulation.step(100000)

def printEnergy(md_out):
    with open(md_out, 'r') as f:
        data = [line.strip().split(',') for line in f]
    for item in data[1:]:
        energy = float(item[4])*unit.kilojoules_per_mole
        energy = energy.in_units_of(unit.kilocalories_per_mole)
        print('Total Energy:', energy)

def plotEnergy(md_out):
    import matplotlib.pyplot as plt
    import numpy as np
    time = []
    energy = []
    with open(md_out, 'r') as f:
        data = [line.strip().split(',') for line in f]
    for item in data[1:]:
        time.append(round(float(item[1]),1))
        e = float(item[4])*unit.kilojoules_per_mole
        e = e.in_units_of(unit.kilocalories_per_mole)/unit.kilocalories_per_mole
        energy.append(e)
    x = np.array(time)
    y = np.array(energy)
    fig, ax = plt.subplots()
    ax.set_xlabel('time (picoseconds)')
    ax.set_ylabel('energy (kcal/mol)')
    ax.plot(x,y)
    plt.show()

path = './etoh_test/sim_openmm/'
md_out = 'etoh_openmm.csv'

#test_openmm_etoh_sim()
#printEnergy(path+md_out)
plotEnergy(path+md_out)
