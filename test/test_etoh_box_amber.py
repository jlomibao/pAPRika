from paprika import amber

def test_amber_etoh_sim():
    path = './etoh_test/sim_amber/'
    topology = 'etoh.prmtop'
    coordinates = 'etoh.rst7'

    sim = amber.Simulation()
    sim.config_pbc_md()
    sim.path = path
    sim.restraint_file = 'no_rest.in'
    sim.executable = 'pmemd'
    sim.prefix = 'etoh_amber'
    sim.inpcrd = coordinates
    sim.ref = coordinates
    sim.topology = topology

    sim.cntrl['ntb'] = 2
    sim.cntrl['nstlim'] = 100000
    sim.cntrl['ntwx'] = 250
    sim.cntrl['ntpr'] = 250
    sim.cntrl['ntwprt'] = 0
    sim.cntrl['cut'] = 9.0
    sim.cntrl['dt'] = 0.002

    sim.run()

def plotEnergy(md_out):
    import re
    import matplotlib.pyplot as plt
    import numpy as np
    time = []
    energy = []
    with open(md_out, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if re.search('NSTEP',line):
                cols = line.strip().split()
                time.append(round(float(cols[5]),1))
            if re.search('Etot',line):
                cols = line.strip().split()
                energy.append(float(cols[2]))
    x = np.array(time[:-2])
    y = np.array(energy[:-2])
    fig, ax = plt.subplots()
    ax.set_xlabel('time (picoseconds)')
    ax.set_ylabel('energy (kcal/mol)')
    ax.plot(x,y)
    plt.show()

path = './etoh_test/sim_amber/'
md_out = 'etoh_amber.out'

#test_amber_etoh_sim()
#plotEnergy(path+md_out)
