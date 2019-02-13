from paprika import amber

def test_amber_etoh_sim():
    path = './etoh_test/sim_amber/'
    topology = 'etoh.prmtop'
    coordinates = 'etoh.rst7'

    sim = amber.Simulation()
    sim.config_pbc_md()
    sim.path = path
    #sim.restraint_file = 'no_rest.in'
    sim.executable = 'pmemd'
    sim.prefix = 'etoh_amber'
    sim.inpcrd = coordinates
    sim.ref = coordinates
    sim.topology = topology

    sim.cntrl['nmropt'] = 0 # No NMR analysis
    sim.cntrl['barostat'] = 2 # Monte Carlo barostat
    sim.cntrl['ntb'] = 2 # periodic boundary / constant pressure
    sim.cntrl['ntf'] = 2 # bond interactions involving H-atoms omitted
    sim.cntrl['ntc'] = 2 # bonds involving hydrogen constrained
    sim.cntrl['ntt'] = 3 # Langevin Dynamics with collision freq of gamma_ln
    sim.cntrl['gamma_ln'] = 1.0 # Collision frequency in ps^-1
    sim.cntrl['nstlim'] = 100000 # number of steps
    sim.cntrl['ntwx'] = 250 # write to coordinate file every ntwx steps
    sim.cntrl['ntpr'] = 250 # print to mdout every ntpr steps
    sim.cntrl['ntwprt'] = 0 # include all atoms of system when writing traj
    sim.cntrl['cut'] = 9.0 # nonbonded cutoff in angstroms
    sim.cntrl['dt'] = 0.002 # timestep in picoseconds

    sim.run()

def plotEnergy(md_out,hist=False):
    import re
    import matplotlib.pyplot as plt
    import numpy as np
    time = []
    energy = []
    count1 = 0
    count2 = 0
    with open(md_out, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if re.search('NSTEP',line):
                count1 +=1
                cols = line.strip().split()
                time.append(round(float(cols[5]),1))
            if re.search('Etot',line):
                count2 +=1
                cols = line.strip().split()
                energy.append(float(cols[2]))
                if (count2!=count1):
                    print(count2,count1,time[-1],energy[-1])

    x = np.array(time[:-2])
    y = np.array(energy[:-2])
    
    fig, ax = plt.subplots()
    ax.set_xlabel('time (picoseconds)')
    ax.set_ylabel('energy (kcal/mol)')
    ax.plot(x,y)
    plt.show()
    
    if hist == True:
        plt.hist(y, 50, alpha=1.0)
        plt.xlabel('energy kcal/mol')
        plt.ylabel('count')
        plt.show()

path = './etoh_test/sim_amber/'
md_out = 'etoh_amber.out'

test_amber_etoh_sim()
#plotEnergy(path+md_out)
