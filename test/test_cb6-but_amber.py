from paprika import amber

def test_amber_cb6but_sim():
    path = './cb6-but_test/guest_inside/'
    topology = 'cb6-but-dum.prmtop'
    coordinates = 'cb6-but-dum.rst7'

    sim = amber.Simulation()
    sim.path = path
    sim.executable = 'pmemd'
    sim.prefix = 'cb6_but_amber'
    sim.inpcrd = coordinates
    sim.topology = topology

    sim.cntrl['nmropt'] = 0 # No NMR analysis
    sim.cntrl['barostat'] = 0 # Monte Carlo barostat
    sim.cntrl['ntb'] = 0 # no periodic boundary
    sim.cntrl['ntp'] = 0 # no pressure scaling
    sim.cntrl['ntf'] = 2 # bond interactions involving H-atoms omitted
    sim.cntrl['ntc'] = 2 # bonds involving hydrogen constrained
    sim.cntrl['nstlim'] = 1000 # number of steps
    sim.cntrl['ntwx'] = 250 # write to coordinate file every ntwx steps
    sim.cntrl['ntpr'] = 250 # print to mdout every ntpr steps
    sim.cntrl['ntwprt'] = 0 # include all atoms of system when writing traj
    sim.cntrl['cut'] = 999.0 # nonbonded cutoff in angstroms
    sim.cntrl['igb'] = 1 
    sim.cntrl['dt'] = 0.002 # timestep in picoseconds

    sim.run()

test_amber_cb6but_sim()
