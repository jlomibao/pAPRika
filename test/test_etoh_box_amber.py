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
    sim.cntrl['nstlim'] = 10000
    sim.cntrl['ntwx'] = 250
    sim.cntrl['ntwprt'] = 0
    sim.cntrl['cut'] = 9.0
    sim.cntrl['dt'] = 0.002

    sim.run()



test_amber_etoh_sim()
