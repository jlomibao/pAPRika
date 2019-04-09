import os

import paprika
from paprika import amber

import logging

logging.getLogger().setLevel(logging.DEBUG)
logging.info('Started logging...')

path = "./windows_apr/"
file_prefix = "cb6-but-dum"

num_steps = 1000000

window_list = os.listdir(path)

for window in window_list:
    sim = amber.Simulation()

    # Production
    sim.executable = "pmemd.cuda"
    sim.path = path + window
    sim.prefix = "amber_prod"
    sim.inpcrd = "amber_equil.rst7"
    sim.ref = file_prefix + ".rst7"
    sim.topology = file_prefix + ".prmtop"
    sim.restraint_file = "rest.in"

    sim.cntrl['imin'] = 0
    sim.cntrl['maxcyc'] = 0
    sim.cntrl['ncyc'] = 0
    sim.cntrl['nstlim'] = num_steps
    sim.cntrl['ntwx'] = 250
    sim.cntrl['ntpr'] = 250
    sim.cntrl['ntwprt'] = 0
    sim.cntrl['ntf'] = 1
    sim.cntrl['nmropt'] = 1
    sim.cntrl['irest'] = 0
    sim.cntrl['ntb'] = 0
    sim.cntrl['cut'] = 999.0
    sim.cntrl['igb'] = 1
    sim.cntrl['dt'] = 0.002
    sim.cntrl['ntr'] = 1
    sim.cntrl['ntp'] = 0
    sim.cntrl['ntc'] = 2
    sim.cntrl['restraint_wt'] = 50.0
    sim.cntrl['restraintmask'] = "'@DUM'"

    sim.run(fail_ok=False)

    # Equilibration
    sim.executable = "pmemd.cuda"
    sim.topology = file_prefix + ".prmtop"
    sim.path = path + window
    sim.restraint_file = "rest.in"
    sim.prefix = "amber_equil"
    sim.inpcrd = "minimize.rst7"
    sim.ref = file_prefix + ".rst7"

    if not os.path.isfile(path + window + "/amber_prod.rst7"):
        sim.run(fail_ok=False)
