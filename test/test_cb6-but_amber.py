from paprika import amber
from paprika import restraints
from paprika.restraints import static_DAT_restraint
import parmed as pmd

def setup_restraints(num_rests=6):
    path = './cb6-but_test/guest_inside/'
    topology = 'cb6-but-dum.prmtop'
    coordinates = 'cb6-but-dum.rst7'

    structure = pmd.load_file(path+topology, path+coordinates, structure=True)

    host = ":CB6"
    guest = ":BUT"

    H = [host+"@C7", host+"@C31", host+"@C19"]
    G = [guest+"@C", guest+"@C3"]
    D = [":DM1", ":DM2", ":DM3"]

    attach_fractions = [1.0]
    windows = [len(attach_fractions), 1, 0]

    static_distance_rest = [D[0], H[0]]
    static_angle_rest_1 = [D[1], D[0], H[0]]
    static_dihedral_rest_1 = [D[2], D[1], D[0], H[0]]
    static_angle_rest_2 = [D[0], H[0], H[1]]
    static_dihedral_rest_2 = [D[1], D[0], H[0], H[1]]
    static_dihedral_rest_3 = [D[0], H[0], H[1], H[2]]

    static_restraint_atoms = [
        static_distance_rest,
        static_angle_rest_1,
        static_dihedral_rest_1,
        static_angle_rest_2,
        static_dihedral_rest_2,
        static_dihedral_rest_3
    ]

    static_restraint_distance_fc = 5.0
    static_restraint_angle_fc = 100.0

    static_restraints = []
    for index, atoms in enumerate(static_restraint_atoms[0:num_rests]):
        this = static_DAT_restraint(
            restraint_mask_list=atoms,
            num_window_list=windows,
            ref_structure=structure,
            force_constant=static_restraint_angle_fc if len(atoms)>2 else static_restraint_distance_fc,
            amber_index=True
        )
        static_restraints.append(this)

    return static_restraints

def test_amber_cb6but_sim(num_rests=0):
    path = './cb6-but_test/guest_inside/'
    topology = 'cb6-but-dum.prmtop'
    coordinates = 'cb6-but-dum.rst7'

    if num_rests > 0:
        rests = setup_restraints(num_rests)
        with open(path+'restraints.{:02d}.in'.format(num_rests), 'w') as f:
            for rest in rests:
                string = restraints.amber_restraint_line(rest,'a000')
                if string is not None:
                    f.write(string)

    sim = amber.Simulation()
    sim.path = path
    if num_rests > 0:
        sim.prefix = 'cb6_but_amber_rest_{:02d}'.format(num_rests)
        sim.restraint_file = 'restraints.{:02d}.in'.format(num_rests)
        sim.cntrl['nmropt'] = 1 # NMR analysis
        sim.cntrl['ntr'] = 1 # Flag for restraining specified atoms
        sim.ref = coordinates
        sim.cntrl['restraint_wt'] = 50.0
        sim.cntrl['restraintmask'] = "'@DUM'"
    else:
        sim.prefix = 'cb6_but_amber'
        sim.cntrl['nmropt'] = 0 # No NMR analysis
    sim.executable = 'pmemd'
    sim.inpcrd = coordinates
    sim.topology = topology

    sim.cntrl['barostat'] = 2 # Monte Carlo barostat
    sim.cntrl['ntb'] = 0 # no periodic boundary
    sim.cntrl['ntp'] = 0 # no pressure scaling
    sim.cntrl['ntf'] = 1 # Complete interaction is calculated
    sim.cntrl['ntc'] = 1 # SHAKE is not performed
    sim.cntrl['nstlim'] = 100000 # number of steps
    sim.cntrl['ntwx'] = 250 # write to coordinate file every ntwx steps
    sim.cntrl['ntpr'] = 250 # print to mdout every ntpr steps
    sim.cntrl['ntwprt'] = 0 # include all atoms of system when writing traj
    sim.cntrl['cut'] = 999.0 # nonbonded cutoff in angstroms
    sim.cntrl['igb'] = 1 # Flag for Born implicit solvent model
    sim.cntrl['dt'] = 0.002 # timestep in picoseconds

    sim.run()

test_amber_cb6but_sim(num_rests=1)
