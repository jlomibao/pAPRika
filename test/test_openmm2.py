from paprika import tleap
from paprika import align
from paprika import dummy
from paprika import restraints
from paprika.restraints import static_DAT_restraint
from paprika.restraints import DAT_restraint
from paprika.restraints import create_window_list
from paprika.restraints import amber_restraint_line
import parmed as pmd
import numpy as np
import pytraj as pt
import shutil
import os
import sys

import logging
logging.basicConfig(filename='test_openmm.log',level=logging.INFO)
logging.getLogger().setLevel(logging.DEBUG)
logging.info('Started logging...')

struct_dir = "cb6-but"
host_frcmod = "cb6.frcmod"
host_mol2 = "cb6.mol2"
guest_mol2 = "but.mol2"
pdb_in = "cb6-but.pdb"
pdb_out = "vac.pdb"
topology_out = "vac.prmtop"
coordinates_out = "vac.rst7"

# Use tleap to create topology and coordinate files

system = tleap.System()
system.output_path = struct_dir
system.pbc_type = None
system.neutralize = False
system.template_lines = [
  "source leaprc.gaff",
  "loadamberparams "+host_frcmod,
  "CB6 = loadmol2 "+host_mol2,
  "BUT = loadmol2 "+guest_mol2,
  "model = loadpdb "+pdb_in,
  "check model",
  "savepdb model "+pdb_out,
  "saveamberparm model "+topology_out+" "+coordinates_out
]

system.build()

# Load in structure via topology and coordinate files output

structure = pmd.load_file(struct_dir+'/'+topology_out,
                          struct_dir+'/'+coordinates_out,
                          structure=True)
host = ":CB6"
guest = ":BUT"

H = [host+"@C7", host+"@C31", host+"@C19"]
G = [guest+"@C", guest+"@C3"]
C1 = [host+"@O", host+"@O6"]
C2 = [host+"@O2", host+"@O8"]
C3 = [host+"@O4", host+"@O10"]
D = [":DM1", ":DM2", ":DM3"]

# Move atoms to align to the z-axis along where the guest will be pulled

aligned_top = "aligned.prmtop"
aligned_coord = "aligned.rst7"

aligned_structure = align.zalign(structure,G[0],G[1])
aligned_structure.save(struct_dir+'/'+aligned_top, overwrite=True)
aligned_structure.save(struct_dir+'/'+aligned_coord, overwrite=True)

# Add dummy atoms

structure = pmd.load_file(struct_dir+'/'+aligned_top,
                          struct_dir+'/'+aligned_coord,
                          structure=True)

structure = dummy.add_dummy(structure, residue_name="DM1", z=-6.0)
structure = dummy.add_dummy(structure, residue_name="DM2", z=-9.0)
structure = dummy.add_dummy(structure, residue_name="DM3", y=2.2, z=-11.2)

dum_top = "aligned_with_dummy.prmtop"
dum_coord = "aligned_with_dummy.rst7"
dum_pdb = "aligned_with_dummy.pdb"

structure.save(struct_dir+'/'+dum_top, overwrite=True)
structure.save(struct_dir+'/'+dum_coord, overwrite=True)
structure.save(struct_dir+'/'+dum_pdb, overwrite=True)

# Write out frcmod and mol2 files for dummy atoms

dum_frcmod = "dummy.frcmod"
dm1_mol2 = "dm1.mol2"
dm2_mol2 = "dm2.mol2"
dm3_mol2 = "dm3.mol2"

dummy.write_dummy_frcmod(filepath=struct_dir+'/'+dum_frcmod)
dummy.write_dummy_mol2(residue_name="DM1", filepath=struct_dir+'/'+dm1_mol2)
dummy.write_dummy_mol2(residue_name="DM2", filepath=struct_dir+'/'+dm2_mol2)
dummy.write_dummy_mol2(residue_name="DM3", filepath=struct_dir+'/'+dm3_mol2)

# Add parameters for dummy atoms and output new topology and coordinate files
file_prefix = "cb6-but-dum"
pdb_final = file_prefix+".pdb"
top_final = file_prefix+".prmtop"
coord_final = file_prefix+".rst7"

system = tleap.System()
system.output_path = struct_dir
system.pbc_type = None
system.neutralize = False

system.template_lines = [
  "source leaprc.gaff",
  "loadamberparams "+host_frcmod,
  "loadamberparams "+dum_frcmod,
  "CB6 = loadmol2 "+host_mol2,
  "BUT = loadmol2 "+guest_mol2,
  "DM1 = loadmol2 "+dm1_mol2,
  "DM2 = loadmol2 "+dm2_mol2,
  "DM3 = loadmol2 "+dm3_mol2,
  "model = loadpdb "+dum_pdb,
  "check model",
  "savepdb model "+pdb_final,
  "saveamberparm model "+top_final+" "+coord_final
]

system.build()


system.build()

structure = pmd.load_file(struct_dir+'/'+top_final, struct_dir+'/'+coord_final)

traj = pt.load(struct_dir+'/'+coord_final,struct_dir+'/'+top_final)
guest_init_dist = pt.distance(traj,D[0]+' '+G[0])[0]
conf_init_dist1 = pt.distance(traj,C1[0]+' '+C1[1])[0]
conf_init_dist2 = pt.distance(traj,C2[0]+' '+C2[1])[0]
conf_init_dist3 = pt.distance(traj,C3[0]+' '+C3[1])[0]

attach_string = "0.00 0.40 0.80 1.60 2.40 4.00 5.50 8.65 11.80 18.10 24.40 37.00 49.60 74.80 100.00"

attach_fractions = [float(i) / 100 for i in attach_string.split()]

pull_distances = np.arange(0.0 + guest_init_dist, 18.0 + guest_init_dist, 1.0)

release_fractions = attach_fractions[::-1]

windows = [len(attach_fractions), len(pull_distances), len(release_fractions)]

# Static Restraints
static_distance_rest = [D[0], H[0]]
#static_angle_rest1 = [D[1], D[0], H[0]]
#static_dihedral_rest1 = [D[2], D[1], D[0], H[0]]
#static_angle_rest2 = [D[0], H[0], H[1]]
#static_dihedral_rest2 = [D[1], D[0], H[0], H[1]]
#static_dihedral_rest3 = [D[0], H[0], H[1], H[2]]
static_restraint_atoms = [
    static_distance_rest,
    #static_angle_rest1,
    #static_dihedral_rest1,
    #static_angle_rest2,
    #static_dihedral_rest2,
    #static_dihedral_rest3
]
static_restraint_distance_fc = 5.0
static_restraint_angle_fc = 100.0

# Guest Restraints
guest_distance_rest = [D[0], G[0]]
guest_angle_rest1 = [D[1], D[0], G[0]]
guest_angle_rest2 = [D[0], G[0], G[1]]
guest_restraint_atoms = [
    guest_distance_rest,
    guest_angle_rest1,
    guest_angle_rest2
]
guest_restraint_targets = [guest_init_dist, 180.0, 180.0]
guest_restraint_target_final = [18.0+guest_init_dist, 180.0, 180.0]
guest_restraint_distance_fc = 5.0
guest_restraint_angle_fc = 100.0

# Create restraint objects
static_restraints = []
for index, atoms in enumerate(static_restraint_atoms):
    this = static_DAT_restraint(restraint_mask_list=atoms,
                               num_window_list=windows,
                               ref_structure=structure,
                               force_constant=static_restraint_angle_fc if len(atoms)>2 else static_restraint_distance_fc,
                               amber_index=True)
    static_restraints.append(this)

guest_restraints = []
for index, atoms in enumerate(guest_restraint_atoms):
    if len(atoms) > 2:
        angle = True
    else:
        angle = False
    this = DAT_restraint()
    this.auto_apr = True
    this.amber_index = True
    this.topology = structure
    this.mask1 = atoms[0]
    this.mask2 = atoms[1]
    if angle:
        this.mask3 = atoms[2]
        this.attach['fc_final'] = guest_restraint_angle_fc
        this.release['fc_final'] = guest_restraint_angle_fc
    else:
        this.attach['fc_final'] = guest_restraint_distance_fc
        this.release['fc_final'] = guest_restraint_distance_fc
    this.attach['target'] = guest_restraint_targets[index]
    this.attach['fraction_list'] = attach_fractions

    this.pull['target_final'] = guest_restraint_target_final[index]
    this.pull['num_windows'] = windows[1]

    this.release['target'] = guest_restraint_targets[index]
    this.release['fraction_list'] = release_fractions

    this.initialize()
    guest_restraints.append(this)

restraints = static_restraints

for i, atom in enumerate(structure.atoms):
    if atom.name == 'DUM':
        atom.mass=0.0

from simtk.openmm.app import *
from simtk.openmm import *
import simtk.openmm.openmm as mm
from simtk.unit import *
from sys import stdout
from mdtraj.reporters import NetCDFReporter

topology = "cb6-but-0m_dum.prmtop"
coordinates = "cb6-but-0m_dum.rst7"
structure.save(struct_dir+'/'+topology, overwrite=True)
structure.save(struct_dir+'/'+coordinates, overwrite=True)

prmtop = AmberPrmtopFile(struct_dir+'/'+topology)
inpcrd = AmberInpcrdFile(struct_dir+'/'+coordinates)

system = prmtop.createSystem(
    nonbondedMethod=NoCutoff,
    constraints=HBonds,
    implicitSolventSaltConc=0.1*moles/liter #may want to use 0.0
)


integrator = LangevinIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds
)

from paprika.openmm import *

opmm_sim = OpenMM_GB_simulation()

OpenMM_GB_simulation.add_openmm_restraints(opmm_sim,system,restraints,"attach",1)

simulation = Simulation(prmtop.topology, system, integrator, mm.Platform.getPlatformByName('Reference'))
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter('cb6-but_output.pdb', 500))
simulation.reporters.append(NetCDFReporter('paprika_openmm_test.nc', 500))
simulation.reporters.append(
    StateDataReporter(
        stdout, 500,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True
    )
)

simulation.step(1000)
