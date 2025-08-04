import hoomd
import gsd.hoomd
import torch
from sim_setup.create_simulation import simulation
from sim_setup.place_mpcd_particles import create_snapshot_gsd_mpcd
from sim_setup.parse_args import (
    task_name,
    collision_period,
    n_sim_steps,
    parallel_plate_speeds,
    init_file_name,
    init_frame_number,
    sidechains_can_bond,
    save_every,
    prod_save_folder,
    reversible_bonds,
)
from sim_setup.make_bonds_custom_action import FormBonds
import time
import os

# ensure that the init_file_name and init_frame_number are defined
if init_file_name == "" or init_frame_number == "":
    raise ValueError("Both init_file_name and init_frame_number must be defined")

start_time = time.time()

# load state from the saved GSD conformations
frames = gsd.hoomd.open(name=init_file_name, mode="r")
frame_number = init_frame_number  # pick which frame to use for the initial state
frame = frames[frame_number]

# create the snapshot with MPCD particles and set the simulation state
snapshot, filled_height, _, _ = create_snapshot_gsd_mpcd(frame)
simulation.create_state_from_snapshot(snapshot)

# use a GPU to run the simulation, if available
has_cuda = torch.cuda.is_available()
if has_cuda:
    torch.set_default_device("cuda")


typeids = torch.tensor(snapshot.particles.typeid)
n_sidechain = torch.sum(typeids == 1)
n_backbone = torch.sum(typeids == 0)


backbone_indices = torch.where(typeids == 0)[0]
sidechain_indices = torch.where(typeids == 1)[0]
bonds = torch.tensor(snapshot.bonds.group, dtype=torch.int)  # Convert bonds to a tensor

# Check if the first particle in the bond is a backbone particle and the second is a sidechain particle
backbone_to_sidechain_mask = torch.isin(bonds[:, 0], backbone_indices) & torch.isin(
    bonds[:, 1], sidechain_indices
)

# Check if the second particle in the bond is a backbone particle and the first is a sidechain particle
sidechain_to_backbone_mask = torch.isin(bonds[:, 1], backbone_indices) & torch.isin(
    bonds[:, 0], sidechain_indices
)

# Get the particle pairs where the masks are true
backbone_to_sidechain_pairs = bonds[backbone_to_sidechain_mask]
sidechain_to_backbone_pairs = bonds[sidechain_to_backbone_mask]

# Create dictionary-like mapping using tensors
backbone_partners1 = dict(
    zip(
        backbone_to_sidechain_pairs[:, 0].tolist(),
        backbone_to_sidechain_pairs[:, 1].tolist(),
    )
)
backbone_partners2 = dict(
    zip(
        sidechain_to_backbone_pairs[:, 1].tolist(),
        sidechain_to_backbone_pairs[:, 0].tolist(),
    )
)

# Merge both dictionaries
backbone_partners = {**backbone_partners1, **backbone_partners2}

# Convert backbone partners to a PyTorch tensor
backbone_partners = list(backbone_partners.keys())
backbone_partners = torch.tensor(backbone_partners, dtype=torch.int)

# this updater will form bonds between side chain particles under certain conditions
form_bonds = FormBonds(
    backbone_partners=backbone_partners,
    has_cuda=has_cuda,
    can_break_bonds=reversible_bonds,
)
check_bond_formation_every = (
    5 * collision_period
)  # lower value = more accurate, but slower to run
bond_form_updater = hoomd.update.CustomUpdater(
    action=form_bonds, trigger=hoomd.trigger.Periodic(check_bond_formation_every)
)
if sidechains_can_bond:
    simulation.operations.updaters.append(bond_form_updater)

# run the simulation
# store bonds and angles in the GSD file
results_filename = f"production_{task_name}"
output_path = f"{results_filename}.gsd"
if os.path.exists(output_path):
    raise FileExistsError(
        f"Output file {output_path} already exists. Aborting to prevent overwrite."
    )

os.makedirs(f"{prod_save_folder}", exist_ok=True)
gsd_writer_prod = hoomd.write.GSD(
    trigger=save_every,
    filename=output_path,
    dynamic=["property", "topology", "momentum", "attribute"],
)


# add a logger to track thermo properties and ensure energy conservation
thermo_peroperties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
simulation.operations.computes.append(thermo_peroperties)
logger = hoomd.logging.Logger()
logger.add(thermo_peroperties)

simulation.operations.writers.append(gsd_writer_prod)
gsd_writer_prod.logger = logger  # log thermo properties to the GSD file

parallel_plate_speeds = torch.tensor(parallel_plate_speeds)

# reactive stage
# for effiency, only update bonds or plate speeds during MPCD collisions
n_collisions = n_sim_steps // collision_period

print("Running reactive stage...")
simulation.run(n_sim_steps)

print("Simulation complete")
gsd_writer_prod.flush()

end_time = time.time()

# end_time - start_time is the total runtime in seconds
total_seconds = end_time - start_time

# Convert to hours, minutes, and seconds
hours, remainder = divmod(total_seconds, 3600)
minutes, seconds = divmod(remainder, 60)

print(f"Runtime: {int(hours):02}:{int(minutes):02}:{int(seconds):02}")
