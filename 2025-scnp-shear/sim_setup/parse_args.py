from torch import pi, set_default_device
import argparse

"""
This script handles the setup of the simulation. It parses the command line arguments and sets up the simulation parameters.
"""

# allow specifying number of backbone and sidechain beads, as well as sidechain locations, from the command line
parser = argparse.ArgumentParser(
    description="Process backbone and sidechain parameters."
)

parser.add_argument(
    "--n_backbone", type=int, default=12, help="Number of backbone elements"
)
parser.add_argument(
    "--backbone_partners",
    type=str,
    default="",
    help="List of index locations of the side chains along the backbone. Default to no side chains. Example: --backbone_partners 3 6 9",
)
parser.add_argument("--fene_k", type=float, default=30, help="FENE bond constant")
parser.add_argument("--fene_r0", type=float, default=1.5, help="FENE bond length")
parser.add_argument("--wca_sigma", type=float, default=1, help="WCA sigma parameter")
parser.add_argument(
    "--wca_epsilon", type=float, default=1, help="WCA epsilon parameter"
)

# https://pubs.acs.org/doi/full/10.1021/acs.macromol.7b01876 investigated kappa between 0 and 50
parser.add_argument(
    "--angle_stiffness",
    type=float,
    default=5,
    help="Bending potential constant. Default to 5.",
)

parser.add_argument(
    "--sidechains_can_bond",
    action="store_true",
    default=False,
    help="Whether sidechain particles can form bonds with each other. Default to False",
)

parser.add_argument(
    "--sidechain_sidechain_fene_k",
    type=float,
    default=30,
    help="FENE bond constant for sidechain-sidechain bonds. Default to 30, the same as the backbone FENE k",
)

# max bond length must be much larger than the capture radius to prevent the sudden force from blowing up the integrator
parser.add_argument(
    "--sidechain_sidechain_fene_r0",
    type=float,
    default=1.5,
    help="Fene r0 between sidechain particles before the bond breaks, as a multiple of wca_sigma. Default to 1.5, the same as the backbone FENE r0",
)

parser.add_argument(
    "--sidechain_equilibrium_bonding_angle",
    type=float,
    default=pi,
    help="Equilibrium bonding angle between a sidechain particle and its associated backbone particle",
)

parser.add_argument(
    "--sidechain_capture_radius",
    type=float,
    default=1.3,
    help="Capture radius for sidechain bonding, as a multiple of wca_sigma",
)

parser.add_argument(
    "--sidechain_sidechain_bond_break_length",
    type=float,
    default=float("inf"),
    help="Length at which a sidechain-sidechain bond will break, as a multiple of wca_sigma. Default is that sidechain-sidechain bonds will never break due to length deviation",
)

parser.add_argument(
    "--sidechain_sidechain_bond_break_angle_tolerance",
    type=float,
    default=pi,
    help="Maximum deviation from the equilibrium angle, in radians, before a sidechain-sidechain bond will break. Default is that sidechain-sidechain bonds will never break due to angle deviation. Must be greater than the capture tolerance, --sidechain_bonding_angle_tolerance",
)


parser.add_argument(
    "--sidechain_bonding_angle_tolerance",
    type=float,
    default=pi / 6,
    help="Tolerance for sidechain bonding angle for initial bond formation, in radians. Default to pi/6 (30 degrees)",
)

parser.add_argument(
    "--task_name",
    type=str,
    default=0,
    help="Task number for the simulation (used for logging purposes)",
)

parser.add_argument(
    "--collision_period",
    type=int,
    default=20,
    help="Number of streaming steps between MPCD collisions",
)

parser.add_argument(
    "--n_sim_steps",
    type=int,
    default=int(1e6),
    help="Number of simulation steps",
)

parser.add_argument(
    "--save_every",
    type=int,
    default=10_000,
    help="Number of steps between saving simulation state",
)

parser.add_argument(
    "--production_save_folder",
    type=str,
    default="production",
    help="Folder to save production data",
)

parser.add_argument(
    "--reversible_bonds",
    action="store_true",
    default=False,
    help="Whether pendant particle bonding is reversible in the simulation. Default to False",
)

parser.add_argument(
    "--n_conformation_gen_steps",
    type=int,
    default=100_000,
    help="If running a script to generate seed conformations, specifies the number of total conformation steps to run",
)


# parallel plate speeds must be defined from a file, since it's a long list that can't be directly passed as a command line argument
def read_floats_from_file(filename):
    if "constant" in filename:
        return [float(filename.split("=")[1])]
    else:
        with open(filename, "r") as file:
            return [float(line.strip()) for line in file]


parser.add_argument(
    "--parallel_plate_speeds",
    type=str,
    default="constant=0",
    help="Path to a txt file containing parallel plate speeds at each simulation step. This should have length n_sim_steps/collision_period, since updates occur after each collision. If a constant speed is desired, use 'constant=X'. The speed should never exceed 1.0, since v_max in each direction should be less than 0.5 as to not exceed the speed of sound. https://github.com/mphowardlab/azplugins/blob/v0.12.0/doc/source/tutorial/01_reverse_perturbation/01_reverse_perturbation_2_mpcd.ipynb",
)

parser.add_argument(
    "--init_file_name",
    type=str,
    default="",
    help="Path to the state GSD file to load. Required for production.",
)

parser.add_argument(
    "--init_frame_number",
    type=int,
    default=0,
    help="Frame number to load from the state GSD file. Required for production.",
)

parser.add_argument(
    "--init_conformations_folder",
    type=str,
    default="conformations_during_init",
    help="Folder to save initial conformations",
)

parser.add_argument(
    "--bath_temp",
    type=float,
    default=1.0,
    help="Temperature of the MPCD fluid bath. Default is 1.0",
)

parser.add_argument(
    "--bead_mass",
    type=float,
    default=5,
    help="Mass of each polymer bead. Default is 5",
)

parser.add_argument(
    "--box_dim",
    type=int,
    default=100,
    help="Edge length of the simulation box. Default is 30",
)

args = parser.parse_args()

# Parse arguments to variables
n_backbone: int = args.n_backbone
backbone_partners: list[str] = args.backbone_partners

if len(backbone_partners) > 0:
    backbone_partners = [int(i) for i in backbone_partners.split(" ")]
    n_sidechain: int = len(backbone_partners)
else:
    n_sidechain: int = 0
fene_k: float = args.fene_k
fene_r0: float = args.fene_r0
wca_sigma: float = args.wca_sigma
wca_epsilon: float = args.wca_epsilon
task_name: float = args.task_name
collision_period: int = args.collision_period
parallel_plate_speeds: float = read_floats_from_file(args.parallel_plate_speeds)

# n_sim_steps: The number of simulation steps to run. This is calculated by dividing the total number of steps by the number of steps per collision.
n_sim_steps: int = args.n_sim_steps

# kappa: The stiffness of the angle potential between backbone and sidechain particles.
kappa: float = args.angle_stiffness

# sidechains_can_bond: A boolean value indicating whether sidechain particles can form bonds with each other.
sidechains_can_bond: bool = args.sidechains_can_bond

sidechain_sidechain_fene_k: float = args.sidechain_sidechain_fene_k

# sidechain_equilibrium_bonding_angle: The equilibrium bonding angle between a sidechain particle and its associated backbone particle.
sidechain_equilibrium_bonding_angle: float = args.sidechain_equilibrium_bonding_angle

# sidechain_capture_radius: The capture radius for sidechain bonding in terms of wca_sigma. Do not multiply by wca_sigma again. This is used to determine if a sidechain particle is within the acceptable range of the equilibrium bonding angle.
sidechain_capture_radius: float = args.sidechain_capture_radius * wca_sigma
sidechain_sidechain_fene_r0: float = args.sidechain_sidechain_fene_r0 * wca_sigma
sidechain_sidechain_bond_break_length: float = (
    args.sidechain_sidechain_bond_break_length * wca_sigma
)
sidechain_sidechain_bond_break_angle_tolerance: float = (
    args.sidechain_sidechain_bond_break_angle_tolerance
)

# sidechain_bonding_angle_tolerance: The tolerance for the sidechain bonding angle. This is used to determine if a sidechain particle is within the acceptable range of the equilibrium bonding angle.
sidechain_bonding_angle_tolerance: float = args.sidechain_bonding_angle_tolerance

init_file_name = args.init_file_name
init_frame_number = args.init_frame_number

save_every = args.save_every
prod_save_folder = args.production_save_folder

init_conformations_folder = args.init_conformations_folder
reversible_bonds = args.reversible_bonds

# if reversible bonds is true, then sidechains can bond must also be true
sidechains_can_bond = sidechains_can_bond or reversible_bonds

# bath_temp: The temperature of the MPCD fluid bath.
bath_temp = args.bath_temp

# bead_mass: The mass of each polymer bead.
bead_mass = args.bead_mass

# box_dim: The edge length of the simulation box.
box_dim = args.box_dim

n_conformation_gen_steps = args.n_conformation_gen_steps

# must use the CPU to build the polymer, because adding and removing bonds must happen on the CPU
set_default_device("cpu")
