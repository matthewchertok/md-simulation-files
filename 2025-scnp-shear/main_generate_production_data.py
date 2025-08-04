import subprocess
import os
import argparse

# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Change to the script directory
os.chdir(script_dir)

# Parse command line arguments
parser = argparse.ArgumentParser(description="Run production data generation")
parser.add_argument(
    "--reversible_bonding",
    type=bool,
    default=False,
    help="Enable reversible bonding (True/False)",
)
args = parser.parse_args()

reversible_bonding = args.reversible_bonding

init_conformations_folder = "conformations_during_init_n=50_100_150"

"""Load pre-defined configurations and run the reactive phase"""
# Load the list of files
# Change working directory up two levels
os.chdir("../..")
file_names = sorted(
    f"../data/{init_conformations_folder}/{name}"
    for name in os.listdir(f"../data/{init_conformations_folder}")
    if name.endswith(".gsd")
)
file_names = file_names * 10  # repeat each file 10 times because we want replicates
file_names.sort(
    key=lambda fname: os.path.basename(fname)
)  # sort so all replicates are together

chain_names = [
    name.split("init_conformations_")[1].split(".gsd")[0] for name in file_names
]

# Define the start and end indices for the job array. Currently set to run all files.
START_IDX = 0
END_IDX = 10800
file_names = file_names[START_IDX:END_IDX]
chain_names = chain_names[START_IDX:END_IDX]


# read a single condition from the slurm script job array
idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
file_name = file_names[idx]
chain_name = chain_names[idx]

# frame numbers 0-9, inclusive, are the replicates for each simulation
frame_number = idx % 10

frac_backbone_beads_with_sidechain = float(chain_name.split("_f=")[1].split(", ")[0])
kappa_val = float(chain_name.split("kappa=")[1].split(", ")[0])
plate_speed_val = chain_name.split("plate_speed=")[1].split("_")[0]
n_backbone = int(chain_name.split("n_backbone=")[1].split("_")[0])

bath_temp = chain_name.split("bath_temp=")[1].split("_")[0]

# run production using the first frame of the conformations GSD
n_sim_steps = int(2e7)

chain_name = f"{chain_name}_replicate{frame_number}"

production_save_folder = (
    "production_reversible_bonds_n=50_100_150"
    if reversible_bonding
    else "production_irreversible_bonds_n=50_100_150"
)
reversible_bonding_arg = (
    "--reversible_bonds" if reversible_bonding else "--sidechains_can_bond"
)

# ensure no chain self-interactions by making the box edges large enough
box_dim = 100 if n_backbone <= 100 else n_backbone

# adjust plate speed to maintain constant shear rate when box size changes
if box_dim > 100:
    if plate_speed_val.startswith("constant="):
        speed = float(plate_speed_val.split("=")[1])
        plate_speed_val = f"constant={round(speed * box_dim/100, 2)}"


command = [
    "python",
    "run_reactive_steps.py",
    "--n_backbone",
    str(n_backbone),
    "--angle_stiffness",
    str(kappa_val),
    "--parallel_plate_speeds",
    plate_speed_val,
    "--n_sim_steps",
    str(n_sim_steps),
    "--task_name",
    f"{chain_name}",
    "--init_file_name",
    file_name,
    "--init_frame_number",
    str(frame_number),
    reversible_bonding_arg,
    "--production_save_folder",
    production_save_folder,
    "--bath_temp",
    str(bath_temp),
    "--box_dim",
    str(box_dim),
]

# Execute the command
result = subprocess.run(command, capture_output=True, text=True)

# Optionally, print the output and any errors
print("Output:", result.stdout, flush=True)
print("Errors:", result.stderr, flush=True)
