import itertools
import subprocess
import numpy as np
import os

"""Create precursor polymer chains and run the simulation."""

# reasoning: https://pubs.rsc.org/en/content/articlepdf/2015/sm/c4sm02475c says 400 beads is a good number
# but I'm using 100 to speed up the simulation by ensuring a box that can produce steady flow without secondary flows
# linker fractions
f = [0.1, 0.2, 0.4]

# blockiness
beta = [0.2, 0.8]

# chain stiffness: https://pubs.acs.org/doi/full/10.1021/acs.macromol.7b01876
kappa = [0, 5, 10]

# parallel plate speed
# shear rate here = 2 * plate_speed / 100
# so these correlate to shear rates of 0, 0.001, 0.005, and 0.01
plate_speed = ["constant=0.0", "constant=0.05", "constant=0.25", "constant=0.5"]

# 2 different backbone lengths, in terms of number of beads

# use 80 because it is divisible by 10 so we can do f=0.1, 0.2, 0.4
n_backbone = [50, 100, 150]

# create 5 chains for each variable combination
chains_per_combo = 5
combinations = list(itertools.product(f, beta, kappa, plate_speed, n_backbone)) * chains_per_combo

# Function to calculate the identity array I based on the selected indices
def calculate_I(beads_with_pendant_group: np.array, n_backbone: int):
    I = np.zeros(n_backbone - 1, dtype=int)
    for j in range(n_backbone - 1):
        if (j in beads_with_pendant_group and j + 1 in beads_with_pendant_group) or (
            j not in beads_with_pendant_group and j + 1 not in beads_with_pendant_group
        ):
            I[j] = 1
    return I


# Function to calculate the beta value based on I, n, and m
def calculate_beta(I: np.array, n_backbone: int, n_sidechain: int):
    f = n_sidechain / n_backbone
    b = 1 / n_backbone * (1 + np.sum(I))
    b_min = abs(2 * (f - 0.5))
    beta = (b - b_min) / (1 - b_min)
    return beta


def pick_backbone_beads_with_pendant_group(f, target_beta, N_backbone=100):
    # Step 1: Calculate b_min
    b_min = abs(2 * (f - 0.5))

    # Step 2: Calculate n_blocks
    b = target_beta * (1 - b_min) + b_min
    n_blocks = round(N_backbone * (1 - b) / 2)

    # Step 3: Calculate n_sidechain
    n_sidechain = int(f * N_backbone)

    # Step 4: Randomly distribute the sidechains across the blocks
    block_sizes = np.ones(n_blocks, dtype=int)  # Start with at least 1 per block
    remaining_sidechains = n_sidechain - n_blocks  # Remaining sidechains to distribute

    if remaining_sidechains > 0:
        # Randomly distribute the remaining sidechains across blocks
        for _ in range(remaining_sidechains):
            block_sizes[np.random.randint(n_blocks)] += 1

    # Step 5: Place the blocks with flexible shifting to allow sliding up to the edges
    total_block_length = sum(block_sizes) + (
        n_blocks - 1
    )  # Including gaps of at least 1
    free_space = N_backbone - total_block_length  # Remaining space to allow shifting

    selected_indices = []
    current_position = 0  # Start placing the blocks from position 0

    for i, size in enumerate(block_sizes):
        # Determine the space left after placing the blocks, and allow sliding within that range
        max_shift = (
            free_space if i == 0 else free_space // (n_blocks - i)
        )  # Larger freedom for sliding for later blocks
        shift = np.random.randint(0, max_shift + 1)

        start_idx = current_position + shift
        block_indices = list(range(start_idx, start_idx + size))
        selected_indices.extend(block_indices)

        # Move to the next block position considering the block size
        current_position = (
            start_idx + size + 1
        )  # Ensure a gap of at least 1 between blocks
        free_space -= (
            shift  # Reduce the remaining free space by the shift applied to this block
        )

    # Return the selected indices sorted
    return np.sort(selected_indices)


# name the combination and chain number
# for example, element at index 0 will be named combo0_chain0
# element at index 71 will be named combo71_chain0
# element at index 72 will be named combo0_chain1
n_unique_combos = len(combinations) // chains_per_combo
chain_names = [
    f"combo{i % n_unique_combos}_chain{i // n_unique_combos}"
    for i in range(len(combinations))
]

# read a single condition from the slurm script job array
idx = int(os.environ["SLURM_ARRAY_TASK_ID"])
input_combo = combinations[idx]
chain_name = chain_names[idx]

# create several examples for each combination
f_val, beta_val, kappa_val, plate_speed_val, n_backbone = input_combo
n_sidechain = int(f_val * n_backbone)

tolerance = 0.001  # Define tolerance

# Optimizing indices to achieve the target beta value
print("Optimizing indices to achieve the target beta value...")
beta_for_optimized_indices = None

# ensure that the beta value is within the tolerance. In the rare case that it is not, pick new indices
while (
    beta_for_optimized_indices is None
    or abs(beta_for_optimized_indices - beta_val) > tolerance
):
    backbone_partners = pick_backbone_beads_with_pendant_group(
        f_val, beta_val, n_backbone
    )
    backbone_partners = sorted(backbone_partners)  # Put indices in ascending order

    i_for_optimized_indices = calculate_I(backbone_partners, n_backbone)
    beta_for_optimized_indices = calculate_beta(
        i_for_optimized_indices, n_backbone, n_sidechain
    )

print(backbone_partners)
print(f"f: {f_val}, beta: {beta_for_optimized_indices}, target beta: {beta_val}")

# run the simulation
backbone_partners_str = " ".join(map(str, backbone_partners))

# command must be a list of strings
# define bath_temp as a variable
bath_temp = 1.0

# ensure no chain self-interactions by making the box edges large enough
box_dim = 100 if n_backbone <= 100 else n_backbone

# adjust plate speed to maintain constant shear rate when box size changes
if box_dim > 100:
    if plate_speed_val.startswith("constant="):
        speed = float(plate_speed_val.split("=")[1])
        plate_speed_val = f"constant={round(speed * box_dim/100, 6)}"
"end if"

command = [
    "python",
    "generate_seed_conformations.py",
    "--n_backbone",
    str(n_backbone),
    "--backbone_partners",
    backbone_partners_str,
    "--angle_stiffness",
    str(kappa_val),
    "--parallel_plate_speeds",
    plate_speed_val,
    "--bath_temp",
    str(bath_temp),
    "--task_name",
    f"{chain_name}_f={f_val}, beta={beta_val}, kappa={kappa_val}, plate_speed={plate_speed_val}_bath_temp={bath_temp}_box_dim={box_dim}_n_backbone={n_backbone}_n_sidechain={n_sidechain}",
    "--init_conformations_folder",
    "conformations_during_init_n=50_100_150",
    "--box_dim",
    str(box_dim)
]

# Execute the command
result = subprocess.run(command, capture_output=True, text=True)

# Optionally, print the output and any errors
print("Output:", result.stdout, flush=True)
print("Errors:", result.stderr, flush=True)
