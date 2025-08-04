import hoomd
from torch.cuda import is_available
from torch import set_default_device
from sim_setup.parse_args import (
    task_name,
    save_every,
    n_conformation_gen_steps,
    init_conformations_folder,
)
from sim_setup.create_simulation import simulation
import time
import os

start_time = time.time()
results_filename = f"init_conformations_{task_name}"

# use a GPU to run the simulation, if available
has_cuda = is_available()
if has_cuda:
    set_default_device("cuda")
print(f"Using CUDA: {has_cuda}")

# bring the simulation to steady state
print("Equilibrating...")
simulation.run(int(1e6))
print("Equilibration complete")

output_dir = init_conformations_folder
os.makedirs(output_dir, exist_ok=True)

# check if file exists and abort if it does
if os.path.exists(f"{output_dir}/{results_filename}.gsd"):
    print(
        f"File {output_dir}/{results_filename}.gsd already exists. Aborting to avoid overwriting."
    )
    exit(1)
gsd_writer_init = hoomd.write.GSD(
    trigger=save_every,
    filename=f"{output_dir}/{results_filename}.gsd",
    dynamic=["property", "topology", "momentum", "attribute"],
)
simulation.operations.writers.append(gsd_writer_init)

# if needed, generate seed conformations
print("Generating seed conformations...")
simulation.run(int(n_conformation_gen_steps))
print("Seed conformations generated")
gsd_writer_init.flush()  # write the data to disk

end_time = time.time()

# end_time - start_time is the total runtime in seconds
total_seconds = end_time - start_time

# Convert to hours, minutes, and seconds
hours, remainder = divmod(total_seconds, 3600)
minutes, seconds = divmod(remainder, 60)

print(f"Runtime: {int(hours):02}:{int(minutes):02}:{int(seconds):02}")
