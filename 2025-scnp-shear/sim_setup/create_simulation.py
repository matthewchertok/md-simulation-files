import hoomd
from .place_mpcd_particles import snapshot, filled_height, density, kT
from .place_md_particles import fene, lj_wca, cosinesq
from .parse_args import (
    collision_period,
    parallel_plate_speeds,
    init_file_name,
)

# create a simulation state if no initial file is provided
should_create_sim_state = len(init_file_name) == 0

# create a simulation. use GPU if available, otherwise use CPU
device = hoomd.device.auto_select()
simulation = hoomd.Simulation(device=device, seed=1)

# if an initial file is provided, we will load the state in a later step
if should_create_sim_state:
    simulation.create_state_from_snapshot(snapshot)

# configure the MPCD integrator, which handles both MPCD and MD particles
forces = [fene, lj_wca, cosinesq]

dt = 0.005  # time step for the simulation, per the hoomd docs (https://hoomd-blue.readthedocs.io/en/latest/tutorial/09-Multiparticle-Collision-Dynamics/02-Diffusion-in-Solution.html)
integrator = hoomd.mpcd.Integrator(
    dt=dt, forces=forces
)  # must specify the forces acting on the MD particles. dt must be small with FENE potential, otherwise the bonds will blow up.

# set up stochastic rotation dynamics (SRD) collisions for the MPCD particles
integrator.collision_method = hoomd.mpcd.collide.StochasticRotationDynamics(
    period=collision_period, angle=130, kT=kT, embedded_particles=hoomd.filter.All()
)

# sort MPCD particles in memory for computational efficiency
integrator.mpcd_particle_sorter = hoomd.mpcd.tune.ParticleSorter(
    trigger=integrator.collision_method.period * 20
)

# define the streaming method based on the parallel plate geometry (this will set the velocity profile according to Navier-Stokes)
plates = hoomd.mpcd.geometry.ParallelPlates(
    separation=filled_height, speed=parallel_plate_speeds[0], no_slip=True
)  # speed creates shear

# Use bounce back integration because there are boundary conditions on the fluid (parallel plates). If there weren't, I would've used Bulk.
integrator.streaming_method = hoomd.mpcd.stream.BounceBack(
    period=integrator.collision_method.period, geometry=plates
)
filler = hoomd.mpcd.fill.GeometryFiller(
    type="solvent", density=density, kT=kT, geometry=plates
)
integrator.virtual_particle_fillers.append(filler)

# use the NVE ensemble on MD particles because MPCD particles act like a temperature bath
# replacement for hoomd.md.methods.ConstantVolume(); this also confines MD particles between the plates
nve = hoomd.mpcd.methods.BounceBack(filter=hoomd.filter.All(), geometry=plates)
integrator.methods.append(nve)

simulation.operations.integrator = integrator