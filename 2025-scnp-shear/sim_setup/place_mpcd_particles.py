from .place_md_particles import frame, mpcd_padding
from hoomd import Snapshot
from hoomd.communicator import Communicator
from torch import tensor, rand, randn, sqrt, set_default_device
from .parse_args import bath_temp


def create_snapshot_gsd_mpcd(gsd_frame):
    """
    Create a snapshot with MPCD particles from a GSD frame.

    Args:
        gsd_frame (GSDFrame): The GSD frame to create the snapshot from.

    Returns:
        snapshot (Snapshot): The created snapshot with MPCD particles.
        filled_height (int): The height of the region of the simulation box filled with MPCD particles.
        density (float): The density of MPCD particles in the simulation box.
        kT (float): The temperature of the MPCD particles.
    """
    # set device to CPU to initialize the simulation, since this must be done on the CPU
    set_default_device("cpu")
    snapshot = Snapshot.from_gsd_frame(gsd_frame, Communicator())

    # this is how hoomd defines the directions
    length = snapshot.configuration.box[0]
    height = snapshot.configuration.box[1]
    width = snapshot.configuration.box[2]

    # according to the hoomd docs:
    # The simulation box needs padded in the y direction to account for
    # the collision cells wrapping through the y periodic boundary.
    # A padding of 4 cells (2 each direction) should be sufficient.
    # the fraction of the box that is filled with MPCD particles
    filled_height = height - mpcd_padding
    density = 5.0  # number of MPCD particles per unit volume

    snapshot.mpcd.types = ["solvent"]
    snapshot.mpcd.N = int(density * length * filled_height * width)

    # randomly place particles within the MPCD region
    snapshot.mpcd.position[:] = tensor(
        [0.5 * length, 0.5 * filled_height, 0.5 * width]
    ) * (
        rand(snapshot.mpcd.N, 3) * 2 - 1
    )  # x2 - 1 to get the range -1 to 1

    # initialize MPCD particle velocities according to the maxwell-boltzmann distribution. Here, we set k_b*T=1
    kT = bath_temp
    scale = sqrt(tensor(kT))
    shape = (snapshot.mpcd.N, 3)
    velocity = randn(*shape) * scale
    velocity -= velocity.mean(dim=0)  # Shift to guarantee the mean velocity is 0
    snapshot.mpcd.velocity[:] = velocity
    return snapshot, filled_height, density, kT


"end def"

snapshot, filled_height, density, kT = create_snapshot_gsd_mpcd(frame)
