import torch
import gsd.hoomd
import hoomd.md as md
from .parse_args import (
    n_backbone,
    backbone_partners,
    n_sidechain,
    fene_k,
    fene_r0,
    wca_sigma,
    wca_epsilon,
    bead_mass,
    kappa,
    sidechain_sidechain_fene_k,
    sidechain_sidechain_fene_r0,
    box_dim,
)
import torch

# define the polymer
frame = gsd.hoomd.Frame()

# give enough space to avoid periodic self-interactions
# box must be a cube to prevent secondary flows from developing from periodic boundary conditions
box_length = box_dim
box_width = box_length

# total padding added to prevent MPCD particles from interacting through geometry
mpcd_padding = 4
box_height = box_length + mpcd_padding

# ParallelPlates puts the plates in the y direction, so the height is the distance between the plates
frame.configuration.box = [box_length, box_height, box_width, 0, 0, 0]

frame.particles.N = n_backbone + n_sidechain  # total number of particles
frame.particles.types = ["backbone", "sidechain"]  # particle types
frame.particles.typeid = [0] * n_backbone + [1] * n_sidechain  # particle type id

frame.bonds.N = n_backbone - 1 + n_sidechain  # number of bonds
frame.bonds.types = [
    "backbone-backbone",
    "backbone-sidechain",
    "sidechain-sidechain",
]  # bond type


# assign bond types. no sidechain-sidechain bonds yet
frame.bonds.typeid = [0] * (n_backbone - 1) + [1] * n_sidechain

# define bonds between backbone particles
frame.bonds.group = []
for backbone_particle in range(n_backbone - 1):
    frame.bonds.group.append([backbone_particle, backbone_particle + 1])

for particle in range(n_sidechain):
    # bond the sidechain with the selected backbone particle
    frame.bonds.group.append([n_backbone + particle, backbone_partners[particle]])


# bond potential is FENE. I will turn off the WCA part and define it separately.
fene = md.bond.FENEWCA()

# reasonable values taken from hoomd documentation
# include WCA parameters for the FENE potential because the pair potential exludes bonded particles
fene.params["backbone-backbone"] = dict(
    k=fene_k, r0=fene_r0, epsilon=wca_epsilon, sigma=wca_sigma, delta=0
)
# define the same parameters for the sidechain-backbone bond
fene.params["backbone-sidechain"] = dict(
    k=fene_k, r0=fene_r0, epsilon=wca_epsilon, sigma=wca_sigma, delta=0
)

fene.params["sidechain-sidechain"] = dict(
    k=sidechain_sidechain_fene_k,
    r0=sidechain_sidechain_fene_r0,
    epsilon=wca_epsilon,
    sigma=wca_sigma,
    delta=0,
)

# define the WCA potential between all non-bonded particles
lj_wca = md.pair.LJ(
    nlist=md.nlist.Cell(buffer=0.4), default_r_cut=2 ** (1 / 6) * wca_sigma
)

# define the WCA parameters between backbone particles
# setting r_cut automatically shifts the potential to zero at the cutoff
lj_wca.params[("backbone", "backbone")] = dict(epsilon=wca_epsilon, sigma=wca_sigma)

# same WCA parameters between backbone and sidechain particles
lj_wca.params[("backbone", "sidechain")] = dict(epsilon=wca_epsilon, sigma=wca_sigma)

# same WCA parameters between sidechain particles
lj_wca.params[("sidechain", "sidechain")] = dict(epsilon=wca_epsilon, sigma=wca_sigma)

# define the angles. ONLY THE BACKBONE PARTICLES GET ANGLE POTENTIALS
frame.angles.N = n_backbone - 2

# assign the angle types
frame.angles.types = [
    "backbone-backbone-backbone",
    "backbone-sidechain-sidechain",
]

# standard backbone-backbone-backbone
# backbone-sidechain-sidechain angles may form during the reactive phase but are not included in the initial configuration
frame.angles.typeid = [0] * (n_backbone - 2)

frame.angles.group = []
# define the angles between consecutive backbone particles
for particle in range(n_backbone - 2):
    frame.angles.group.append([particle, particle + 1, particle + 2])


# define the angle potential. This is directly from Bruce's code.
def bend_potential_and_torque(theta, kappa=kappa, constant=1.0):
    """Calculates bending potential energy and torque."""
    U = kappa * (constant + torch.cos(theta))  # bending potential
    T = kappa * torch.sin(theta)  # torque
    return U, T


cosinesq = md.angle.Table(width=5000)

# tabulate potential and torque for angles between 0 and pi radians (see HOOMD 2.9 docs)
x = torch.linspace(0, torch.pi, 5000)
angle_potentials, torque = bend_potential_and_torque(x)

# assume all angles are the same
cosinesq.params["backbone-backbone-backbone"] = dict(U=angle_potentials, tau=torque)
cosinesq.params["backbone-sidechain-sidechain"] = dict(U=angle_potentials, tau=torque)
# Bruce didn't use dihedrals, so I won't either

# %%
# build the polymer
# first, place the backbone particles on a diagonal line
# if there is an associated side chain, place it perpendicular to the backbone (it will equilibrate to the correct angle)

initial_spacing = 0.9 * wca_sigma

# empty array to store the positions
# for some reason, there is a fene bond error if the dtype is anything other than float32
frame.particles.position = torch.zeros((frame.particles.N, 3), dtype=torch.float32)

# place the backbone particles in the x-z plane
direction = torch.tensor([1, 0, 1], dtype=torch.float32) / torch.norm(
    torch.tensor([1, 0, 1], dtype=torch.float32)
)

# place the sidechain particle perpendicular to the backbone
sidechain_direction = torch.tensor([0, 1, 0], dtype=torch.float32) / torch.norm(
    torch.tensor([0, 1, 0], dtype=torch.float32)
)

# place the particles. First check if the first particle has a sidechain
if len(backbone_partners) > 0:
    if 0 in backbone_partners:
        frame.particles.position[n_backbone + backbone_partners.index(0)] = (
            frame.particles.position[0] + initial_spacing * sidechain_direction
        )

# then place the rest of the backbone particles and sidechain particles
for particle in range(1, n_backbone):
    frame.particles.position[particle] = (
        frame.particles.position[particle - 1] + initial_spacing * direction
    )
    if len(backbone_partners) > 0:
        if particle in backbone_partners:
            frame.particles.position[n_backbone + backbone_partners.index(particle)] = (
                frame.particles.position[particle]
                + initial_spacing * sidechain_direction
            )

# center the polymer
center_of_mass = torch.mean(frame.particles.position, dim=0)
frame.particles.position -= center_of_mass

frame.particles.mass = torch.ones(frame.particles.N, dtype=torch.float32) * bead_mass
