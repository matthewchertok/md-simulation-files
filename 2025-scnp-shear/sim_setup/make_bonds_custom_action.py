import hoomd
import torch
from numba import cuda
import numpy as np
from .parse_args import (
    sidechain_capture_radius,
    sidechain_equilibrium_bonding_angle,
    sidechain_bonding_angle_tolerance,
)


class FormBonds(hoomd.custom.Action):
    """
    Custom action class for forming bonds.

    Args:
        backbone_partners (list): List of backbone partners.

    Attributes:
        backbone_partners (list): List of backbone partners.

    Methods:
        act(timestep): Perform the bond formation action.

    """

    def __init__(self, backbone_partners, has_cuda, can_break_bonds):
        self.backbone_partners = backbone_partners
        self.has_cuda = has_cuda
        self.can_break_bonds = can_break_bonds

    def act(self, _):
        if self.has_cuda:
            with self._state.gpu_local_snapshot as snapshot:
                new_snapshot = self._activate_bonds(
                    snapshot,
                    has_cuda=True,
                    backbone_partners=self.backbone_partners,
                    can_break_bonds=self.can_break_bonds,
                )
        else:
            with self._state.cpu_local_snapshot as snapshot:
                new_snapshot = self._activate_bonds(
                    snapshot,
                    has_cuda=False,
                    backbone_partners=self.backbone_partners,
                    can_break_bonds=self.can_break_bonds,
                )

        if new_snapshot is not None:
            self._state.set_snapshot(new_snapshot)

    def _activate_bonds(
        self, local_snapshot, has_cuda, backbone_partners, can_break_bonds
    ):
        """
        Update the bonds in the simulation based on sidechain-sidechain interactions.

        Args:
            simulation (Simulation): The hoomd simulation object.
            sim_step (int): The current simulation step. Required for printing debug messages.
            active_bonds (torch.Tensor): The bonds that are currently active in the simulation.
            backbone_partners (torch.Tensor): The indices of the backbone particles that are bonded to each sidechain particle.

        Returns:
            snapshot (Snapshot): The updated snapshot of the simulation state.
        """
        newly_active, newly_broken = self._compute_sidechain_bonds(
            local_snapshot,
            has_cuda=has_cuda,
            backbone_partners=backbone_partners,
            can_break_bonds=can_break_bonds,
        )
        snapshot = None
        # break bonds if reversible
        if newly_broken.size(0) > 0:
            snapshot = self._state.get_snapshot()

            # Remove bonds where either particle matches any particle in newly_broken
            # need to convert to cpu when working with snapshots.
            newly_broken_np = newly_broken.cpu().numpy()
            broken_particles = np.unique(newly_broken_np.flatten())

            prev_bonds_group = snapshot.bonds.group.copy()
            prev_bond_typeids = snapshot.bonds.typeid.copy()

            bonds_to_break_mask = (
                np.isin(prev_bonds_group[:, 0], newly_broken_np[:, 0])
                & np.isin(prev_bonds_group[:, 1], newly_broken_np[:, 1])
            ) | (
                np.isin(prev_bonds_group[:, 0], newly_broken_np[:, 1])
                & np.isin(prev_bonds_group[:, 1], newly_broken_np[:, 0])
            )
            bonds_to_keep_mask = ~bonds_to_break_mask

            snapshot.bonds.N = np.count_nonzero(bonds_to_keep_mask)
            snapshot.bonds.group[:] = prev_bonds_group[bonds_to_keep_mask]
            snapshot.bonds.typeid[:] = prev_bond_typeids[bonds_to_keep_mask]

            # Remove angles where any two elements match both elements in any bond in newly_broken
            prev_angles_group = snapshot.angles.group.copy()
            prev_angle_typeids = snapshot.angles.typeid.copy()

            angles_to_keep = ~(
                (
                    np.isin(prev_angles_group[:, 0], broken_particles)
                    & np.isin(prev_angles_group[:, 1], broken_particles)
                )
                | (
                    np.isin(prev_angles_group[:, 1], broken_particles)
                    & np.isin(prev_angles_group[:, 2], broken_particles)
                )
            )

            snapshot.angles.N = np.count_nonzero(angles_to_keep)
            snapshot.angles.group[:] = prev_angles_group[angles_to_keep]
            snapshot.angles.typeid[:] = prev_angle_typeids[angles_to_keep]

        # add bonds if there are any new bonds
        if newly_active.size(0) > 0:
            if snapshot is None:
                snapshot = self._state.get_snapshot()
            prev_bonds_group = snapshot.bonds.group.copy()
            prev_angles_group = snapshot.angles.group.copy()
            prev_bond_typeids = snapshot.bonds.typeid.copy()
            prev_angle_typeids = snapshot.angles.typeid.copy()
            backbone_partners_array = backbone_partners.cpu().numpy()

            # Add bonds and angles for newly active bonds
            sidechain1 = newly_active[:, 0].cpu().numpy()
            sidechain2 = newly_active[:, 1].cpu().numpy()

            new_bonds = np.column_stack((sidechain1, sidechain2))
            sidechain_sidechain_bond_typeid = 2
            new_bond_typeids = np.full(
                new_bonds.shape[0], sidechain_sidechain_bond_typeid
            )

            n_backbone = np.count_nonzero(snapshot.particles.typeid == 0)
            backbone1 = backbone_partners_array[sidechain1 - n_backbone]
            backbone2 = backbone_partners_array[sidechain2 - n_backbone]

            new_angles1 = np.column_stack((backbone1, sidechain1, sidechain2))
            new_angles2 = np.column_stack((backbone2, sidechain2, sidechain1))
            backbone_sidechain_sidechain_angle_typeid = 1
            new_angle_typeids = np.full(
                new_angles1.shape[0], backbone_sidechain_sidechain_angle_typeid
            )

            # Append new bonds and angles
            snapshot.bonds.N += new_bonds.shape[0]
            snapshot.bonds.group[:] = np.vstack([prev_bonds_group, new_bonds])
            snapshot.bonds.typeid[:] = np.concatenate(
                [prev_bond_typeids, new_bond_typeids]
            )

            snapshot.angles.N += 2 * new_angles1.shape[0]
            snapshot.angles.group[:] = np.vstack(
                [prev_angles_group, new_angles1, new_angles2]
            )
            snapshot.angles.typeid[:] = np.concatenate(
                [prev_angle_typeids, new_angle_typeids, new_angle_typeids]
            )
        """END ADD OR REMOVE BONDS"""
        return snapshot

    def _compute_sidechain_bonds(
        self, local_snapshot, has_cuda, backbone_partners, can_break_bonds
    ):
        """
        Update the bonds in the system based on the given snapshot and active bonds.

        Args:
            local_snapshot (Snapshot): The current snapshot of the system.
            has_cuda (bool): Whether or not the simulation is running on a GPU.
            backbone_partners (Tensor): Indices of the backbone particles bonded to each sidechain particle.
            can_break_bonds (bool): Flag to determine if bonds can be broken.

        Returns:
            tuple: Two tensors containing indices of newly active bonds and bonds to be broken.
        """
        # Get box dimensions
        box_size = torch.tensor(
            [
                local_snapshot.global_box.Lx,
                local_snapshot.global_box.Ly,
                local_snapshot.global_box.Lz,
            ],
            dtype=torch.float32,
        )

        # Get existing bonds
        existing_bonds = torch.tensor(local_snapshot.bonds.group, dtype=torch.int64)

        # Get all positions and particle IDs from the snapshot
        all_positions = torch.tensor(
            local_snapshot.particles.position, dtype=torch.float32
        )
        particle_ids = torch.tensor(local_snapshot.particles.typeid, dtype=torch.int64)

        # Reorder positions according to reverse tag, which puts them in the same order as the original snapshot
        if not has_cuda:
            all_rtags = torch.tensor(local_snapshot.particles.rtag, dtype=torch.int64)
        else:
            # For GPU, use copy to get writable array
            rtag_numba_array = cuda.as_cuda_array(local_snapshot.particles.rtag)
            all_rtags = torch.tensor(rtag_numba_array.copy_to_host(), dtype=torch.int64)

        all_positions = all_positions[all_rtags]
        particle_ids = particle_ids[all_rtags]

        # Identify sidechain particles
        sidechain_particle_indices = torch.where(particle_ids == 1)[0]
        sidechain_positions = all_positions[sidechain_particle_indices]
        n_sidechain = sidechain_particle_indices.size(0)

        # Compute distances between sidechain particles
        displacement_vectors = (
            sidechain_positions[:, None, :] - sidechain_positions[None, :, :]
        )
        displacement_vectors -= box_size * torch.round(displacement_vectors / box_size)
        distances = torch.norm(displacement_vectors, dim=2)

        # Get all unique pairs of sidechain particles (upper triangle)
        i, j = torch.triu_indices(n_sidechain, n_sidechain, offset=1)
        sidechain_pairs = torch.stack([i, j], dim=1)
        pair_distances = distances[i, j]

        # Filter pairs within the capture radius
        valid_length_mask = pair_distances <= sidechain_capture_radius

        # Map sidechain indices to global indices in all_positions
        sidechain_pairs_global = torch.stack(
            [
                sidechain_particle_indices[sidechain_pairs[:, 0]],
                sidechain_particle_indices[sidechain_pairs[:, 1]],
            ],
            dim=1,
        )

        # Filter pairs with valid distances
        valid_pairs_global = sidechain_pairs_global[valid_length_mask]
        valid_pair_indices = sidechain_pairs[valid_length_mask]

        # Get backbone partners for each sidechain particle
        backbone_indices = backbone_partners  # Indices into all_positions

        # Compute bond vectors for angle calculations
        backbone1_indices = backbone_indices[valid_pair_indices[:, 0]]
        backbone2_indices = backbone_indices[valid_pair_indices[:, 1]]

        bond_backbone1_sidechain1 = (
            sidechain_positions[valid_pair_indices[:, 0]]
            - all_positions[backbone1_indices]
        )
        bond_sidechain1_sidechain2 = (
            sidechain_positions[valid_pair_indices[:, 1]]
            - sidechain_positions[valid_pair_indices[:, 0]]
        )
        bond_sidechain2_backbone2 = (
            all_positions[backbone2_indices]
            - sidechain_positions[valid_pair_indices[:, 1]]
        )

        # Apply minimum image convention
        bond_backbone1_sidechain1 -= box_size * torch.round(
            bond_backbone1_sidechain1 / box_size
        )
        bond_sidechain1_sidechain2 -= box_size * torch.round(
            bond_sidechain1_sidechain2 / box_size
        )
        bond_sidechain2_backbone2 -= box_size * torch.round(
            bond_sidechain2_backbone2 / box_size
        )

        # Compute angles
        angle1 = torch.acos(
            torch.sum(-bond_backbone1_sidechain1 * bond_sidechain1_sidechain2, dim=1)
            / (
                torch.norm(bond_backbone1_sidechain1, dim=1)
                * torch.norm(bond_sidechain1_sidechain2, dim=1)
            )
        )
        angle2 = torch.acos(
            torch.sum(-bond_sidechain1_sidechain2 * bond_sidechain2_backbone2, dim=1)
            / (
                torch.norm(bond_sidechain1_sidechain2, dim=1)
                * torch.norm(bond_sidechain2_backbone2, dim=1)
            )
        )

        # Check validity of angles
        angle_diff1 = torch.abs(angle1 - sidechain_equilibrium_bonding_angle)
        angle_diff2 = torch.abs(angle2 - sidechain_equilibrium_bonding_angle)
        valid_angle_mask = (angle_diff1 <= sidechain_bonding_angle_tolerance) & (
            angle_diff2 <= sidechain_bonding_angle_tolerance
        )

        # Apply valid angle mask
        valid_pairs_global = valid_pairs_global[valid_angle_mask]
        valid_pair_indices = valid_pair_indices[valid_angle_mask]

        # Enforce mutual closest partner using vectorized operations
        # Compute distances between all valid pairs
        valid_distances = pair_distances[valid_length_mask][valid_angle_mask]

        # Get sidechain particle indices for valid pairs
        indices_i = valid_pair_indices[:, 0]
        indices_j = valid_pair_indices[:, 1]

        # Create a matrix to hold distances (initialized with infinity)
        n_sidechain = sidechain_particle_indices.size(0)
        distance_matrix = torch.full((n_sidechain, n_sidechain), float("inf"))

        # Fill in the distances for valid pairs
        distance_matrix[indices_i, indices_j] = valid_distances
        distance_matrix[indices_j, indices_i] = valid_distances  # Symmetric matrix

        # Find the closest partner for each particle
        _, closest_indices = distance_matrix.min(dim=1)

        # Check for mutual closest partners
        indices = torch.arange(n_sidechain)
        mutual_mask = (closest_indices[closest_indices] == indices) & (
            indices != closest_indices
        )

        # Extract mutually closest pairs
        mutual_pairs = torch.stack(
            [indices[mutual_mask], closest_indices[mutual_mask]], dim=1
        )
        mutual_pairs = mutual_pairs[
            mutual_pairs[:, 0] < mutual_pairs[:, 1]
        ]  # Avoid duplicates

        # Map local sidechain indices to global indices
        new_bonds = torch.stack(
            [
                sidechain_particle_indices[mutual_pairs[:, 0]],
                sidechain_particle_indices[mutual_pairs[:, 1]],
            ],
            dim=1,
        )

        # Remove bonds involving particles already bonded
        existing_sidechain_bonds = existing_bonds[
            torch.isin(existing_bonds[:, 0], sidechain_particle_indices)
            & torch.isin(existing_bonds[:, 1], sidechain_particle_indices)
        ]
        existing_bonded_particles = existing_sidechain_bonds.flatten().unique()
        mask = ~torch.isin(new_bonds[:, 0], existing_bonded_particles) & ~torch.isin(
            new_bonds[:, 1], existing_bonded_particles
        )
        new_bonds = new_bonds[mask]

        # Compute bonds to break if allowed
        if can_break_bonds:
            sidechain_bonds = existing_sidechain_bonds

            if sidechain_bonds.size(0) == 0:
                newly_broken_bonds = torch.empty((0, 2), dtype=torch.int64)
            else:
                # Compute bond vectors and lengths
                bond_vecs = (
                    all_positions[sidechain_bonds[:, 0]]
                    - all_positions[sidechain_bonds[:, 1]]
                )
                bond_vecs -= box_size * torch.round(bond_vecs / box_size)
                bond_lengths = torch.norm(bond_vecs, dim=1)

                # Get backbone partners
                sc1_indices = torch.searchsorted(
                    sidechain_particle_indices, sidechain_bonds[:, 0]
                )
                sc2_indices = torch.searchsorted(
                    sidechain_particle_indices, sidechain_bonds[:, 1]
                )
                bb1_indices = backbone_indices[sc1_indices]
                bb2_indices = backbone_indices[sc2_indices]

                # Compute angle vectors
                vec_bb1_sc1 = (
                    all_positions[sidechain_bonds[:, 0]] - all_positions[bb1_indices]
                )
                vec_sc1_sc2 = (
                    all_positions[sidechain_bonds[:, 1]]
                    - all_positions[sidechain_bonds[:, 0]]
                )
                vec_sc2_bb2 = (
                    all_positions[bb2_indices] - all_positions[sidechain_bonds[:, 1]]
                )

                vec_bb1_sc1 -= box_size * torch.round(vec_bb1_sc1 / box_size)
                vec_sc1_sc2 -= box_size * torch.round(vec_sc1_sc2 / box_size)
                vec_sc2_bb2 -= box_size * torch.round(vec_sc2_bb2 / box_size)

                # Compute angles
                angle1 = torch.acos(
                    torch.sum(vec_bb1_sc1 * -vec_sc1_sc2, dim=1)
                    / (torch.norm(vec_bb1_sc1, dim=1) * torch.norm(vec_sc1_sc2, dim=1))
                )
                angle2 = torch.acos(
                    torch.sum(vec_sc1_sc2 * -vec_sc2_bb2, dim=1)
                    / (torch.norm(vec_sc1_sc2, dim=1) * torch.norm(vec_sc2_bb2, dim=1))
                )

                # Compute break probabilities
                p_break_r = 1.0 / (1.0 + torch.exp(-30.0 * (bond_lengths - 1.4)))
                angle_weight_factor = 1e-3
                p_break_theta1 = (
                    1.0
                    / (1.0 + torch.exp(-3.0 * (angle1 - (7.0 * np.pi / 12.0))))
                    * angle_weight_factor
                )
                p_break_theta2 = (
                    1.0
                    / (1.0 + torch.exp(-3.0 * (angle2 - (7.0 * np.pi / 12.0))))
                    * angle_weight_factor
                )

                # Decide which bonds to break
                random_values = torch.rand(p_break_r.size(0))
                break_mask = (
                    (random_values < p_break_r)
                    | (random_values < p_break_theta1)
                    | (random_values < p_break_theta2)
                )

                newly_broken_bonds = sidechain_bonds[break_mask]
        else:
            newly_broken_bonds = torch.empty((0, 2), dtype=torch.int64)

        return new_bonds, newly_broken_bonds

    """END ACTIVATE BONDS"""


# end FormBonds
