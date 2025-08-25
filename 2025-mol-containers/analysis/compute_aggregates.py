import argparse
from joblib import Parallel, delayed
import numpy as np
from scipy.spatial.distance import cdist

import functools
print = functools.partial(print, flush=True)

def gen_mol_dict(data_file):
    '''
        Returns a dictionary of molecules from a sys.data file, where the keys
        are the molecule IDs and the values are the atom IDs.
    '''

    with open(data_file, 'r') as handle:
        lines = handle.readlines()
    atom_line = 0
    bond_line = 0
    for idx, line in enumerate(lines):
        if 'Atoms' in line.strip():
            atom_line = idx
        if 'Bonds' in line.strip():
            bond_line = idx
        if 'Velocities' in line.strip():
            bond_line = idx
        if atom_line != 0 and bond_line != 0:
            break
    coord_start = atom_line + 2
    coord_end = bond_line - 1
    coord_lines = lines[coord_start:coord_end]

    mol_dict = {}
    for line in coord_lines:
        vals = line.strip().split()
        if int(vals[1]) not in mol_dict:
            mol_dict[int(vals[1])] = []
        mol_dict[int(vals[1])].append(int(vals[0]))

    return mol_dict

def get_frame_length(traj_file):

    step_counter = 0
    frame_length = 0
    with open(traj_file, 'r') as f:
        while step_counter < 2:
            if 'TIMESTEP' in f.readline().strip():
                step_counter += 1
            frame_length += 1
    frame_length -= 1
    
    return frame_length

def compute_aggregate(frame_id, lines, mol_dict, cutoff=3.0):
    '''
        Returns the IDs of molecules that are contained in the same aggregate
        as the start_molecule, per the specified cutoff.
    '''

    atom_to_mol_dict = {}
    included_atoms_unsorted = []
    for mol, atoms in mol_dict.items():
        for atom in atoms:
            atom_to_mol_dict[atom] = mol
            included_atoms_unsorted.append(atom)

    coords = []
    included_atoms = []
    for line in lines[9:]:
        vals = line.strip().split()
        if int(vals[0]) in included_atoms_unsorted:
            coord = [float(vals[2]), float(vals[3]), float(vals[4])]
            coords.append(coord)
            included_atoms.append(int(vals[0]))
    coords = np.array(coords)

    container_molecules = [k for k,v in mol_dict.items() if len(v) > 75]

    largest_aggregate = 0
    for container_molecule in container_molecules:

        molecules = [container_molecule]
        added_molecules = True
        while added_molecules:

            atoms = []
            for mol in molecules:
                for atom in mol_dict[mol]:
                    atoms.append(included_atoms.index(atom))
            atoms = np.array(atoms)
            atoms = atoms.astype(int)
            agg_coords = coords[atoms]

            dist = cdist(coords, agg_coords)
            within_cut = np.where(np.min(dist, axis=1) < cutoff, 1, 0).reshape(-1)

            added_molecules = False
            for atom in range(within_cut.shape[0]):
                if atom_to_mol_dict[included_atoms[atom]] not in molecules and within_cut[atom] == 1.0:
                    molecules.append(atom_to_mol_dict[included_atoms[atom]])
                    added_molecules = True

        container_counter = 1
        for mol in container_molecules.copy():
            if mol != container_molecule and mol in molecules:
                container_molecules.remove(mol)
                container_counter += 1

        surf_number = len(molecules) - container_counter
        if surf_number > largest_aggregate:
            largest_aggregate = surf_number

    print(f'Completed frame: {frame_id}.')

    return frame_id, largest_aggregate

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sys_path')
    parser.add_argument('--traj_file')
    parser.add_argument('--out', default='test.npy')
    args = parser.parse_args()
    sys_path = args.sys_path
    traj_file = args.traj_file

    mol_dict = gen_mol_dict(data_file=sys_path)
    mol_dict = {k:v for k,v in mol_dict.items() if len(v) > 10}

    frame_length = get_frame_length(traj_file=traj_file)
    traj_lines = open(traj_file, 'r').readlines()
    num_frames = int(len(traj_lines) / frame_length)
    print(f'Computing for {num_frames} number of frames...')

    results = Parallel(n_jobs=-1)(
        delayed(compute_aggregate)(
            lines=traj_lines[idx * frame_length: (idx+1) * frame_length], 
            mol_dict=mol_dict, 
            cutoff=3.0,
            frame_id=idx
        ) for idx in range(num_frames)
    )

    results = np.array(results)
    results = results[np.argsort(results[:,0])]
    np.save(f'./data/{args.out}', results)

    print('Completed.')
