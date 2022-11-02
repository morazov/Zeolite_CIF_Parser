import numpy as np
import pandas as pd
from fractions import Fraction
import os

# This script converts all of the CIF files found on:
# 'https://america.iza-structure.org/IZA-SC/ftc_table.php'
# from a compact CIF format (primitive cell + symmetry operations)
# to a more useful XYZ format that can be used directly to visualize the full
# unit cell, with all symmetry operations applied to the primitive unit cell,
# and duplicate generated locations from multiple symmetry operations removed.

# Please cite the crystallographic data source if you wish to use the data:
# Ch. Baerlocher and L.B. McCusker
# Database of Zeolite Structures: http://www.iza-structure.org/databases/


def symmetry_generate(original_positions, symmetry_matrix, offset_vector):
    """Generates equivalent positions, given a symmetry matrix and offset

    Args:
        original_positions (numpy.array): fractional coordinate vectors of
        points in a primitive unit cell; shape is 3 by N array, where N is
        number of points
        symmetry_matrix (numpy.array): transformation matrix that generates
        equivalent position; shape is 3 by 3 array
        offset_vector (numpy.array): offset component of the tranformation;
        shape is 3 by 1 array

    Returns:
        (numpy.array): symmetry equivalent positions in unit cell in fractional
        coordinates; shape = 3 by n array
    """
    new_positions = np.dot(symmetry_matrix, original_positions) + offset_vector
    # Translate generated positions that are outside of the unit cell back into
    # unit cell
    new_positions[new_positions < 0] += 1
    new_positions[new_positions > 1] -= 1
    return new_positions


def parse_sym_str(symmetry_string, separator=','):
    """Converts a symmetry operation string into matrix/vector format

    The CIF files found in the IZA database represent symmetry operations with
    symmetry vector strings, but the format is not consistent for all files.
    This code parses a string of type "1/2-y,+x,-1/4+z" or "-y+1/2,+x,+z-1/4"
    into a corresponding coordinate rotation matrix and constant offset vector.
    This covers all formats used in the IZA database, as of October, 2022.

    Args:
        symmetry_string (str): The symmetry string to be parsed
        separator (str, optional): Separator character for different
        dimensions. Defaults to ','.

    Returns:
        (list): Contains rotation matrix (0th element) and offset vector
        (1st element), both being numpy arrays.
    """

    xyz = symmetry_string.split(separator)  # (x, y, z) transforms
    M = np.zeros((3, 3))  # rotation matrix
    offset = np.zeros((3, 1))  # offset vector

    for row in range(3):
        x = xyz[row].find('x')  # index of 'x' in transform
        y = xyz[row].find('y')  # index of 'y' in transform
        z = xyz[row].find('z')  # index of 'z' in transform

        if x > 0:
            M[row, :] += int(xyz[row][x-1]+'1')*np.array([1, 0, 0])
        elif x == 0:
            M[row, :] += np.array([1, 0, 0])

        if y > 0:
            M[row, :] += int(xyz[row][y-1]+'1')*np.array([0, 1, 0])
        elif y == 0:
            M[row, :] += np.array([0, 1, 0])

        if z > 0:
            M[row, :] += int(xyz[row][z-1]+'1')*np.array([0, 0, 1])
        elif z == 0:
            M[row, :] += np.array([0, 0, 1])

        xyz_loc = np.array([x, y, z])
        xyz_loc = xyz_loc[xyz_loc > -1]
        first_var_loc = min(xyz_loc)
        last_var_loc = max(xyz_loc)
        if first_var_loc > 1:
            offset[row] += float(Fraction(xyz[row][:first_var_loc-1]))
        if last_var_loc < len(xyz[row])-1:
            offset[row] += float(Fraction(xyz[row][last_var_loc+1:]))

    return [M, offset]


def fractional_to_cartesian_matrix(a, b, c, alpha, beta, gamma):
    """Fractional to Cartesian coordinates transformation matrix.

    Generates a matrix that takes fractional coordinates from a unit cell to
    Cartesian coordinates, given unit cell dimensions. E.g.: X=A*x, where
    X is the Cartesian coordinate vector for the point x in fractional
    coordinates, and A is the transformation matrix.

    Args:
        a (float): cell length in a dimension
        c (float): cell length in b dimension
        b (float): cell length in c dimension
        alpha (float): angle between a and c axes, in degrees
        beta (float):  angle between b and c axes, in degrees
        gamma (float): angle between a and b axes, in degrees

    Returns:
        (numpy.array): 3x3 transformation matrix
    """
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    gamma = np.deg2rad(gamma)
    volume = a*b*c*np.sqrt(1-np.cos(alpha)**2 - np.cos(beta)**2
                           - np.cos(gamma)**2
                           + 2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    return np.array([[a, b*np.cos(gamma), c*np.cos(beta)],
                     [0, b*np.sin(gamma), c*(np.cos(alpha)-np.cos(beta) *
                                             np.cos(gamma))/np.sin(gamma)],
                     [0, 0, volume/(a*b*np.sin(gamma))]])


# Process all files generated by IZA_CIFs_scrape.py script to convert from
# CIF format to XYZ format (site and atom types preserved)
path = "./CIFs"
dir_list = os.listdir(path)
for file in dir_list:
    # Read file and store individual lines for parsing
    with open(f"./CIFs/{file}", 'r') as f:
        lines = f.readlines()

    # Initialize counters to help parse file
    count = 0  # Counter for line number entering into a "loop_" statement
    subcount = 0  # Counter for line number in a given "loop_" statement
    sym_count = 0  # Counter for number of symmetry operations

    # Initialize dataframe for symmetry operator strings, matrix, and offset
    sym_block = pd.DataFrame(columns=['sym_str', 'M', 'offset'])

    # Initialize dataframe to store fractional coordinates for primitive cell
    data_block = pd.DataFrame(columns=['Site', 'Atom', 'x', 'y', 'z'])

    # Find the unit cell dimensions (a, b, c, alpha, beta, gamma)
    for line in lines:
        if line.find('_cell_length_a', 0, 25) > -1:
            a = float(line.split()[-1].split('(')[0])
            b = float(lines[subcount+1].split()[-1].split('(')[0])
            c = float(lines[subcount+2].split()[-1].split('(')[0])
            alpha = float(lines[subcount+3].split()[-1].split('(')[0])
            beta = float(lines[subcount+4].split()[-1].split('(')[0])
            gamma = float(lines[subcount+5].split()[-1].split('(')[0])
            count = subcount + 6
            subcount = 0
            break
        subcount += 1

    # Find the start of the symmetry operation strings
    for line in lines[count:]:
        subcount += 1
        if line.find('_symmetry_equiv_pos_as_xyz', 0, 30) > -1:
            count += subcount
            subcount = 0
            break

    # Parse all symmetry operations,storing all results in sym_block
    for line in lines[count:]:
        subcount += 1
        # Distinguish between symmetry string formats
        if line[0] == "'":
            sym_str = line[1:-2]
            M, offset = parse_sym_str(sym_str)
            sym_block.loc[sym_count] = [sym_str, M, offset]
            sym_count += 1
        elif len(line.replace(' ', '')) > 1:
            sym_str = line[:-1].lower()
            M, offset = parse_sym_str(sym_str)
            sym_block.loc[sym_count] = [sym_str, M, offset]
            sym_count += 1
        else:
            count += subcount
            subcount = 0
            break

    # Find the start of the atom site fractional coordinates
    for line in lines[count:]:
        subcount += 1
        if line.find('_atom_site_fract_z', 0, 20) > -1:
            count += subcount
            subcount = 0
            break

    # Store the fractional coordinates of all unique sites in primitive cell
    for line in lines[count:]:
        if len(line) > 1:
            data_block.loc[subcount] = line.split()
        subcount += 1

    # Convert strings stored in pandas dataframe to numpy float array
    data_numpy = data_block[['x', 'y', 'z']].astype(float).to_numpy()

    # Apply symmetry tranform for each equivalent position
    all_atoms = np.array([np.transpose(symmetry_generate(np.transpose(
        data_numpy), M, offset)) for M, offset in zip(sym_block['M'],
                                                      sym_block['offset'])])

    # Reshape result into 3xN array
    all_atoms = np.transpose(all_atoms.reshape(-1, 3))

    # Extend framework up to 2x2x2 unit cells, as long as size increase is
    # under threshold. Note: some of the atoms will be removed later because
    # multiple symmetry operations can point to the same position in the
    # unit cell, so the atom threshold may be fine-tuned for specific
    # frameworks of interest, to show more unit cells. Here a 2000 atom
    # threshold is used as a starting point.
    atom_threshold = 2000
    number_atoms_in_cell = np.unique(all_atoms, axis=1).shape[1]

    if number_atoms_in_cell*2 < atom_threshold:
        all_atoms = np.concatenate(
            (all_atoms, (all_atoms+np.array([1, 0, 0]).reshape(3, 1))), axis=1)
    if number_atoms_in_cell*4 < atom_threshold:
        all_atoms = np.concatenate(
            (all_atoms, (all_atoms+np.array([0, 1, 0]).reshape(3, 1))), axis=1)
    if number_atoms_in_cell*8 < atom_threshold:
        all_atoms = np.concatenate(
            (all_atoms, (all_atoms+np.array([0, 0, 1]).reshape(3, 1))), axis=1)

    # Calculate the transform matrix to take fractional coordinate to Cartesian
    cartesian_transform = fractional_to_cartesian_matrix(
        a, b, c, alpha, beta, gamma)

    # Apply the Cartesian transform on the unit cell
    unit_cell = pd.DataFrame(
        np.transpose(np.dot(cartesian_transform, np.concatenate(
            [np.identity(3), np.ones((3, 1))], axis=1))),
        columns=['x', 'y', 'z'])

    # Apply the Cartesian transform on atom positions
    cartesian_df = pd.DataFrame(np.transpose(
        np.dot(cartesian_transform, all_atoms)), columns=['x', 'y', 'z'])

    # Add Site and Atom labels to each position
    cartesian_df[['Site', 'Atom']] = pd.concat(
        [data_block.loc[:, ['Site', 'Atom']]]
        * (len(cartesian_df) // len(data_block)), ignore_index=True)

    # Find duplicate atoms (same cartesian position generated by diffent
    # symmetry operations). Note: pandas dropduplicates() will not spot
    # all instances of duplicates because of roundoff errors. Need to
    # check distances. A good cutoff is the H-H bond distance, as it is
    # the shortest possible bond, and anything shorter is not feasible.
    threshold = 0.74**2
    # Get site labels of all unique sites in primitive cell
    sites = data_block['Site']
    for site in sites:
        # Find all instances of a given site and get their coordinates
        site_indices = cartesian_df['Site'] == site
        site_xyzs = cartesian_df.loc[site_indices, ['x', 'y', 'z']].to_numpy()
        site_rows = np.arange(site_xyzs.shape[0])  # enumerate the instances
        disgard = []  # accumulate duplicates
        for row in site_rows:
            # Only check instanes if they are not duplicates
            if row not in disgard:
                # Get the distances between an instance and all other instances
                distance_sq = np.sum((site_xyzs-site_xyzs[row, :])**2, axis=1)
                # Check for duplicates, but do not disgard self
                too_close = distance_sq < threshold
                too_close[row] = False
                too_close, = np.where(too_close)
                disgard += too_close.tolist()
        cartesian_df.drop(
            site_indices[site_indices].iloc[disgard].index, inplace=True)

    # Store unit cell info and the number of cells shown
    unit_cell.to_csv(f"./UnitCells/{file[:-4]}.csv", index=False)
    with open(f"./UnitCells/{file[:-4]}.csv", 'a') as cell:
        cell.write(f"Number of atoms in unit cell: {number_atoms_in_cell}")
    # Store XYZ data
    cartesian_df.to_csv(f"./XYZs/{file[:-4]}.csv", index=False)
