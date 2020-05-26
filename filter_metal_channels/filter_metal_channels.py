# -*- coding: utf-8 -*-

from __future__ import absolute_import

import logging
from operator import itemgetter

import numpy as np
from matminer.featurizers.base import MultipleFeaturizer
from matminer.featurizers.composition import (ElementProperty, Stoichiometry, ValenceOrbital)
from matminer.featurizers.structure import (ChemicalOrdering, DensityFeatures, MaximumPackingEfficiency,
                                            SiteStatsFingerprint, StructuralHeterogeneity, StructureComposition)
from mendeleev import element
from numba import jit
from pymatgen import Structure
from pymatgen.analysis import dimensionality
from six.moves import zip

from .mathutils import cross, dot, norm

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


def get_metallic_substructure(  # pylint:disable=dangerous-default-value, too-many-arguments
    s: Structure,
    return_non_metal_structure: bool = True,
    white_list: list = ['O', 'H'],
    allow_expansion: bool = False,
    primitive: bool = True,
    other_metals: list = ['Al', 'In', 'Ga'],
) -> Structure:
    """
    By default, it will create a Niggli reduced structure which is better for be because of the shape and it is
    also faster if it is smaller ...

    Args:
        s (Structure): pymatgen structure object
        return_non_metal_structure (bool): if true, also return non-metallic substructure
        white_list (list): elements that will be removed from the non-metall substructure
        allow_expansion (bool): if true, then attempt expansion
        primitive (bool):if true, we use the niggli reduced cell.
        other_metals (list): list of strings that are checked for in additions to transition metals and lanthanides

    Returns:
        structure object with only metallic substructure
    """
    crystal = s.copy()

    logger.debug('the formula of the crystal is %s', crystal.formula)
    if primitive:
        crystal = crystal.get_primitive_structure()
    if allow_expansion:
        crystal.make_supercell((2, 2, 2))

    elements_to_delete = [
        str(e)
        for e in crystal.species
        if not (e.is_transition_metal or e.is_rare_earth_metal or e.symbol in other_metals)
    ]

    logger.debug('the elements to delete are %s', elements_to_delete)

    if return_non_metal_structure:
        crystal2 = s.copy()
        metal_elements = (list(set(crystal2.symbol_set) - set(elements_to_delete)) + white_list)
        logger.debug('the metal elements are %s', metal_elements)
        crystal2.remove_species(metal_elements)
        crystal.remove_species(list(set(elements_to_delete)))
        return crystal, crystal2

    crystal.remove_species(list(set(elements_to_delete)))
    return crystal


def points_in_cylinder(pt1: np.ndarray, pt2: np.ndarray, r: float, q: np.ndarray):
    """
    https://stackoverflow.com/questions/47932955/how-to-check-if-a-3d-point-is-inside-a-cylinder

    ToDo: Fully vectorize it

    Args:
        pt1: landmark 1
        pt2: landmark 2
        r: radius of cylinder
        q: test point

    Returns:

    """
    vec = pt2 - pt1

    return len(
        np.where(dot(q - pt1, vec) >= 0 and dot(q - pt2, vec) <= 0 and norm(cross(q - pt1, vec)) <= r * norm(vec))[0])  # pylint:disable=chained-comparison


def points_in_cylinder_vectorized(pt1: np.ndarray, pt2: np.ndarray, r: float, q: np.ndarray) -> int:
    """
    Vectorized version of the check for the number of points from q in a cylinder spanned between
    pt1 and pt2 with radius r.

    Args:
        pt1 (np.ndarray): boundary point for cylinder
        pt2 (np.ndarray): boundary point for cylinder
        r (float): radius of cylinder
        q (np.ndarray): list of coordinates that should be checked for presence in cylinder

    Returns:
        int: number of points in this cylinder
    """
    vec = pt2 - pt1
    vec = vec.reshape(1, -1)
    pt1 = pt1.reshape(1, -1)
    pt2 = pt2.reshape(1, -1)

    return len(
        np.where((np.einsum('ij, ij->i', q - pt1, vec) >= 0) & (np.einsum('ij, ij->i', q - pt2, vec) <= 0) &
                 (np.linalg.norm(np.cross(q - pt1, vec), axis=1) <= r * norm(vec[0])))[0])


def get_dimensionality(s: Structure, method='Gorai') -> int:
    """
    Wrapper for the pymatgen dimensionality function.

    Args:
        s (pymatgen.Structure): Structure object
        method (str): method used for the dimensionality check. Currently only Gorai algorithm implemented

    Returns:

    """
    if method == 'Gorai':
        dim = dimensionality.get_dimensionality_gorai(s)
    else:
        raise NotImplementedError

    return dim


@jit
def get_points_within_tolerance(array: np.ndarray, tolerance: float) -> list:
    """
    Iterates trough points in array and find pairs (Euclidean pairwise distnace)

    Uses Numba.

    Args:
        array (np.ndarray):
        tolerance (float):

    Returns:
        list of tuples
    """
    res = []
    for i, coord0 in enumerate(array):
        for j, coord1 in enumerate(array):
            if i < j:
                if np.abs(coord0 - coord1) < tolerance:
                    res.append((i, j))
    return res


@jit
def check_parallel(vectors: list, threshold=0.1) -> bool:
    """

    Checks if all vectors in a list are parallel.

    Uses Numba.

    Args:
        vectors (list): list of np.ndarrays
        threshold (float): threhold (area) for parellelilty check

    Returns:
        boolean.
    """
    for i, vector0 in enumerate(vectors):
        for j, vector1 in enumerate(vectors):
            if i < j:
                if (dot(vector0, vector1) / (norm(vector0) * norm(vector1)) - 1 > threshold):
                    return False
    return True


@jit
def check_colinear(a: np.ndarray, b: np.ndarray, c: np.ndarray, tolerance: float = 0.05) -> bool:
    """
    Calculates the area of a triangle spanned by three points.
    If it is below the tolerance, the points are colinear.

    Uses Numba.

    Args:
        a (np.array):
        b (np.array):
        c (np.array):
        tolerance (float): tolerance factor for colinearity check

    Returns:
        boolean
    """
    a_copy = a / norm(a)
    c_copy = c / norm(c)
    b_copy = b / norm(b)
    area = norm(cross((a_copy - b_copy), (a_copy - c_copy)))
    return bool(area < tolerance)


def _find_landmarks(
    metal_radii: np.ndarray,
    metal_substructure: Structure,
    metal_coords: np.ndarray,
    tolerance_factor_2: float,
):
    """

    Args:
        metal_radii (np.ndarray):
        metal_substructure (pymatgen.Structure):
        metal_coords (np.ndarray):
        tolerance_factor_2 (float):

    Returns:

    """
    metal_tolerance = (metal_radii[metal_substructure[0].species_string] * tolerance_factor_2)
    landmarks_z = get_points_within_tolerance(metal_coords[:, 2], metal_tolerance)
    landmarks_y = get_points_within_tolerance(metal_coords[:, 1], metal_tolerance)
    landmarks_x = get_points_within_tolerance(metal_coords[:, 0], metal_tolerance)

    landmarks_z.sort(key=itemgetter(1))
    landmarks_y.sort(key=itemgetter(1))
    landmarks_x.sort(key=itemgetter(1))

    landmarks_z = set(landmarks_z)
    landmarks_y = set(landmarks_y)
    landmarks_x = set(landmarks_x)

    return landmarks_x, landmarks_y, landmarks_z


def _find_chain_landmarks(landmarks_z: list, landmarks_y: list, landmarks_x: list, metal_coords: np.ndarray) -> set:
    """
    Finds chain (3 metal atoms in a UC) of landmarks by intersecting sets.

    Args:
        landmarks_z (list): landmarks in z direction
        landmarks_y (list): landmarks in y direction
        landmarks_x (list): landmarks in x direction
        metal_coords (np.ndarray): metal coordinates

    Returns:
        set of chained landmark tuples.
    """
    landmarks = []
    for landmark_set in [landmarks_z, landmarks_y, landmarks_x]:
        for _, landmark0 in enumerate(landmark_set):
            for _, landmark1 in enumerate(landmark_set):
                if landmark0[1] == landmark1[0]:
                    #  tolerance = metal_radii[metal_substructure[
                    #      landmark0[1]].species_string] * tolerance_factor

                    a = metal_coords[landmark0[0]]
                    b = metal_coords[landmark0[1]]
                    c = metal_coords[landmark1[1]]

                    if check_colinear(a, b, c):
                        landmarks.append((landmark0[0], landmark1[1]))

    return set(landmarks)


def _get_metal_radii_dict(metal_substructure: Structure) -> dict:
    """
    Gets a dictionary with metal element strings as keys and their VdW radii as value to avoid frequent DB lookup.s

    Args:
        metal_substructure (pymatgen.Structure): Structure object containing only the metallic elements

    Returns:

    """
    metal_radii = {}
    for symbol in metal_substructure.symbol_set:
        metal_radii[symbol] = element(symbol).vdw_radius / 100
    return metal_radii


def _find_good_landmarks(
    landmarks: list,
    mc: np.ndarray,
    metal_substructure: Structure,
    metal_radii: dict,
    non_metal_coord: np.ndarray,
    tolerance_factor_1: float,
) -> list:
    """
    Iterates over landmarks and finds the ones that contain no non-metal in the cylinder spanned by them.

    Args:
        landmarks (list): list of tuples
        mc (np.ndarray): metal coordinates
        metal_substructure (pymatgen.Structure): pymatgen structure object with only metals
        metal_radii (dict): dictionary with the metal radii. Element strings as keys and radii as values
        non_metal_coord (np.ndarray): Non-metal coordinates
        tolerance_factor_1 (float): Tolerance factor (multiple of VdW radius)

    Returns:
        list of tuples that contain no non-metal between them.
    """
    good_landmarks = []
    metal_coords = mc.copy()
    for landmark in landmarks:
        p1 = metal_coords[landmark[0]]
        p2 = metal_coords[landmark[1]]

        connecting_vector = p2 - p1

        p2 += 5.0 * connecting_vector
        p1 -= 5.0 * connecting_vector

        tolerance = (metal_radii[metal_substructure[landmark[0]].species_string] * tolerance_factor_1)

        pic = points_in_cylinder_vectorized(p1, p2, tolerance, non_metal_coord)

        if pic == 0:
            good_landmarks.append(landmark)

    return good_landmarks


def _check_non_metal_in_cylinder(
    p1: np.ndarray,
    p2: np.ndarray,
    metal_substructure: Structure,
    metal_radii: dict,
    non_metal_coord: np.ndarray,
    tolerance_factor_1: float,
) -> int:
    """
    Check if the cylinder spanned between p1 and p2 with radius tolerance_factor_1 contains any non-metallic elements.

    Args:
        p1 (np.ndarray): boundary point for cylinder
        p2 (np.ndarray): boundary point for cylinder
        metal_substructure (pymatgen.Structure): pymatgen structure with only the metallic part
        metal_radii (dict): dictionary with the metal radii
        non_metal_coord (np.ndarray): coordinates of non-metal elements
        tolerance_factor_1 (float):

    Returns:

    """

    tolerance = metal_radii[metal_substructure[0].species_string] * tolerance_factor_1

    pic = points_in_cylinder_vectorized(p1, p2, tolerance, non_metal_coord)

    return pic


def _get_distance_checked_landmarks(
    good_landmarks: list,
    metal_radii: dict,
    metal_substructure: Structure,
    tolerance_factor_1: float,
    direction_dict: dict = None,
) -> list:
    """
    Loops over landmarks and checks that the maximum distance in a chain formed by them is smaller
    than the threshold.

    Args:
        good_landmarks (list): list of tuples of the landmarks
        metal_radii (dict): dictionary with the radius for each metal
        metal_substructure (pymatgen.Structure): pymatgen.Structure object with only the metallic parts
        tolerance_factor_1 (float): tolerance for the distance threshold
        direction_dict (dict): dictionary with the direction (value) each landmark (key) is parallel to

    Returns:
        list of landmarks that fullfill the distance criterion.

    """
    distance_checked_landmarks = []

    if direction_dict:
        a = metal_substructure.lattice.a
        b = metal_substructure.lattice.b
        c = metal_substructure.lattice.c

        lattice = np.array([a, b, c])

        metal_coord = np.abs(metal_substructure.cart_coords)

        for landmark in good_landmarks:
            direction = direction_dict[landmark]

            if isinstance(direction, (np.ndarray, list)):
                max_distance = (metal_radii[metal_substructure[landmark[0]].species_string] * tolerance_factor_1)

                direction = direction_dict[landmark]

                maximum_point_translation = lattice * direction

                # now, project the relevant coordinates
                coord_0 = metal_coord[landmark[0]] * direction
                coord_1 = metal_coord[landmark[1]] * direction

                distance_0 = np.min([norm(coord_0), norm(coord_1)])
                distance_1 = np.min([
                    norm(maximum_point_translation - coord_0),
                    norm(maximum_point_translation - coord_1),
                ])

                distance_2 = metal_substructure.get_distance(landmark[0], landmark[1])

                relevant_distance = np.max([distance_0, distance_1, distance_2])

                if relevant_distance < max_distance:
                    logger.debug('%s is within distance threshold of %s', landmark, max_distance)
                    distance_checked_landmarks.append(landmark)
                else:
                    logger.debug(
                        'the distance was too large between %s for a threshold of %s',
                        landmark,
                        max_distance,
                    )

    else:
        for landmark in good_landmarks:
            max_distance = (metal_radii[metal_substructure[landmark[0]].species_string] * tolerance_factor_1)

            distance = metal_substructure.get_distance(landmark[0], landmark[1])
            if distance < max_distance:
                logger.debug('%s is within distance threshold of %s', landmark, max_distance)
                distance_checked_landmarks.append(landmark)
            else:
                logger.debug(
                    'the distance was too large between %s for a threshold of %s',
                    landmark,
                    max_distance,
                )

    return distance_checked_landmarks


@jit
def _check_parallel_to_translation(landmark: tuple, metal_coords: np.ndarray, threshold: float = 0.2) -> np.ndarray:
    """
    Checks if the vector spanned by a landmark tuple is parallel to the translation direactions of
    the lattice

    Args:
        landmark (tuple): metal atom indices
        metal_coords (np.ndarray): coordinates of the metal substructures
        threshold (float): threshold for the paralellity check

    Returns:
        False if not parallel, the direction it is parallel to as np.ndarray if parallel
    """

    directions = [
        np.array([0, 0, 1]),
        np.array([0, 1, 0]),
        np.array([1, 0, 0]),
        np.array([1, 1, 1]),
        np.array([1, 1, 0]),
        np.array([1, 0, 1]),
        np.array([0, 1, 1]),
    ]

    vector0 = metal_coords[landmark[0]] - metal_coords[landmark[1]]

    for direction in directions:
        cp = cross(vector0, direction)

        if norm(cp) < threshold:
            return direction

    return False


def _get_landmarks_parallel_to_translation(landmarks: list, metal_coords: np.ndarray, threshold: float = 1.0):
    """
    Get landmarks that have same translation symmetry as the lattice

    Args:
        landmarks (list): list of landmarks, i.e metal atoms indices
        metal_coords (np.ndarray): coordinates of the metal substructures
        threshold (float): threshold (area) for the parallelity check

    Returns:
        List of tuples with good landmarks
        A dictionary if the direction a given tuple is parallel to
    """
    good_landmarks = []
    direction_dict = {}
    for landmark in landmarks:
        direction = _check_parallel_to_translation(landmark, metal_coords, threshold)

        if isinstance(direction, (np.ndarray, list)):
            good_landmarks.append(landmark)
            direction_dict[landmark] = direction

    return good_landmarks, direction_dict


def get_channels(  # pylint:disable=too-many-arguments, too-many-locals, too-many-statements, too-many-branches
    metal_substructure: Structure,
    non_metal_substructure: Structure,
    tolerance_factor_1: float = 2.0,
    tolerance_factor_2: float = 1.0,
    max_distance_factor: float = 3.0,
    mixed_metal: bool = True,
    all_in_channel: bool = False,
) -> int:
    """
    Assumption: Metal type is not important, i.e. will also detect mixed metal channels.

    Args:
        metal_substructure (Structure):
        non_metal_substructure (Structure):
        tolerance_factor_1 (float): factor by which the vdw radius of the element
            is multiplied build the cylinder of exclusion
        tolerance_factor_2 (float): factor by which the vdw radius of the element
            is multiplied to build the cylinder for metal overlap
        mixed_metal (bool): switch. If true, then also mixed metal channels are considered as channels.
        all_in_channel (bool): only true if ALL metals in the structure are part of a channel
        max_distance_factor (float): muliple of the VdW radius that is allowed as maximum distance between the
            metals of the chain

    Returns:
        empty list if no channel was found
        list of tuples of the metal atoms indices of the outer boundaries of the chain, if there is a chain

    """

    metal_coords = metal_substructure.cart_coords
    non_metal_coord = non_metal_substructure.cart_coords

    metal_radii = _get_metal_radii_dict(metal_substructure)

    a = metal_substructure.lattice.a
    b = metal_substructure.lattice.b
    c = metal_substructure.lattice.c

    logger.debug('Working on a cell with lattice constants a %s, b %s and c %s', a, b, c)
    directions = [
        np.array([0, 0, c]),
        np.array([0, b, 0]),
        np.array([a, 0, 0]),
        np.array([a, b, c]),
        np.array([a, b, 0]),
        np.array([a, 0, c]),
        np.array([0, b, c]),
    ]

    logger.debug('I found %s metal sites', len(metal_coords))
    logger.debug('I found %s non-metal sites', len(non_metal_coord))

    if len(metal_coords) == 0:
        logger.info('I found nothing that I consider to be a metal in this structure')
        return []

    if len(non_metal_coord) == 0:
        logger.error('I found nothing that is not a metal in this structure. '
                     'I am confused and will not investigate this structure any further and simply return 0')
        return []

    if not mixed_metal:
        # Basically check if the two elements of the tuples are the same element
        raise NotImplementedError

    if mixed_metal:  # pylint:disable=too-many-nested-blocks
        landmarks_x, landmarks_y, landmarks_z = _find_landmarks(metal_radii, metal_substructure, metal_coords,
                                                                tolerance_factor_2)

        # Now, we should find out if there is a "chain" of landmarks
        landmarks = _find_chain_landmarks(landmarks_z, landmarks_y, landmarks_x, metal_coords)

        logger.debug('the chained landmarks are %s', landmarks)

        joined_landmarks = landmarks_x | landmarks_y | landmarks_z

        if not (joined_landmarks and landmarks):
            # this is a case where we only have two elements in one structure
            # maybe we can generalize and only use this case, but I am not sure performance - wise because
            # the number of chained ones is generally lower
            logger.info('this is a case, where there is no chain')
            good_landmarks = _find_good_landmarks(
                joined_landmarks,
                metal_coords,
                metal_substructure,
                metal_radii,
                non_metal_coord,
                tolerance_factor_1,
            )

            if good_landmarks:
                logger.debug('checking now %s for translational symmetry', good_landmarks)
                good_landmarks, direction_dict = _get_landmarks_parallel_to_translation(good_landmarks, metal_coords)

        elif landmarks:
            logger.info('this is a case, where I found a chain within the UC')
            good_landmarks = _find_good_landmarks(
                landmarks,
                metal_coords,
                metal_substructure,
                metal_radii,
                non_metal_coord,
                tolerance_factor_1,
            )

            if good_landmarks:
                logger.debug('checking now %s for translational symmetry', good_landmarks)
                good_landmarks, direction_dict = _get_landmarks_parallel_to_translation(good_landmarks, metal_coords)

        logger.debug('the good landmarks are %s', good_landmarks)

        # if the user provides a max distance we use this for one more filtering step
        if max_distance_factor and good_landmarks:
            logger.debug('filtering for distance threshold')
            good_landmarks = _get_distance_checked_landmarks(
                good_landmarks,
                metal_radii,
                metal_substructure,
                max_distance_factor,
                direction_dict,
            )

        if not good_landmarks:
            logger.debug('checking now for pairs of metals with translational symmetry')
            good_landmarks = _find_good_landmarks(
                joined_landmarks,
                metal_coords,
                metal_substructure,
                metal_radii,
                non_metal_coord,
                tolerance_factor_1,
            )

            if good_landmarks:
                logger.debug('filtering %s for translational symmetry', good_landmarks)
                good_landmarks, direction_dict = _get_landmarks_parallel_to_translation(good_landmarks, metal_coords)

                if max_distance_factor and good_landmarks:
                    logger.debug('filtering for distance threshold')
                    good_landmarks = _get_distance_checked_landmarks(
                        good_landmarks,
                        metal_radii,
                        metal_substructure,
                        max_distance_factor,
                        direction_dict,
                    )

        if good_landmarks:
            good_sites = len(list(set(sum(good_landmarks, ()))))
            logger.debug('the number of good sites are %s', good_sites)
            if all_in_channel:
                if len(metal_coords) == good_sites:
                    return good_landmarks
                else:
                    return []
            else:
                return good_landmarks
        else:

            # Now we still to make sure that we do not have a case
            # where the replication is due the fact that there is no other element in front or
            # behind the metal. To do this, we use the translation symmetry of the lattice
            # to generate a cylinder for each metal element in each direction
            logger.debug('Now checking if there are empty connections between the metals and its replicas')
            good_metals = []
            for index, metal in enumerate(metal_coords):
                for i, direction in enumerate(directions):
                    p1 = metal - direction
                    p2 = metal + direction

                    pic = _check_non_metal_in_cylinder(
                        p1,
                        p2,
                        metal_substructure,
                        metal_radii,
                        non_metal_coord,
                        tolerance_factor_1,
                    )

                    if pic == 0:
                        good_metals.append((index, norm(direction)))

            if good_metals:
                # If the user provides a distance factor, we use this for one more filtering step
                if max_distance_factor:
                    distance_checked_landmarks = []
                    for metal in good_metals:
                        max_distance = (metal_radii[metal_substructure[0].species_string] * max_distance_factor)

                        if metal[1] < max_distance:
                            distance_checked_landmarks.append(metal)

                    if distance_checked_landmarks:
                        return [(i[0], i[0]) for i in distance_checked_landmarks]
                    else:
                        return []
                else:
                    return [(i[0], i[0]) for i in good_metals]
            else:
                return []


def get_metal_metadata(metal_substructure, index):
    symbol = metal_substructure[index].species_string
    el = element(symbol)

    element_properties = {
        'metal_electron_configuration': el.econf,
        'metal_electronegativity_pauling': el.en_pauling,
        'metal_number_electrons': el.electrons,
        'metal_group': el.group,
        'metal_block': el.block,
        'metal_dipole_polarizability': el.dipole_polarizability,
        'metal_electron_affinity': el.electron_affinity,
        'metal_mendeelev_number': el.mendeleev_number,
        'metal_metallic_radius': el.metallic_radius,
        'metal_metallic_radius_c12': el.metallic_radius_c12,
        'metal_oxidation_states': el.oxistates,
        'metal_vdw_radius': el.vdw_radius,
    }

    return element_properties


def get_structure_properties(structure: Structure, mode: str = 'all') -> dict:

    if mode == 'all':
        featurizer = MultipleFeaturizer([
            SiteStatsFingerprint.from_preset('CoordinationNumber_ward-prb-2017'),
            StructuralHeterogeneity(),
            ChemicalOrdering(),
            DensityFeatures(),
            MaximumPackingEfficiency(),
            SiteStatsFingerprint.from_preset('LocalPropertyDifference_ward-prb-2017'),
            StructureComposition(Stoichiometry()),
            StructureComposition(ElementProperty.from_preset('magpie')),
            StructureComposition(ValenceOrbital(props=['frac'])),
        ])
    else:
        # Calculate only those which do not need a Voronoi tesselation
        featurizer = MultipleFeaturizer([
            DensityFeatures(),
            StructureComposition(Stoichiometry()),
            StructureComposition(ElementProperty.from_preset('magpie')),
            StructureComposition(ValenceOrbital(props=['frac'])),
        ])

    X = featurizer.featurize(structure)

    matminer_dict = dict(list(zip(featurizer.feature_labels(), X)))

    matminer_dict['volume'] = structure.volume
    return matminer_dict


def get_channel_dimensionality(landmarks, metall_substructure, tolerance_factor=2):
    """
    Experimental. Will return np.nan if the channel is formed by only on metal in the UC.
    Otherwise it just checks if a landmark is part of channel intwo directions

    Args:
        landmarks
        metall_substructure:


    Returns:

    """
    s = metall_substructure.copy()

    flattened_landmarks = set(sum(landmarks, ()))
    planes = []

    metal_radii_dict = _get_metal_radii_dict(s)

    if len(flattened_landmarks) > 1:
        indices_to_delete = list(set(s.symbol_set) - flattened_landmarks)
        s.remove_sites(indices_to_delete)

        landmarks_x, landmarks_y, landmarks_z = _find_landmarks(metal_radii_dict, s, s.cart_coords, tolerance_factor)

        landmarks_xz = landmarks_x | landmarks_z
        if landmarks_xz:
            planes.append('xz')

        landmarks_xy = landmarks_x | landmarks_y
        if landmarks_xy:
            planes.append('xy')

        landmarks_yz = landmarks_y | landmarks_z
        if landmarks_yz:
            planes.append('yz')

        if len(landmarks_x) > 1:
            planes.append('xx')

        if len(landmarks_y) > 1:
            planes.append('yy')

        if len(landmarks_z) > 1:
            planes.append('zz')

    return planes
