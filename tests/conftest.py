#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os

import pytest
from pymatgen import Structure

THIS_DIR = os.getcwd()


@pytest.fixture(scope='module')
def get_channel_structure_objects():
    s1 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'ZSTU-1.cif'))
    s2 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'RuTBPZn.cif'))
    s3 = Structure.from_file(os.path.join(THIS_DIR, 'structures', '1832956.cif'))
    s4 = Structure.from_file(os.path.join(THIS_DIR, 'structures', '893545.cif'))
    # s5 = Structure.from_file(os.path.join(THIS_DIR, 'structures',
    #                                      '775691.cif')) they are quite fa

    s8 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AHOKIR01_clean_min_charges.cif'))

    s9 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'RuTBPZn_Al.cif'))
    # s10 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'pegcip.cif'))
    # s10 = Structure.from_file(os.path.join(THIS_DIR, 'structures',
    #                                      'AMUCOB.cif'))
    return [s1, s2, s3, s4, s8, s9]


@pytest.fixture(scope='module')
def get_channels_far():
    s5 = Structure.from_file(os.path.join(THIS_DIR, 'structures', '775691.cif'))
    s6 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'mof-806.cif'))
    s7 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'vilwix_zn.cif'))
    s9 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'xempek.cif'))
    return [s5, s6, s7, s9]


@pytest.fixture(scope='module')
def get_only_metal_structures_channel():
    s1 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'RuTBPZn_only_metal.cif'))
    s2 = Structure.from_file(os.path.join(THIS_DIR, 'structures', '1832956_only_metal.cif'))

    return s1, s2


@pytest.fixture(scope='module')
def get_all_structures():
    structures = []
    for structure in os.listdir(os.path.join(THIS_DIR, 'structures')):
        structure_path = os.path.join(THIS_DIR, 'structures', structure)
        s = Structure.from_file(structure_path)
        structures.append(s)
    return structures


@pytest.fixture(scope='module')
def get_only_metal_structures_no_channel():
    raise NotImplementedError


@pytest.fixture(scope='module')
def get_no_channel_structures():  # pylint:disable=too-many-locals
    s1 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AMILUE_clean_min_charges.cif'))
    s2 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'ABUWOJ_clean_min_charges.cif'))
    s3 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'ATOXEN_clean_min_charges.cif'))
    s4 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AFITIT_clean_min_charges.cif'))
    s6 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BUVXOG_clean_min_charges.cif'))
    s7 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BOWQAG_clean_min_charges.cif'))
    s8 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AROFET_clean_min_charges.cif'))
    s9 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AVAQIX_clean_min_charges.cif'))
    s10 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BETGAK_clean_min_charges.cif'))
    s11 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BUVXOG_clean_min_charges.cif'))
    s12 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BONWAD_clean_min_charges.cif'))
    s13 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BOJCIN_clean_min_charges.cif'))
    s14 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'ACOLIP_clean_min_charges.cif'))
    s15 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AMIKOX_clean_min_charges.cif'))
    s16 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BARZAW_clean_min_charges.cif'))
    s17 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'ANUGIA_clean_min_charges.cif'))
    s18 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'AVIMOI_clean_min_charges.cif'))
    s19 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BENXUP_clean_min_charges.cif'))
    s20 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BERGAI_clean_min_charges.cif'))
    s21 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'BONWIL_clean_min_charges.cif'))

    return [
        s1,
        s2,
        s3,
        s4,
        s6,
        s7,
        s8,
        s9,
        s10,
        s11,
        s12,
        s13,
        s14,
        s15,
        s16,
        s17,
        s18,
        s19,
        s20,
        s21,
    ]


@pytest.fixture(scope='module')
def get_channel_no_tm():
    s1 = Structure.from_file(os.path.join(THIS_DIR, 'structures', 'vilwix.cif'))
    return s1
