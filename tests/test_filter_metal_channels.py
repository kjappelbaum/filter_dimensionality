#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import, print_function

import os

import numpy as np

from filter_metal_channels import filter_metal_channels

THIS_DIR = os.getcwd()


def test_get_metal_substructure(get_five_structure_objects):
    s1, s2, s3, _, _ = get_five_structure_objects

    s1_metall = filter_metal_channels.get_metallic_substructure(s1, return_non_metal_structure=False)
    s2_metall = filter_metal_channels.get_metallic_substructure(s2, return_non_metal_structure=False)
    s3_metall = filter_metal_channels.get_metallic_substructure(s3, return_non_metal_structure=False)

    s1_metall1, s1_nonmetall = filter_metal_channels.get_metallic_substructure(s1, return_non_metal_structure=True)

    for site in s1_metall.species:
        assert site.is_transition_metal

    for site in s2_metall.species:
        assert site.is_transition_metal

    for site in s3_metall.species:
        assert site.is_transition_metal

    for site in s1_metall1.species:
        assert site.is_transition_metal

    for site in s1_nonmetall.species:
        assert not site.is_transition_metal


def test_points_in_cylinder():
    p1 = np.array([0, 0, 0])
    p2 = np.array([0, 0, 1])
    probe = np.array([0.2, 0.2, 0.8])
    assert filter_metal_channels.points_in_cylinder(p1, p2, 1, probe) == 1

    probe = np.array([1, 1, 0.8])
    assert filter_metal_channels.points_in_cylinder(p1, p2, 0.5, probe) == 0


def test_points_in_cylinder_vectorized():
    p1 = np.array([0, 0, 0])
    p2 = np.array([0, 0, 1])
    probe = np.array([[0.2, 0.2, 0.8], [0.1, 0.1, 0.6], [0.1, 0.1, 0.1], [0.2, 0.1, 0.6]])
    assert filter_metal_channels.points_in_cylinder_vectorized(p1, p2, 1, probe) == 4

    probe = np.array([[5, 0.2, 0.8], [0.1, 0.1, 0.6], [0.1, 0.1, 0.1], [0.2, 0.1, 0.6]])
    assert filter_metal_channels.points_in_cylinder_vectorized(p1, p2, 1, probe) == 3


def test_get_dimensionality(get_only_metal_structures_channel):
    s1, s2 = get_only_metal_structures_channel

    dim = filter_metal_channels.get_dimensionality(s1)
    assert dim == 1

    dim = filter_metal_channels.get_dimensionality(s2)
    assert dim == 1


def test_random_structures(get_all_structures):
    # As a quick check for me if anything weird happens in any part of the code
    for i, s in enumerate(get_all_structures):
        print(f'working on {i}')
        s1_metall1, s1_nonmetall = filter_metal_channels.get_metallic_substructure(s, return_non_metal_structure=True)

        num_channels = filter_metal_channels.get_channels(s1_metall1, s1_nonmetall)

        print(f'found {num_channels} channels in {i}')


def test_get_channels(get_channel_structure_objects):
    for i, s in enumerate(get_channel_structure_objects):
        print(i)
        s1_metall1, s1_nonmetall = filter_metal_channels.get_metallic_substructure(s, return_non_metal_structure=True)

        num_channels = filter_metal_channels.get_channels(s1_metall1, s1_nonmetall)

        print(f'found {num_channels} channels')
        assert len(num_channels) >= 1


def test_get_channel_no_tm(get_channel_no_tm):
    s1 = get_channel_no_tm

    s1_metall1, s1_nonmetall = filter_metal_channels.get_metallic_substructure(s1, return_non_metal_structure=True)

    num_channels = filter_metal_channels.get_channels(s1_metall1, s1_nonmetall)

    assert len(num_channels) == 0


def test_get_channel_no_channel(get_no_channel_structures):

    for i, s in enumerate(get_no_channel_structures):
        print(i)
        s1_metall1, s1_nonmetall = filter_metal_channels.get_metallic_substructure(s, return_non_metal_structure=True)

        num_channels = filter_metal_channels.get_channels(s1_metall1, s1_nonmetall)

        print(f'found {num_channels} channels')

        assert len(num_channels) < 1
