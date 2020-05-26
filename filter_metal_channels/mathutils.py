# -*- coding: utf-8 -*-

from __future__ import absolute_import

import math

import numpy as np
from numba import double, jit


@jit(nopython=True)
def norm(vec):
    """ Calculate the norm of a 3d vector. """
    return math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])


@jit(nopython=True)
def dot(vec1, vec2):
    """ Calculate the dot product of two 3d vectors. """
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]


@jit
def cross(vec1, vec2):
    """ Calculate the cross product of two 3d vectors. """
    result = np.zeros(3)
    return cross_(vec1, vec2, result)


@jit(nopython=True)
def cross_(vec1, vec2, result):
    """ Calculate the cross product of two 3d vectors. """
    a1, a2, a3 = double(vec1[0]), double(vec1[1]), double(vec1[2])
    b1, b2, b3 = double(vec2[0]), double(vec2[1]), double(vec2[2])
    result[0] = a2 * b3 - a3 * b2
    result[1] = a3 * b1 - a1 * b3
    result[2] = a1 * b2 - a2 * b1
    return result
