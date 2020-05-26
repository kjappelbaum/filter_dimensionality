#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import absolute_import

import logging
import math
import os
import re
from glob import glob
from pathlib import Path

import numpy as np
from crossref.restful import Works
from pyscopus import Scopus

from ccdc import io
from numba import double, jit

__author__ = 'Kevin M. Jablonka'
_license__ = 'MIT'
__version__ = '0.1.0'
__email__ = 'kevin.jablonka@epfl.ch'
__status__ = 'Dev'

scopus = Scopus('97d9d219e9efce0823123c557ec179a9')

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)


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


def get_structure_list(directory: str, extension: str = 'cif') -> list:
    """

    Args:
        directory:
        extension:

    Returns:

    """
    logger.info('getting structure list')
    if extension:
        structure_list = glob(os.path.join(directory, ''.join(['*.', extension])))
    else:
        structure_list = glob(os.path.join(directory, '*'))
    return structure_list


def get_bibliometric_information(cif_path):
    """
    Assumes that the filename is the CSD-key.
    Assumes you are in the EPFL VPN.
    Uses my scopus API key.

    Args:
        cif_path:

    Returns:

    """

    stem = Path(cif_path).stem.upper()

    deposition_number = None
    title = None
    citations = None
    url = None
    abstract = None
    funding = None
    doi_ccsd = None
    doi_paper = None
    pages = None
    journal = None
    affilations = None
    year = None
    disorder_csd = None
    authors = None
    uv_regex_result = None
    photo_regex_result = None
    electronic_regex_result = None
    csd_remarks = None
    csd_has_disorder = None
    chemical_name = None
    formula = None

    try:
        csd_reader = io.EntryReader('CSD')

        reader = csd_reader.entry(stem)
        disorder_csd = reader.disorder_details
        doi_paper = reader.publication.doi
        csd_has_disorder = reader.has_disorder
        doi_ccsd = reader.doi

        deposition_number = reader.ccdc_number
        chemical_name = reader.chemical_name
        formula = reader.formula

        uv_regex_list = ['uv', 'uv-vis', 'vis']
        uv_regex = re.compile('|'.join(uv_regex_list), re.IGNORECASE)
        uv_regex_result = False

        photo_regex_list = ['photo', 'absorp', 'light', 'lumin']
        photo_regex = re.compile('|'.join(photo_regex_list), re.IGNORECASE)
        photo_regex_result = False

        electronic_regex_list = ['electronic', 'conduc']
        electronic_regex = re.compile('|'.join(electronic_regex_list), re.IGNORECASE)
        electronic_regex_result = False

        if doi_paper:
            logger.info('found DOI %s', doi_paper)
            works = Works()
            query_res = works.doi(doi_paper)
            if query_res:

                if 'title' in query_res.keys():
                    if len(query_res['title']) > 0:
                        title = query_res['title'][-1]
                    else:
                        title = query_res['title']

                if 'is-referenced-by-count' in query_res.keys():
                    citations = query_res['is-referenced-by-count']

                if 'link' in query_res.keys():
                    if len(query_res['link']) > 0:
                        url = query_res['link'][0]['URL']
                    else:
                        url = query_res['link']

                if 'abstract' in query_res.keys():
                    abstract = query_res['abstract']

                funding_sublis = []
                if 'funder' in query_res.keys():
                    for f in query_res['funder']:
                        funding_sublis.append(f['name'])

                    funding = funding_sublis
                else:
                    funding = None

                if 'page' in query_res.keys():
                    pages = query_res['page']

                if 'container-title' in query_res.keys():
                    if len(query_res['container-title']) > 0:
                        journal = query_res['container-title'][0]
                    else:
                        journal = query_res['container-title']

                if 'author' in query_res.keys():
                    author_sublist = []
                    for a in query_res['author']:
                        author_sublist.append(a['family'])
                    authors = author_sublist

                    affiliation_sublist = []
                    for a in query_res['author']:
                        affiliation_sublist.append(a['affiliation'])

                if abstract is None:
                    search_result = scopus.search(doi_paper)
                    if len(search_result) > 0:
                        try:
                            abstract = scopus.retrieve_abstract(search_result['scopus_id'].values[0])['abstract']
                        except Exception:
                            abstract = None

        if abstract:
            electronic_regex_result = re.findall(electronic_regex, abstract)
            if len(electronic_regex_result) > 0:
                electronic_regex_result = True

            uv_regex_result = re.findall(uv_regex, abstract)
            if len(uv_regex_result) > 0:
                uv_regex_result = True

            photo_regex_result = re.findall(photo_regex, abstract)
            if len(photo_regex_result) > 0:
                photo_regex_result = True

    except Exception:
        logger.info('Could not retrieve CSD info')

        result_dict = {
            'deposition_number': deposition_number,
            'title': title,
            'csd_abbrv': stem,
            'citations': citations,
            'url': url,
            'abstract': abstract,
            'funding': funding,
            'doi_ccsd': doi_ccsd,
            'doi_paper': doi_paper,
            'formula': formula,
            'pages': pages,
            'journal': journal,
            'affilations': affilations,
            'year': year,
            'disorder_csd': disorder_csd,
            'authors': authors,
            'remarks': csd_remarks,
            'csd_has_disorder': csd_has_disorder,
            'chemical_name': chemical_name,
            'uv_regex_result': uv_regex_result,
            'photo_regex_result': photo_regex_result,
            'electronic_regex_result': electronic_regex_result,
        }

    else:
        result_dict = {
            'deposition_number': deposition_number,
            'title': title,
            'csd_abbrv': stem,
            'citations': citations,
            'url': url,
            'abstract': abstract,
            'funding': funding,
            'doi_ccsd': doi_ccsd,
            'doi_paper': doi_paper,
            'formula': formula,
            'pages': pages,
            'journal': journal,
            'affilations': affilations,
            'year': year,
            'disorder_csd': disorder_csd,
            'authors': authors,
            'remarks': csd_remarks,
            'csd_has_disorder': csd_has_disorder,
            'chemical_name': chemical_name,
            'uv_regex_result': uv_regex_result,
            'photo_regex_result': photo_regex_result,
            'electronic_regex_result': electronic_regex_result,
        }

    return result_dict


def isnotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False
