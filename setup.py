# -*- coding: utf-8 -*-
from __future__ import absolute_import

from setuptools import setup

import versioneer

setup(
    name='filter_metal_channels',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=['filter_metal_channels'],
    url='',
    license='MIT',
    entry_points={
        'console_scripts': ['screen_metal_channels=filter_metal_channels.run_screening:main'],
    },
    install_requires=[
        'numpy',
        'pymatgen',
        'tqdm',
        'pandas',
        'mendeleev',
        'matminer',
        'numba',
    ],
    extras_require={
        'testing': ['pytest', 'pytest-cov<2.6'],
        'docs': ['sphinx-rtd-theme', 'sphinxcontrib-bibtex'],
        'pre-commit': ['pre-commit', 'yapf', 'prospector', 'pylint', 'versioneer'],
        'bibliography': ['crossrefapi', 'pyscopus', 'ccdc'],
    },
    author='Kevin M. Jablonka, Andres Adolfo Ortega Guerrero, Maria Fumanal Quintana, Berend Smit',
    author_email='kevin.jablonka@epfl.ch',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Development Status :: 1 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    description='tool to find MOFs with metal channels',
)
