Background
-----------

Channel detection algorithm version 1
**************************************
* The first step is to separate the structures into a metallic and a non-metallic
  substructure

  .. warning::

    Metallic currently means all transition metals and lanthanides. If needed, we could also
    give an option to specify this in more detail

* For efficiency reasons, we first check if there any "chains" with at least three metal
  sites that share a coordinates.

* For those chains, we check if there is any other element in a cylinder connecting the
  outer landmarks of the channel and if the points are co-linear and if they have translational
  symmetry. If selected, we also check the distance maximum distance between metals atoms of the
  chain

* If we do not find such a chain, we look for structures with two or one metal atoms with
  translational symmetry that do not have a non-whitelisted element within the channels
  and that fulfill the distance criteria (if selected)

* For the ones for which we detected channels, we run additional featurizations and sanity
  checks.


The algorithm has the following parameters:

* the radius of the cylinder for metal chain detection (as a fraction of the VdW radius)

* the radius of the cylinder for detection of other elements in the metal chain
  (as a fraction of the VdW radius)

* distance between two elements of the chain (as a fraction of the VdW radius)

* number of elements in the chain

Most of the routines are based on set intersections to identify atoms that share a dimension.

The output of the find main :code:`get_channels()` function is a list of tuples (which is of course
empty if there is no channel) and which contains the indices of the atoms spanning the channel if
there seems to be a channel.
The number of channels per unit cell is simply the length of this list.

.. note::

    For geometric reasons (and also speed if it is possible to reduce huge cells) we work with Niggli
    reduced structures


Featurization and checks
-------------------------

In any mode you will receive at least the following information for each structure that seems to have channels:

*   Stochiometric compostion descriptors as :math:`p` norms of the elemental fraction (cp. Ward PRB 2011)

    * 0-norm

    * 10-norm

    * 2-norm

    * 3-norm

    * 5-norm

    * 7-norm

*   All kinds of statistics derived from atomic properties

    * avg_dev AtomicWeight

    * avg_dev Column

    * avg_dev CovalentRadius

    * avg_dev Electronegativity

    * avg_dev GSbandgap

    * avg_dev GSmagmom

    * avg_dev GSvolume_pa

    * avg_dev MeltingT

    * avg_dev MendeleevNumber

    * avg_dev NUnfilled

    * avg_dev NValence

    * avg_dev NdUnfilled

    * avg_dev NdValence

    * avg_dev NfUnfilled

    * avg_dev NfValence

    * avg_dev NpUnfilled

    * avg_dev NpValence

    * avg_dev NsUnfilled

    * avg_dev NsValence

    * avg_dev Number

    * avg_dev Row

    * avg_dev SpaceGroupNumber

    * density

    * frac d valence electrons

    * frac f valence electrons

    * frac p valence electrons

    * frac s valence electrons

    * maximum AtomicWeight

    * maximum Column

    * maximum CovalentRadius
    * maximum Electronegativity
    * maximum GSbandgap
    * maximum GSmagmom
    * maximum GSvolume_pa
    * maximum MeltingT
    * maximum MendeleevNumber
    * maximum NUnfilled
    * maximum NValence
    * maximum NdUnfilled
    * maximum NdValence
    * maximum NfUnfilled
    * maximum NfValence
    * maximum NpUnfilled
    * maximum NpValence
    * maximum NsUnfilled
    * maximum NsValence
    * maximum Number
    * maximum Row
    * maximum SpaceGroupNumber
    * mean AtomicWeight
    * mean Column
    * mean CovalentRadius
    * mean Electronegativity
    * mean GSbandgap
    * mean GSmagmom
    * mean GSvolume_pa
    * mean MeltingT
    * mean MendeleevNumber
    * mean NUnfilled
    * mean NValence
    * mean NdUnfilled
    * mean NdValence
    * mean NfUnfilled
    * mean NfValence
    * mean NpUnfilled
    * mean NpValence
    * mean NsUnfilled
    * mean NsValence
    * mean Number
    * mean Row
    * mean SpaceGroupNumber
    * minimum AtomicWeight
    * minimum Column
    * minimum CovalentRadius
    * minimum Electronegativity
    * minimum GSbandgap
    * minimum GSmagmom
    * minimum GSvolume_pa
    * minimum MeltingT
    * minimum MendeleevNumber
    * minimum NUnfilled
    * minimum NValence
    * minimum NdUnfilled
    * minimum NdValence
    * minimum NfUnfilled
    * minimum NfValence
    * minimum NpUnfilled
    * minimum NpValence
    * minimum NsUnfilled
    * minimum NsValence
    * minimum Number
    * minimum Row
    * minimum SpaceGroupNumber
    * mode AtomicWeight
    * mode Column
    * mode CovalentRadius
    * mode Electronegativity
    * mode GSbandgap
    * mode GSmagmom
    * mode GSvolume_pa
    * mode MeltingT
    * mode MendeleevNumber
    * mode NUnfilled
    * mode NValence
    * mode NdUnfilled
    * mode NdValence
    * mode NfUnfilled
    * mode NfValence
    * mode NpUnfilled
    * mode NpValence
    * mode NsUnfilled
    * mode NsValence
    * mode Number
    * mode Row
    * mode SpaceGroupNumber
    * name
    * num_channels
    * packing fraction
    * range AtomicWeight
    * range Column
    * range CovalentRadius
    * range Electronegativity
    * range GSbandgap
    * range GSmagmom
    * range GSvolume_pa
    * range MeltingT
    * range MendeleevNumber
    * range NUnfilled
    * range NValence
    * range NdUnfilled
    * range NdValence
    * range NfUnfilled
    * range NfValence
    * range NpUnfilled
    * range NpValence
    * range NsUnfilled
    * range NsValence
    * range Number
    * range Row
    * range SpaceGroupNumber
    * vpa

Additionally, you will find the number of channel pairs in the unit cell and the index pairs of atoms
in the unit cell that the code thinks form a channel.

We then also provide some metal properties like:

* metal_electron_configuration

* metal_electronegativity_pauling

* metal.en_pauling

* metal_number_electrons

* metal_group

* metal_block

* metal_dipole_polarizability

* metal_electron_affinity

* metal_mendeelev_number

* metal_metallic_radius

* metal_metallic_radius_c12

* metal_oxidation_states

* metal_vdw_radius

Also, we provide some simple sanity checks on the structures,
which will lead to the following columns:

* hydrogens (the code checks with C with less then two non-carbon neighbors also have hydrogens of neighbors,
  this leads to false positives)

* unbound (this checks if there is unbound solvent, this is simply based on a check if there is something that is
  far from everything else)

* clashing (this is a simple check if there are clashing atoms)

* cif_error (this is a sanity check if pymatgen can read the cif file)

It will also try to get some bibliometric information using the crossref and scopus apis
after scraping the Web-CSD site:

* CSD DOI
* Paper DOI
* Paper abstract, authors, title, affiliation funding, pages, year, journal

If we found an abstract, we also parse if for some regexes:

* for UV/VIS
* for electronic/conduc
* for photo/luminescence

to quickly get an impression where there is useful experimental data.
