Quickstart
-----------

Installation
............

The package automatically installs a range of dependencies, we recommend to use it in a virtual environment, e.g. a
conda environment, i.e.

::

    conda create --name filter_metal_channels python=3.7
    pip install git+https://github.com/kjappelbaum/filter_dimensionality.git


We cannot automatically install the CSD API. You need to download the CSD System with
your license and then download the latest python API and install it with something like

::

    conda install -c /home/kevin/Downloads/ccdc_conda_channel  Pillow six lxml numpy matplotlib
    conda install -c /home/kevin/Downloads/ccdc_conda_channel  csd-python-api

The most recent conda installer is available at `http://cdn.ccdc.cam.ac.uk/2019/API_3_2.1.0/csd-python-api-2.1.0-linux-64-py3.7-conda.zip?Expires=1577833200&Signature=puSXxvRDtTFtK~3xvAu8cnj-zqBwcofw7Gb-r6H8FnYjt88FgOaKT~AAkF731di9lvAx2qrC~893gV6xflmRrZaef6AeIIPx4RAzWT6mIPJDOFWnNf9eFefbv9eibrCmcLXteiqqqOVtSjjTnOR-zzJaRbWbRANNrpq~q4fEDEm3wF7ZtxOgZB-uP~r7Ors3u1wgWEzto1wAeTJ3EgzHYcmlqDfhYJec8YN5iHsjaCk~94UpF1RuRzbtuU3NVNFXf770qIUKQr5tF~-BeXmx6EKqU-9Xj-K2qcJgqt-OWGhfSZa5PrrOqYFAjokZWed70s2O-CrZqsChg21O-7WhXg__&Key-Pair-Id=APKAJKVKVM6B4553JOZA`

Then, they also want you to add a path. It is best to just follow the `installation instructions <https://downloads.ccdc.cam.ac.uk/documentation/API/installation_notes.html>`_.
For Linux, it is

::

    export CSDHOME=/home/kevin/CCDC/CSD_2019

Note that this Python API works only in Python 3.7 and it is hard to fix it since they hardcoded a lot of library paths.


Usage
......

Console script
```````````````
The easiest way to use the package is via the console script :code:`screen_metal_channels` that is automatically installed.
It can take the following arguments:

* :code:`-f`: path to the folder that contains the structures (required)
* :code:`-ext`: extension of the files in the folder (default: :code:`.cif`)
* :code:`-out`: filename for the output (:code:`.csv` file with features)
* :code:`-njobs`: maximum number of workers (default :code:`2`)
* :code:`-featuremode`: allows to select if you want to compute all features (also those that require a Voronoi tesselation)
  or only cheap ones.

I.e. it should be generally enough to simply run

::

    screen_metal_channels -f tests/structures

You should then be informed about the progress of the screening by a progress bar.


API
````
The API deserves no design award. Still you can initialize a :code:`ScreenMetalChannel` class object with a list
of structure paths or use the constructor method :code:`.from_folder()` to create this object from a folder.

The :code:`run()` function of the class object then runs the screening and takes the same arguments as the console script,
except for the output name as it simply returns a pandas DataFrame object.
