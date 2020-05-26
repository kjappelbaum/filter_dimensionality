# filter_metal_channels

Experimental tool to filter out structures with 'metal channels' from structural databases.

If uses a pretty naive geometric analysis and I tried to implement some naive optimizations using numba's jit and vectorization.
If the structures come from the CSD, one can also use this library to parse for some bibliometric information.

## Usage

## Data

Some results for parsing the CSD release from May 2019 are in `data`.
