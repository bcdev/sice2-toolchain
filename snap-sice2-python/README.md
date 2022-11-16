SICE2 SNAP plugins README
=========================

### Contents of SICE2 SNAP plugins delivery v1.0:
- esa_snappy nbm
- esa_snappy_conf.bash, esa_snappy_conf.bat (from this package)
- sice2 nbm
- README draft
- example graphs
- example gpt bat
- example gpt bash

### Installation
- if not already done, install latest SNAP release 9.0
- in SNAP desktop, update all modules to latest available versions
- install esa_snappy and sice nbms via SNAP desktop plugin manager
- install Python >= 3.7 (Windows), >= 3.8 (Linux), e.g. miniconda. See more details below.
- install all additionally required Python libraries. See imports in pySICEv21 code to get an idea what is needed.
- install new 'esa_snappy':
  - make sure environment variable SNAP_HOME is set to SNAP installation directory, e.g. SNAP_HOME=/home/jason/snap_9.0.0
  - check if old 'snappy' was installed in SNAP earlier: if directory <home>/.snap/snap-python/snappy exists, remove it
  - configure new 'esa_snappy': save esa_snappy_conf.bash or esa_snappy_conf.bat to arbitrary directory and run there:
  - --> on Windows: esa_snappy_conf.bat </path/to/python_executable>
  - --> on Linux: esa_snappy_conf.bash </path/to/python_executable>

### Run the plugin
- after successful installation, you should see two operators in SNAP gpt help:
  - gpt -h py_sice2_v21_op
  - gpt -h py_sice2_v21_tifdirs_op
  - --> The first operator can be regarded as the 'main' SICE2 SNAP operator. It expects an OLCI L1 product as input.
Subsets in other formats as SAFE are also accepted. As we know from S3Snow, an IdePix cloud mask product can be optionally be provided 
if given on the raster as the OLCI L1b. The output is a single product containing all the variables which pySICEv21 writes into single
tif files. As for most SNAP operators, the user can select from many possible output formats. 
This operator is probably the one to be considered for external users.
  - --> The second operator basically mirrors the functionality of the pySICEv21 pure Python version. As pySICEv21, this operator takes
as input single tif files which were generated before by the pySICE preprocessing chain. However, different from pySICEv21,
the output is written again into a single target product containing all variables. If needed, a split into single products can be
done easily with the SNAP subset operator ('band subset').
  - the two operators are also available in SNAP desktop via menu 'Optical --> Thematic Land Processing'

### Performance
- py_sice2_v21_op: Processing of a full OLCI L1 FR scene (4865x4091 pixel) with all spectral variables being written: ~75min. Thus, it is recommended to create
appropriate subsets in advance for the regions of interest.
- py_sice2_v21_tifdirs_op: Processing of Greenland mosaics for one day: ~7min (comparable with performance of pySICEv21 pure Python)

### Verification
- for two days of Greenland mosaics, the SNAP plugin results were compared with pySICEv21 pure Python for all output variables.
The SNAP plugin results reproduce the pure Python results (within numerical tolerances). 

### Known issues and current limitations
- As reported, a major (blocking) issue was that the 'snappy' module which comes with the latest SNAP 9 release does not allow the use
of Python 'xarray' module due to a name conflict. This issue was overcome for SICE2 with the introduction of a renamed module 'esa_snappy'.
For varios reasons this is provided as an external plugin (see above) rather than a SNAP core component. I.e., this approach ensures
backward compatibility if a user needs to switch back to the original 'snappy' if he has existing Python code which imports snappy.  
- The use of xarray in pySICEv21 requires a Python version >= 3.7. On Windows, the SNAP plugins were run successfully with 3.7 and 3.9.
On Linux, the 'jpy' Java/Python bridge module used by snappy and esa_snappy does not seem to work with Python 3.7 (several jpy versions
have been tested. There seems to be an issue with the underlying native C code. For the moment, we recommend to use Python 3.9 on Linux
which has been tested successfully.  
- the SNAP plugins cannot yet be used on Mac OS platforms (no jpy wheels generated yet)
- there is no help documentation yet available in SNAP desktop (will be added if the final layout of the operstors is agreed within the team)
- In SNAP 9 there is a performance issue with the tile-based computation. In an unpredictable manner, a given tile is often computed
multiple times. This issue will be adressed with high priority in the upcoming SNAP maintenance/evolution phase (kick-off Nov 16/17th).
A current workaround is to use the SICE plugins not from SNAP desktop, but ONLY from the command line with a gpt graph (as illustrated in the examples).
- IdePix cloud mask usage when processing from L1: although we reached significant improvements during the S3Snow project, the detection of
clouds over snow remains extremely difficult and is still a challenge. Thus, the cloud mask should always be interpreted with care.

*Have fun!*
