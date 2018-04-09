clark-biom
===========
.. image:: https://img.shields.io/travis/smdabdoub/clark-biom.svg?style=plastic
    :target: https://travis-ci.org/smdabdoub/clark-biom
    :alt: Travis CI build status

Create BIOM-format tables (http://biom-format.org) from CLARK output 
(http://clark.cs.ucr.edu/) for use with downstream tools such as
PhyloToAST (http://phylotoast.org).

Installation
------------

From PyPI:

.. code-block:: bash

    $ pip install clark-biom

From GitHub:

.. code-block:: bash

    $ pip install git+http://github.com/smdabdoub/clark-biom.git

From source:

.. code-block:: bash

    $ python setup.py install


Requirements
------------

- biom-format >= 2.1.5
- h5py >= 2.5.0 [optional]

Documentation
-------------

The program takes as input, one or more files output from CLARK's 
estimate_abundance tool. Each file is parsed and the counts for each OTU 
(operational taxonomic unit) are recorded, along with database ID (e.g. NCBI), 
and lineage. The extracted data are then stored in a BIOM table where each count
is linked to the Sample and OTU it belongs to. Sample IDs are extracted from the
input filenames (everything up to the '.' preceeding the extension).

The BIOM format currently has two major versions. Version 1.0 uses the 
JSON (JavaScript Object Notation) format as a base. Version 2.x uses the
HDF5 (Hierarchical Data Format v5) as a base. The output format can be
specified with the --fmt option. Note that a tab-separated (tsv) output
format is also available. The resulting file will not contain most of the
metadata, but can be opened by spreadsheet programs.

Version 2 of the BIOM format is used by default for output, but requires the
Python library 'h5py'. If the library is not installed, clark-biom will 
automatically switch to using version 1.0. Note that the output can 
optionally be compressed with gzip (--gzip) for version 1.0 and TSV files. 
Version 2 files are automatically compressed.

Currently the taxonomy for each OTU ID is stored as row metadata in the BIOM
table using the seven-level format used by QIIME and metaphlan: k__K, p__P, ... 
s__S. If you would like another format supported, please file an issue or send a
pull request (note the contribution guidelines).
::

    usage: clark-biom.py [-h] [-o OUTPUT_FP] [--fmt {hdf5,json,tsv}] [--gzip]
                         [--version] [-v]
                         clark_abd_tbl [clark_abd_tbl ...]

Usage examples
--------------

1. Basic usage with default parameters::

    $ clark-biom S1.txt S2.txt

  This produces a compressed BIOM 2.1 file: table.biom
  with sample IDs: S1, S2.

2. BIOM v1.0 output::

    $ clark-biom S1.txt S2.txt --fmt json

  Produces a BIOM 1.0 file: table.biom

3. Compressed TSV output::

    $ clark-biom S1.txt S2.txt --fmt tsv --gzip -o table.tsv

  Produces a TSV file: table.tsv.gz


Program arguments
-----------------

positional arguments::

    clark_abd_tbls          Abundance table files from estimate_abundance.sh

optional arguments::
    
      -o OUTPUT_FP, --output_fp OUTPUT_FP
                            Path to the BIOM-format file. By default, the table
                            will be in the HDF5 BIOM 2.x format. Users can 
                            output to a different format using the --fmt option.
                            The output can also be gzipped using the --gzip 
                            option. Default path is: ./table.biom
      --fmt {hdf5,json,tsv}
                            Set the output format of the BIOM table. Default is
                            HDF5.
      --gzip                Compress the output BIOM table with gzip. HDF5 BIOM
                            (v2.x) files are internally compressed by default,
                            so this option is ignored when specifying --fmt 
                            hdf5.
      --version             Print program's version number and exit
      -v, --verbose         Print status messages during program execution.
      -h, --help            Print this help message and exit
