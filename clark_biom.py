#!/usr/bin/env python
# coding: utf-8
"""
Description
-----------
Create BIOM-format tables (http://biom-format.org) from 
CLARK output (http://clark.cs.ucr.edu/).
"""
from __future__ import absolute_import, division, print_function

import argparse
from collections import OrderedDict
import csv
from datetime import datetime as dt
from gzip import open as gzip_open
import os.path as osp
import sys
from textwrap import dedent as twdd

from biom.table import Table
import numpy as np

try:
    import h5py
    HAVE_H5PY = True
except ImportError:
    HAVE_H5PY = False

__author__ = "Shareef M. Dabdoub"
__copyright__ = "Copyright 2018, Shareef M. Dabdoub"
__credits__ = ["Shareef M. Dabdoub", "Sukirth Ganesan", "Purnima Kumar"]
__license__ = ""
__url__ = "http://github.com/smdabdoub/clark-biom"
__maintainer__ = "Shareef M. Dabdoub"
__email__ = "dabdoub.2@osu.edu"
__version__ = '1.0.0'


field_names = ["Name", "TaxID", "Lineage", "Count",
                "Proportion_All(%)", "Proportion_Classified(%)"]
ranks = ["k", "p", "c", "o", "f", "g", "s"]



def tax_fmt(lineage, name):
    """
    Create a string representation of a taxonomic hierarchy.

    :type lineage: str
    :param lineage: Semi-colon separated list of taxnonmic levels
    :type name: int
    :param name: The standard "Scientific Name" of an organism/level.
    
    >>> tax_fmt("Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus", "Bacillus subtilis")
    ['k__Bacteria', 'p__Firmicutes', 'c__Bacilli', 'o__Bacillales', 'f__Bacillaceae', 'g__Bacillus', 's__subtilis']
    """
    lineage = lineage.split(';')
    tax = [r + "__" + l for r, l in zip(ranks, lineage)]

    if len(tax) == len(ranks)-1:
        name = name.split(' ')
        if lineage[-1] == name[0] and len(name) > 1:
            sp = ' '.join(name[1:])
            tax.append("s__" + sp)

    return tax


def parse_clark_abundance_tbl(data, store_pct=False):
    """
    Parse a single output file from estimate_abundance.sh. Return a list
    of counts at each of the acceptable taxonomic levels, and a list of 
    NCBI IDs and a formatted string representing their taxonomic hierarchies.

    :type data: list of dicts
    :param data: Each entry contains a dict keyed on the global 'field_names'
                 that corresponds to the columns of the CLARK abundance table
                 format.
    """
    # map between NCBI taxonomy IDs and the string rep. of the hierarchy
    taxa = OrderedDict()
    # the master collection of read counts (keyed on NCBI ID)
    counts = OrderedDict()

    for entry in data:
        if entry['TaxID'] == "UNKNOWN":
            continue
        taxa[entry['TaxID']] = tax_fmt(entry['Lineage'], entry['Name'])

        if store_pct:
            counts[entry['TaxID']] = float(entry['Proportion_Classified(%)'])
        else:
            counts[entry['TaxID']] = int(entry['Count'])
    
    return counts, taxa


def process_samples(clark_abd_fps, store_pct=False):
    """
    Parse all clark abundance tables into sample counts dict
    and store global taxon id -> taxonomy data
    """
    taxa = OrderedDict()
    sample_counts = OrderedDict()
    for clark_fp in clark_abd_fps:
        if not osp.isfile(clark_fp):
            raise RuntimeError("ERROR: File '{}' not found.".format(clark_fp))

        # use the clark abundance table filename as the sample ID
        sample_id = osp.splitext(osp.split(clark_fp)[1])[0]

        with open(clark_fp, "rt") as cf:
            try:
                cdr = csv.DictReader(cf, fieldnames=field_names)
                data = [entry for entry in cdr][1:]
            except OSError as oe:
                raise RuntimeError("ERROR: {}".format(oe))

        scounts, staxa = parse_clark_abundance_tbl(data, store_pct=store_pct)

        # update master records
        taxa.update(staxa)
        sample_counts[sample_id] = scounts

    return sample_counts, taxa


def create_biom_table(sample_counts, taxa):
    """
    Create a BIOM table from sample counts and taxonomy metadata.

    :type sample_counts: dict
    :param sample_counts: A dictionary of dictionaries with the first level
                          keyed on sample ID, and the second level keyed on
                          taxon ID with counts as values.
    :type taxa: dict
    :param taxa: A mapping between the taxon IDs from sample_counts to the
                 full representation of the taxonomy string. The values in
                 this dict will be used as metadata in the BIOM table.
    :rtype: biom.Table
    :return: A BIOM table containing the per-sample taxon counts and full
             taxonomy identifiers as metadata for each taxon.
    """
    data = [[0 if taxid not in sample_counts[sid] else sample_counts[sid][taxid] 
              for sid in sample_counts] 
                for taxid in taxa]
    data = np.array(data, dtype=int)
    tax_meta = [{'taxonomy': taxa[taxid]} for taxid in taxa]
    
    gen_str = "clark-biom v{} ({})".format(__version__, __url__)

    return Table(data, list(taxa), list(sample_counts), tax_meta, 
                 type="OTU table", create_date=str(dt.now().isoformat()),
                 generated_by=gen_str, input_is_dense=True)


def write_biom(biomT, output_fp, fmt="hdf5", gzip=False):
    """
    Write the BIOM table to a file.

    :type biomT: biom.table.Table
    :param biomT: A BIOM table containing the per-sample OTU counts and metadata
                  to be written out to file.
    :type output_fp str
    :param output_fp: Path to the BIOM-format file that will be written.
    :type fmt: str
    :param fmt: One of: hdf5, json, tsv. The BIOM version the table will be
                output (2.x, 1.0, 'classic').
    """
    opener = open
    mode = 'w'
    if gzip and fmt != "hdf5":
        if not output_fp.endswith(".gz"):
            output_fp += ".gz"
        opener = gzip_open
        mode = 'wt'

    # HDF5 BIOM files are gzipped by default
    if fmt == "hdf5":
        opener = h5py.File

    with opener(output_fp, mode) as biom_f:
        if fmt == "json":
            biomT.to_json(biomT.generated_by, direct_io=biom_f)
        elif fmt == "tsv":
            biom_f.write(biomT.to_tsv())
        else:
            biomT.to_hdf5(biom_f, biomT.generated_by)

    return output_fp


def write_otu_file(otu_ids, fp):
    """
    Write out a file containing only the list of OTU IDs from the CLARK
    results. One line per ID.

    :type otu_ids: list or iterable
    :param otu_ids: The OTU identifiers that will be written to file.
    :type fp: str
    :param fp: The path to the output file.
    """
    fpdir = osp.split(fp)[0]

    if not fpdir == "" and not osp.isdir(fpdir):
        raise RuntimeError("Specified path does not exist: {}".format(fpdir))

    with open(fp, 'wt') as outf:
        outf.write('\n'.join(otu_ids))


def handle_program_options():
    descr = """\
    Create BIOM-format tables (http://biom-format.org) from CLARK output 
    (http://clark.cs.ucr.edu/).

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

    Usage examples
    --------------

    1. Basic usage with default parameters::

        $ clark-biom S1.csv S2.csv

      This produces a compressed BIOM 2.1 file: table.biom
      with sample IDs: S1, S2.

    2. Process multiple samples from multiple groups::

        $ clark-biom groupA/*.csv groupB/*.csv -o groupsAB.biom

    3. BIOM v1.0 output::

        $ clark-biom S1.csv S2.csv --fmt json

      Produces a BIOM 1.0 file: table.biom

    4. Compressed TSV output::

        $ clark-biom S1.csv S2.csv --fmt tsv --gzip -o table.tsv

      Produces a TSV file: table.tsv.gz


    Program arguments
    -----------------"""

    parser = argparse.ArgumentParser(description=twdd(descr),
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('clark_abd_tbls', nargs='+', metavar="TABLE-FILE",
                        help="Result file from estimate_abundance.sh.")
    parser.add_argument('-o', '--output_fp', default="table.biom",
                        metavar="COMBINED-OUTPUT-FILE",
                        help="Path to the BIOM-format file. By default, the "
                        "table will be in the HDF5 BIOM 2.x format. Users can "
                        "output to a different format using the --fmt option. "
                        "The output can also be gzipped using the --gzip"
                        "option. Default path is: ./table.biom")
    parser.add_argument('--otu-fp', dest="otu_fp", metavar="OTU-FILE",
                        help="Create a file containing just (NCBI) OTU IDs "
                        "for use with a service such as phyloT "
                        "(http://phylot.biobyte.de/) to generate phylogenetic "
                        "trees for use in downstream analysis such as "
                        "UniFrac, iTol (itol.embl.de), or PhyloToAST "
                        "(phylotoast.org).")
    parser.add_argument('--fmt', default="hdf5", 
                        choices=["hdf5", "json", "tsv"],
                        help="Set the output format of the BIOM table. "
                              "Default is HDF5.")
    parser.add_argument('--store-pct', dest="store_pct", action='store_true',
                        help="Record the relative abundances "
                             "('Proportion_Classified' column) instead of "
                             "the raw count ('Count' column) data.")
    parser.add_argument('--gzip', action='store_true',
                        help="Compress the output BIOM table with gzip. "
                              "HDF5 BIOM (v2.x) files are internally "
                              "compressed by default, so this option "
                              "is not needed when specifying --fmt hdf5.")


    parser.add_argument('--version', action='version',                    
             version="clark-biom version {}, {}".format(__version__, __url__))
    parser.add_argument('-v', '--verbose', action='store_true',
                        help="Prints status messages during program "
                             "execution.")

    return parser.parse_args()


def main():
    args = handle_program_options()

    if args.fmt == 'hdf5' and not HAVE_H5PY:
        args.fmt = 'json'
        msg = """\
        Library 'h5py' not found, unable to write BIOM 2.x (HDF5) files.
        Defaulting to BIOM 1.0 (JSON)."""
        print(twdd(msg))

    # load all abundance table files and parse them
    sample_counts, taxa = process_samples(args.clark_abd_tbls, 
                                          store_pct=args.store_pct)

    # create new BIOM table from sample counts and taxon ids
    # add taxonomy strings to row (taxon) metadata
    biomT = create_biom_table(sample_counts, taxa)

    out_fp = write_biom(biomT, args.output_fp, args.fmt, args.gzip)

    if args.otu_fp:
        try:
            write_otu_file(list(taxa), args.otu_fp)
        except RuntimeError as re:
            msg = "ERROR creating OTU file: \n\t{}"
            sys.exit(msg.format(re))

    if args.verbose:
        print("".format(out_fp))
        table_str = """\
        BIOM-format table written to: {out_fp}
        Table contains {rows} rows (OTUs) and {cols} columns (Samples)
        and is {density:.1%} dense.""".format(out_fp=out_fp, 
                                              rows=biomT.shape[0], 
                                              cols=biomT.shape[1],
                                              density=biomT.get_table_density())
        print(twdd(table_str))


if __name__ == '__main__':
    main()

