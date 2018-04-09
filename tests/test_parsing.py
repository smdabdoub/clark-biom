#!/usr/bin/env python
# coding: utf-8
import csv
import importlib
import io
import os
import tempfile
from textwrap import dedent as twdd
import unittest

import clark_biom as cb


def prep_clark_input(data):
    cdr = csv.DictReader(data, fieldnames=cb.field_names, delimiter=",")
    return [entry for entry in cdr][1:]


class clark_biom_Test(unittest.TestCase):
    def setUp(self):
        self.sample_clark_rep = prep_clark_input(io.StringIO(twdd(u"""\
            Name,TaxID,Lineage,Count,Proportion_All(%),Proportion_Classified(%)
            Achromobacter xylosoxidans,85698,Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Alcaligenaceae;Achromobacter,82,0.00142317,0.124620061
            Acinetobacter baumannii,470,Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter,356,0.00617862,0.541033435
            Actinomyces cardiffensis,181487,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,1,1.74E-05,0.001519757
            Actinomyces dentalis,272548,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,15,0.000260335,0.022796353
            Actinomyces georgiae,52768,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,5,8.68E-05,0.007598784
            Actinomyces gerencseriae,52769,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,12,0.000208268,0.018237082
            Actinomyces israelii,1659,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,93,0.00161408,0.141337386
            Actinomyces johnsonii,544581,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,1,1.74E-05,0.001519757
            Actinomyces massiliensis,461393,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,12,0.000208268,0.018237082
            Actinomyces meyeri,52773,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,81,0.00140581,0.123100304
            UNKNOWN,UNKNOWN,UNKNOWN,658,92.0161,-
            """)))

    def run_parse_clark_report(self, manual):
        counts, _ = cb.parse_clark_abundance_tbl(self.sample_clark_rep)

        # Check that there are no differences in taxonomy IDs
        self.assertTrue(len(set(counts).symmetric_difference(set(manual))) == 0)

        # Check that the counts are equivalent
        self.assertTrue(all([manual[tax_id] == counts[tax_id] for tax_id in counts]))


    def test_parse_clark_report_S(self):
        manual = {'85698': 82, '470': 356, '181487': 1, '272548': 15,
                   '52768': 5, '52769': 12, '1659': 93, '544581': 1,
                   '461393': 12, '52773': 81}

        self.run_parse_clark_report(manual)


    def test_parse_clark_report_taxonomy(self):
        manual = {"85698": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", 
                            "o__Burkholderiales", "f__Alcaligenaceae", "g__Achromobacter", 
                            "s__xylosoxidans"],
                  "470": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria",
                          "o__Pseudomonadales", "f__Moraxellaceae", "g__Acinetobacter", 
                          "s__baumannii"],
                  "181487": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                             "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                             "s__cardiffensis"],
                  "272548": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                             "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                             "s__dentalis"],
                  "52768": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                            "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                            "s__georgiae"],
                  "52769": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                            "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                            "s__gerencseriae"],
                  "1659": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                           "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                           "s__israelii"],
                  "544581": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                             "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                             "s__johnsonii"],
                  "461393": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                             "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                             "s__massiliensis"],
                  "52773": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria",
                            "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces",
                            "s__meyeri"]
                 }

        _, taxa = cb.parse_clark_abundance_tbl(self.sample_clark_rep)

        # Check that the manually assigned taxa are eqivalent to the parsed
        self.assertTrue(all([manual[tax_id] == taxa[tax_id] for tax_id in manual]))


    def tearDown(self):
        pass



if __name__ == '__main__':
    unittest.main()