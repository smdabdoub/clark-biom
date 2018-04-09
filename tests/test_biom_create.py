#!/usr/bin/env python
# coding: utf-8
from collections import OrderedDict
import csv
import importlib
import io
import os
import tempfile
from textwrap import dedent as twdd
import unittest

import clark_biom as cb
from tests.test_parsing import prep_clark_input



class clark_biom_Test(unittest.TestCase):
    def setUp(self):
        self.krepA = prep_clark_input(io.StringIO(twdd(u"""\
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
        
        self.krepB = prep_clark_input(io.StringIO(twdd(u"""\
            Name,TaxID,Lineage,Count,Proportion_All(%),Proportion_Classified(%)
            Achromobacter xylosoxidans,85698,Bacteria;Proteobacteria;Betaproteobacteria;Burkholderiales;Alcaligenaceae;Achromobacter,10,0.00142317,0.003241491
            Acinetobacter baumannii,470,Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Moraxellaceae;Acinetobacter,200,0.00617862,0.064829822
            Actinomyces viscosus,1656,Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Actinomycetaceae;Actinomyces,5,8.68E-05,0.001620746
            Aggregatibacter actinomycetemcomitans,714,Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae;Aggregatibacter,212,0.0036794,0.068719611
            Aggregatibacter aphrophilus,732,Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae;Aggregatibacter,2630,0.0456454,0.852512156
            Aggregatibacter segnis,739,Bacteria;Proteobacteria;Gammaproteobacteria;Pasteurellales;Pasteurellaceae;Aggregatibacter,1,1.74E-05,0.000324149
            Agrobacterium fabrum,1176649,Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Agrobacterium,1,1.74E-05,0.000324149
            Agrobacterium tumefaciens,358,Bacteria;Proteobacteria;Alphaproteobacteria;Rhizobiales;Rhizobiaceae;Agrobacterium,9,0.000156201,0.002917342
            Alloscardovia omnicolens,419015,Bacteria;Actinobacteria;Actinobacteria;Bifidobacteriales;Bifidobacteriaceae;Alloscardovia,1,1.74E-05,0.000324149
            Anaerococcus prevotii,33034,Bacteria;Firmicutes;Tissierellia;Tissierellales;Peptoniphilaceae;Anaerococcus,3,5.21E-05,0.000972447
            Arsenicicoccus sp. oral taxon 190,1658671,Bacteria;Actinobacteria;Actinobacteria;Micrococcales;Intrasporangiaceae;Arsenicicoccus,5,8.68E-05,0.001620746
            Atopobium parvulum,1382,Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Atopobiaceae;Atopobium,7,0.00012149,0.002269044
            Atopobium rimae,1383,Bacteria;Actinobacteria;Coriobacteriia;Coriobacteriales;Atopobiaceae;Atopobium,1,1.74E-05,0.000324149
            UNKNOWN,UNKNOWN,UNKNOWN,3085,92.0161,-
            """)))

        # parse the sample reports
        self.taxa = OrderedDict()
        self.sample_counts = OrderedDict()

        countsA, taxaA = cb.parse_clark_abundance_tbl(self.krepA)
        countsB, taxaB = cb.parse_clark_abundance_tbl(self.krepB)
        self.taxa.update(taxaA)
        self.taxa.update(taxaB)
        self.sample_counts["A"] = countsA
        self.sample_counts["B"] = countsB

        # create the BIOM table from the sample counts and taxa
        self.biomT = cb.create_biom_table(self.sample_counts, self.taxa)

    def test_sample_ids(self):
        biom_sample_ids = self.biomT.ids(axis="sample")

        self.assertEqual(len(biom_sample_ids), 2)
        self.assertEqual(set(biom_sample_ids).symmetric_difference(["A","B"]),
                         set())

    def test_observation_ids(self):
        biom_otu_ids = self.biomT.ids(axis="observation")
        otu_ids = {'1176649', '1382', '1383', '1656', '1658671', '1659',
                   '181487', '272548', '33034', '358', '419015', '461393',
                   '470', '52768', '52769', '52773', '544581', '714',
                   '732', '739', '85698'}

        self.assertEqual(len(biom_otu_ids), len(otu_ids))
        self.assertEqual(set(biom_otu_ids).symmetric_difference(otu_ids), 
                         set())

    def test_sample_counts(self):
        countsA = {'85698': 82, '470': 356, '181487': 1, '272548': 15,
                   '52768': 5, '52769': 12, '1659': 93, '544581': 1,
                   '461393': 12, '52773': 81}

        countsB = {'85698': 10, '470': 200, '1656': 5, '714': 212, '732': 2630, '739': 1,
                   '1176649': 1, '358': 9, '419015': 1, '33034': 3, '1658671': 5,
                   '1382': 7, '1383': 1}

        self.assertTrue(all([self.biomT.get_value_by_ids(otu_id, "A") 
                             == countsA[otu_id] for otu_id in countsA]))

        self.assertTrue(all([self.biomT.get_value_by_ids(otu_id, "B") 
                             == countsB[otu_id] for otu_id in countsB]))



    def tearDown(self):
        pass



if __name__ == '__main__':
    unittest.main()
            