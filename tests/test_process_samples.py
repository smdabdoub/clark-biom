#!/usr/bin/env python
# coding: utf-8
import os, os.path as osp
import tempfile
from textwrap import dedent as twdd
import unittest

import clark_biom as cb
from tests.test_parsing import prep_clark_input



class clark_biom_Test(unittest.TestCase):
    def setUp(self):
        self.crepA = twdd(u"""\
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
            """).encode("utf-8")

        self.crepB = twdd(u"""\
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
            """).encode("utf-8")

        # create temp files containing the above clark results
        tempf_crepA = tempfile.NamedTemporaryFile(delete=False)
        tempf_crepA.write(self.crepA)
        tempf_crepA.close()

        tempf_crepB = tempfile.NamedTemporaryFile(delete=False)
        tempf_crepB.write(self.crepB)
        tempf_crepB.close()

        self.fps = [tempf_crepA.name, tempf_crepB.name]
        self.fnames = [osp.split(fp)[1] for fp in self.fps]

        self.sample_counts, self.taxa = cb.process_samples(self.fps, 
                                                           store_pct=False)
        #TODO: test relative abundance parsing
        # self.sample_counts_pct, self.taxa_pct = cb.process_samples(self.fps, 
        #                                                         store_pct=True)

    def test_sample_ids(self):
        proc_snames = set(self.sample_counts)

        self.assertEqual(proc_snames.symmetric_difference(self.fnames), set())


    def test_sample_counts(self):
        countsA = {'85698': 82, '470': 356, '181487': 1, '272548': 15,
                   '52768': 5, '52769': 12, '1659': 93, '544581': 1,
                   '461393': 12, '52773': 81}

        countsB = {'85698': 10, '470': 200, '1656': 5, '714': 212, '732': 2630, '739': 1,
                   '1176649': 1, '358': 9, '419015': 1, '33034': 3, '1658671': 5,
                   '1382': 7, '1383': 1}

        proc_countsA = self.sample_counts[self.fnames[0]]
        proc_countsB = self.sample_counts[self.fnames[1]]


        self.assertTrue(all([countsA[otu_id] == proc_countsA[otu_id]
                               for otu_id in countsA]))
        self.assertTrue(all([countsA[otu_id] == proc_countsA[otu_id]
                               for otu_id in proc_countsA]))

        self.assertTrue(all([countsB[otu_id] == proc_countsB[otu_id]
                               for otu_id in countsB]))
        self.assertTrue(all([countsB[otu_id] == proc_countsB[otu_id]
                               for otu_id in proc_countsB]))


    def test_taxa(self):
        manual = {"85698": ["k__Bacteria", "p__Proteobacteria", "c__Betaproteobacteria", "o__Burkholderiales", "f__Alcaligenaceae", "g__Achromobacter", "s__xylosoxidans"],
                  "470": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pseudomonadales", "f__Moraxellaceae", "g__Acinetobacter", "s__baumannii"],
                  "181487": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__cardiffensis"],
                  "272548": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__dentalis"],
                  "52768": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__georgiae"],
                  "52769": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__gerencseriae"],
                  "1659": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__israelii"],
                  "544581": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__johnsonii"],
                  "461393": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__massiliensis"],
                  "52773": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__meyeri"],
                  "1656": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Actinomycetales", "f__Actinomycetaceae", "g__Actinomyces", "s__viscosus"],
                  "714": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Aggregatibacter", "s__actinomycetemcomitans"],
                  "732": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Aggregatibacter", "s__aphrophilus"],
                  "739": ["k__Bacteria", "p__Proteobacteria", "c__Gammaproteobacteria", "o__Pasteurellales", "f__Pasteurellaceae", "g__Aggregatibacter", "s__segnis"],
                  "1176649": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Rhizobiaceae", "g__Agrobacterium", "s__fabrum"],
                  "358": ["k__Bacteria", "p__Proteobacteria", "c__Alphaproteobacteria", "o__Rhizobiales", "f__Rhizobiaceae", "g__Agrobacterium", "s__tumefaciens"],
                  "419015": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Bifidobacteriales", "f__Bifidobacteriaceae", "g__Alloscardovia", "s__omnicolens"],
                  "33034": ["k__Bacteria", "p__Firmicutes", "c__Tissierellia", "o__Tissierellales", "f__Peptoniphilaceae", "g__Anaerococcus", "s__prevotii"],
                  "1658671": ["k__Bacteria", "p__Actinobacteria", "c__Actinobacteria", "o__Micrococcales", "f__Intrasporangiaceae", "g__Arsenicicoccus", "s__sp. oral taxon 190"],
                  "1382": ["k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia", "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium", "s__parvulum"],
                  "1383": ["k__Bacteria", "p__Actinobacteria", "c__Coriobacteriia", "o__Coriobacteriales", "f__Atopobiaceae", "g__Atopobium", "s__rimae"]
                 }

        self.assertTrue(all([manual[tax_id] == self.taxa[tax_id] 
                               for tax_id in manual]))



    def tearDown(self):
        for fp in self.fps:
          os.unlink(fp)
          assert not os.path.exists(fp)



if __name__ == '__main__':
    unittest.main()
            