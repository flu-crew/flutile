#!/usr/bin/env python3

import flutile.functions as f
import unittest
import datetime
import smof


class TestIndexedAADiff(unittest.TestCase):
    def test_gapped_indices(self):
        self.assertEqual(f.gapped_indices(""), [])
        self.assertEqual(f.gapped_indices("ATG"), ["1", "2", "3"])
        self.assertEqual(f.gapped_indices("--ATG"), ["-2", "-1", "1", "2", "3"])
        self.assertEqual(
            f.gapped_indices("--ATG--"), ["-2", "-1", "1", "2", "3", "3+1", "3+2"]
        )
        self.assertEqual(
            f.gapped_indices("--A--T-G--"),
            ["-2", "-1", "1", "1+1", "1+2", "2", "2+1", "3", "3+1", "3+2"],
        )
        self.assertEqual(f.gapped_indices("-"), ["-1"])
        self.assertEqual(f.gapped_indices("---"), ["-3", "-2", "-1"])

    #  def test_indexed_aa_diff_table(self):
    #      #  --GA--T--   index sequence
    #      #  AAT-GGTTT   reference sequence
    #      #  AATTGGATT   a compared sequence
    #      #  -------------------------------
    #      #  Site   Ref   S1
    #      #  -2     A
    #      #  -1     A
    #      #  1      T
    #      #  2      -     T
    #      #  2+1    G     T
    #      #  2+2    G
    #      #  3      T     G
    #      #  3+1    T
    #      #  3+2    T
    #      simple = [("Ind", "--GA--T--"), ("Ref", "AAT-GGTTT"),  ("S1", "AATTGGATT")]
    #      self.assertEqual(f.inexed_aa_diff_table(simple),
    #          (   ["Site", "Ref", "S1"],
    #            [
    #              [ "-2" , "A", ""  ],
    #              [ "-1" , "A", ""  ],
    #              [ "1"  , "T", ""  ],
    #              [ "2"  , "-", "T" ],
    #              [ "2+1", "G", "T" ],
    #              [ "2+2", "G", ""  ],
    #              [ "3"  , "T", "A" ],
    #              [ "3+1", "T", ""  ]
    #              [ "3+2", "T", ""  ]
    #            ]
    #          )


class TestParsers(unittest.TestCase):
    def test_aadiff(self):
        simple_seqs = [("A", "AAAAA"), ("B", "TAAAA"), ("C", "TAAAA"), ("D", "TAAAC")]
        self.assertEqual(
            list(f.aadiff_table(simple_seqs)),
            [
                ["site", "A", "B", "C", "D"],
                ["1", "A", "T", "T", "T"],
                ["5", "A", "", "", "C"],
            ],
        )

    def test_pident(self):
        self.assertEqual(f.pident("AAAAAAAAAA", "AAAAAAAAAA"), 1.0)
        self.assertEqual(f.pident("AAAAAAAAAA", "AAAAAAAAAT"), 0.9)
        self.assertEqual(f.pident("AAAAAA--AAAA", "AAAAAG--AAAT"), 0.8)
        # No ungapped aligned region
        self.assertEqual(f.pident("----AAAA", "AAAA----"), 0)
        # An error is raised if the sequences are not of equal length
        self.assertRaises(f.InputError, f.pident, "AAAAAAAAAA", "AAA")

    def test_parseOutDate(self):
        self.assertEqual(
            f.parseOutDate("ladida 2019-01-01 fodidu"), datetime.date(2019, 1, 1)
        )

        # extracts the FIRST date that is observed (is this really what I want?)
        self.assertEqual(
            f.parseOutDate("ladida 2019-01-01 fodidu 2018-02-02"),
            datetime.date(2019, 1, 1),
        )

    def test_components(self):
        self.assertEqual(f.components([]), [])
        self.assertEqual(f.components([(1, 2)]), [{1, 2}])
        self.assertEqual(f.components([(1, 2), (2, 3), (4, 5)]), [{1, 2, 3}, {4, 5}])
        self.assertEqual(
            f.components([(1, 2), (2, 3), (4, 5), (1, 5)]), [{1, 2, 3, 4, 5}]
        )
        self.assertEqual(
            f.components([(1, 2), (2, 3), (4, 5), (5, 1)]), [{1, 2, 3, 4, 5}]
        )
        self.assertEqual(
            f.components([(1, 2), (6, 7), (2, 3), (4, 5), (5, 1), (7, 8)]),
            [{1, 2, 3, 4, 5}, {6, 7, 8}],
        )

    def test_represent(self, maxDiff=300):
        fasta = [
            ("A|2019-06-01", "AAAAA"),
            ("B|2019-06-06", "AAAAA"),
            ("C|2019-06-06", "TTAAA"),
            ("D|2018-06-06", "AAAAA"),
        ]
        g, a = f.represent(fasta, max_day_sep=5, min_pident_sep=1.0, same_state=False)
        b = fasta[1:]
        self.assertEqual(sorted(a), b)

    def test_represent_states(self, maxDiff=300):
        fasta = [
            ("A|Iowa|2019-06-01", "AAAAA"),
            ("B|Iowa|2019-06-06", "AAAAA"),
            ("C|Iowa|2019-06-06", "TTAAA"),
            ("D|Iowa|2018-06-06", "AAAAA"),
        ]
        g, a = f.represent(fasta, max_day_sep=5, min_pident_sep=1.0, same_state=True)
        b = fasta[1:]
        self.assertEqual(sorted(a), b)

        fasta = [
            ("A|Iowa|2019-06-01", "AAAAA"),
            ("B|Nebraska|2019-06-06", "AAAAA"),
            ("C|Iowa|2019-06-06", "TTAAA"),
            ("D|Iowa|2018-06-06", "AAAAA"),
        ]
        g, a = f.represent(fasta, max_day_sep=5, min_pident_sep=1.0, same_state=True)
        self.assertEqual(sorted(a), fasta)

    def test_ha_range_map(self):
        # H1   | H2  | H3  | H4  | H5  | H6  | H7  | H8  | H9  | H10 | H11 | H12 | H13 | H14 | H15 | H16 | H17 | H18 |
        # ---  | --  | --  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
        #      |     | 1   |     |     |     |     |     |     |     |     |     |     |     |     |     |     |     |
        #      |     | 2   |     |     |     |     |     |     |     |     |     |     | 1   |     |     |     |     |
        #      |     | 3   |     |     |     |     |     |     |     |     |     |     | 2   |     |     |     |     |
        #      |     | 4   |     |     |     |     |     |     |     |     |     |     | 3   |     |     |     |     |
        #      |     | 5   | 1   |     |     |     |     |     |     |     |     |     | 4   |     |     |     |     |
        #      |     | 6   | 2   |     |     |     |     |     |     |     |     |     | 5   |     |     |     |     |
        #      |     | 7   | 3   |     |     |     |     |     |     |     |     |     | 6   |     |     |     |     |
        #      |     | 8   | 4   |     |     |     |     |     |     |     |     |     | 7   |     |     |     |     |
        #      |     | 9   | 5   |     |     |     |     |     |     |     |     |     | 8   |     |     |     |     |
        #      |     | 10  | 6   |     |     |     | 1   |     | 1   |     | 1   |     | 9   |     |     |     |     |
        #  1   | 1   | 11  | 7   |  1  | 1   | 1   | 2   | 1   | 2   | 1   | 2   | 1   | 10  | 1   | 1   | 1   |  1  |
        #  2   | 2   | 12  | 8   |  2  | 2   | 2   | 3   | 2   | 3   | 2   | 3   | 2   | 11  | 2   | 2   | 2   |  2  |
        #  3   | 3   | 13  | 9   |  3  | 3   | 3   | 4   | 3   | 4   | 3   | 4   | 3   | 12  | 3   | 3   | 3   |  3  |
        #  4   | 4   | 14  | 10  |  4  | 4   | 4   | 5   | 4   | 5   | 4   | 5   | 4   | 13  | 4   | 4   | 4   |  4  |
        #  5   | 5   | 15  | 11  |  5  | 5   | 5   | 6   | 5   | 6   | 5   | 6   | 5   | 14  | 5   | 5   | 5   |  5  |
        #  6   | 6   | 16  | 12  |  6  | 6   | 6   | 7   | 6   | 7   | 6   | 7   | 6   | 15  | 6   | 6   | 6   |  6  |
        #  7   | 7   | 17  | 13  |  7  | 7   | 7   | 8   | 7   | 8   | 7   | 8   | 7   | 16  | 7   | 7   | 7   |  7  |
        #  8   | 8   | 18  | 14  |  8  | 8   | 8   | 9   | 8   | 9   | 8   | 9   | 8   | 17  | 8   | 8   | 8   |  8  |
        #  ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
        #  507 | 505 | 509 | 507 | 510 | 509 | 503 | 508 | 500 | 504 | 506 | 506 | 505 | 510 | 511 | 505 | 504 |  -  |
        #  508 | 506 | 510 | 508 | 511 | 510 | 504 | 509 | 501 | 505 | 507 | 507 | 506 | 511 | 512 | 506 | 505 |  -  |
        #  509 | 507 | -   | -   | 512 | 511 | -   | 510 | 502 | -   | 508 | 508 | 507 | -   | -   | 507 | 506 |  -  |
        #  -   | -   | -   | -   | -   | -   | -   | -   | -   | -   | 509 | -   | 508 | -   | -   | 508 | -   |  -  |
        #  510 | 508 | 511 | 509 | 513 | 512 | 505 | 511 | 503 | 506 | 510 | 509 | 509 | 512 | 513 | 509 | 507 |  -  |
        #  511 | 509 | 512 | 510 | 514 | 513 | 506 | 512 | 504 | 507 | 511 | 510 | 510 | 513 | 514 | 510 | 508 |  -  |
        #  ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
        #  545 | 543 | 546 | 544 | 548 | 547 | 540 | 546 | -   | 541 | 545 | 544 | 544 | 547 | 548 | 544 | 542 |  -  |
        #  546 | 544 | 547 | 545 | 549 | 548 | 541 | 547 | -   | 542 | 546 | 545 | 545 | 548 | 549 | 545 | 543 |  -  |
        #  547 | 545 | 548 | 546 | 550 | 549 | 542 | 548 | -   | 543 | 547 | 546 | 546 | 549 | 550 | 546 | 544 |  -  |
        #  548 | 546 | 549 | 547 | 551 | 550 | 543 | 549 | -   | 544 | 548 | 547 | 547 | 550 | 551 | 547 | 545 |  -  |
        #  549 | 547 | 550 | 548 | 552 | 551 | 544 | 550 | -   | 545 | 549 | 548 | 548 | 551 | 552 | 548 | 546 |  -  |

        # trivial test that positions map to themselves
        self.assertEqual(f.map_ha_range(start=1, end=3, subtype1=1, subtype2=1), (1, 3))
        # x | x
        self.assertEqual(
            f.map_ha_range(start=512, end=512, subtype1=5, subtype2=6), (511, 511)
        )

        #  H3  | H4
        #  --  | ---
        #  1   |
        #  2   |
        #  3   |
        #  4   |
        #  5   | 1
        #  6   | 2
        self.assertEqual(
            f.map_ha_range(start=1, end=3, subtype1=3, subtype2=4), (None, None)
        )
        self.assertEqual(f.map_ha_range(start=3, end=5, subtype1=3, subtype2=4), (1, 1))
        self.assertEqual(f.map_ha_range(start=1, end=2, subtype1=4, subtype2=3), (5, 6))

        #  H12 | H13
        #  --- | ---
        #  508 | 507
        #  -   | 508
        #  509 | 509
        self.assertEqual(
            f.map_ha_range(start=508, end=509, subtype1=12, subtype2=13),
            (507, 509),
        )  # H12->H13
        self.assertEqual(
            f.map_ha_range(start=507, end=509, subtype1=13, subtype2=12),
            (508, 509),
        )  # H13->H12

        #  H1  | H2
        #  --- | ---
        #  508 | 506
        #  509 | 507
        #  -   | -
        self.assertEqual(
            f.map_ha_range(start=508, end=509, subtype1=1, subtype2=2), (506, 507)
        )

        #  H1  | H2
        #  --- | ---
        #  -   | -
        #  510 | 508
        #  511 | 509
        self.assertEqual(
            f.map_ha_range(start=510, end=511, subtype1=1, subtype2=2), (508, 509)
        )

        #  H10 | H11
        #  --- | ---
        #  -   | 509
        #  506 | 510
        #  507 | 511
        self.assertEqual(
            f.map_ha_range(start=506, end=507, subtype1=10, subtype2=11), (510, 511)
        )
        self.assertEqual(
            f.map_ha_range(start=510, end=511, subtype1=11, subtype2=10), (506, 507)
        )
        self.assertEqual(
            f.map_ha_range(start=509, end=511, subtype1=11, subtype2=10),
            (506, 507),
        )

        # H12 | H13
        # --- | ---
        # 507 | 506
        # 508 | 507
        # -   | 508
        self.assertEqual(
            f.map_ha_range(start=507, end=508, subtype1=12, subtype2=13), (506, 507)
        )
        self.assertEqual(
            f.map_ha_range(start=506, end=507, subtype1=13, subtype2=12), (507, 508)
        )
        self.assertEqual(
            f.map_ha_range(start=506, end=508, subtype1=13, subtype2=12),
            (507, 508),
        )

        # === End Studies ===
        # H1   | H2  | H3  | H4  | H5  | H6  | H7  | H8  | H9  | H10 | H11 | H12 | H13 | H14 | H15 | H16 | H17 | H18 |
        # ---  | --  | --  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
        #  549 | 547 | 550 | 548 | 552 | 551 | 544 | 550 | -   | 545 | 549 | 548 | 548 | 551 | 552 | 548 | 546 |  -  |
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=2), (547, 547)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=3), (550, 550)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=4), (548, 548)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=5), (552, 552)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=6), (551, 551)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=7), (544, 544)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=8), (550, 550)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=9), (None, None)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=10), (545, 545)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=11), (549, 549)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=12), (548, 548)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=13), (548, 548)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=14), (551, 551)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=15), (552, 552)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=16), (548, 548)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=17), (546, 546)
        )
        self.assertEqual(
            f.map_ha_range(start=549, end=549, subtype1=1, subtype2=18), (None, None)
        )

    def test_ungap_indices(self):
        self.assertEqual(f.ungap_indices(start=1, end=1, fasta="A"), (1, 1))
        self.assertEqual(f.ungap_indices(start=1, end=1, fasta="A-"), (1, 1))
        self.assertEqual(f.ungap_indices(start=1, end=1, fasta="-A"), (2, 2))
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta="GATACA-"), (2, 3))
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta="-GATACA"), (3, 4))
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta="G-ATACA"), (3, 4))
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta="GA-TACA"), (2, 4))
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta="GAT-ACA"), (2, 3))
        # alternative gap characters work
        self.assertEqual(f.ungap_indices(start=2, end=3, fasta=".GA_TACA"), (3, 5))
        # end studies
        self.assertEqual(f.ungap_indices(start=4, end=6, fasta="-GA-TACA"), (6, 8))
        self.assertEqual(f.ungap_indices(start=4, end=600, fasta="-GA-TACA"), (6, 8))
        self.assertEqual(f.ungap_indices(start=500, end=600, fasta="-GA-TACA"), (9, 9))

    def test_parse_motif(self):
        self.assertEqual(
            f.parse_motif(motif_str=" 1 ", subtype="H1"), ("H1:1", [(1, 1)])
        )
        self.assertEqual(
            f.parse_motif(motif_str="x = 1 ", subtype="H1"), ("x", [(1, 1)])
        )
        self.assertEqual(
            f.parse_motif(motif_str="x= 1 , 3", subtype="H1"),
            ("x", [(1, 1), (3, 3)]),
        )
        self.assertEqual(
            f.parse_motif(motif_str=" x=1, 3 -5", subtype="H1"),
            ("x", [(1, 1), (3, 5)]),
        )

    def test_concat(self):
        self.assertEqual(f.concat([]), [])
        self.assertEqual(f.concat([[1]]), [1])
        self.assertEqual(f.concat([[1, 2], [3, 4, 5], []]), [1, 2, 3, 4, 5])

    def test_unconcat(self):
        self.assertEqual(f.unconcat([], [1, 1]), [])
        self.assertEqual(f.unconcat(["1234", "567"], [1, 1]), ["1234", "567"])
        self.assertEqual(f.unconcat(["1234", "567"], [2]), ["1234567"])
        self.assertEqual(f.unconcat(["1", "234", "567"], [2, 1]), ["1234", "567"])
        # it isn't entirely obvious what to do if the widths don't sum to the expected size
        # perhaps this should raise an error?
        self.assertEqual(f.unconcat(["1", "234", "567"], [2]), ["1234"])
        # and if there are extra width terms, I just ignore them
        self.assertEqual(f.unconcat(["1", "234", "567"], [1, 2, 10]), ["1", "234567"])
        self.assertEqual(f.unconcat(["1", "234", "567"], [1, 10]), ["1", "234567"])
        # 0's do nothing
        self.assertEqual(
            f.unconcat(["1", "234", "567"], [0, 1, 10, 0]), ["1", "234567"]
        )
        # negative values do nothing
        self.assertEqual(
            f.unconcat(["1", "234", "567"], [-42, 1, 10, 0]), ["1", "234567"]
        )

    def test_map_dna2dna(self):
        fna1 = [smof.FastaEntry(header="A", seq="ATGTTTAAATTTAAA")]
        aln1 = [
            smof.FastaEntry(header="ref", seq="MGYGY"),
            smof.FastaEntry(header="A", seq="MFKFK"),
        ]

        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 1)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "ATG")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 2)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "ATGTTT")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 5)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "ATGTTTAAATTTAAA")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(2, 5)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "TTTAAATTTAAA")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 6)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "ATGTTTAAATTTAAA")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 60)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "ATGTTTAAATTTAAA")],
        )
        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(50, 60)], fna=fna1, aln=aln1)[0]
            ],
            [("A", "")],
        )

        fna2 = [
            smof.FastaEntry(header="A", seq="ATGTTTAAATTTAAA"),
            smof.FastaEntry(header="B", seq="ATGTTTAAATTTAAA"),
        ]
        aln2 = [
            smof.FastaEntry(header="ref", seq="MGYGY"),
            smof.FastaEntry(header="A", seq="MFKFK"),
            smof.FastaEntry(header="B", seq="MFKFK"),
        ]

        self.assertEqual(
            [
                smof.to_pair(x)
                for x in f.map_dna2dna(bounds=[(1, 1)], fna=fna2, aln=aln2)[0]
            ],
            [("A", "ATG"), ("B", "ATG")],
        )
        self.assertEqual(
            [
                [smof.to_pair(x) for x in xs]
                for xs in f.map_dna2dna(bounds=[(1, 1), (2, 2)], fna=fna2, aln=aln2)
            ],
            [[("A", "ATG"), ("B", "ATG")], [("A", "TTT"), ("B", "TTT")]],
        )


if __name__ == "__main__":
    unittest.main()
