#!/usr/bin/env python3

import flutile.main as f
import unittest
import datetime as dt
import sys


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
        self.assertEqual(f.pident("AAAAAAAAAA", "AAAAAAAAAA"), 100)
        self.assertEqual(f.pident("AAAAAAAAAA", "AAAAAAAAAT"), 90)
        self.assertEqual(f.pident("AAAAAA--AAAA", "AAAAAG--AAAT"), 80)
        # No ungapped aligned region
        self.assertEqual(f.pident("----AAAA", "AAAA----"), 0)
        # An error is raised if the sequences are not of equal length
        self.assertRaises(f.InputError, f.pident, "AAAAAAAAAA", "AAA")

    def test_parseOutDate(self):
        self.assertEqual(
            f.parseOutDate("ladida 2019-01-01 fodidu"), dt.date(2019, 1, 1)
        )

        # extracts the FIRST date that is observed (is this really what I want?)
        self.assertEqual(
            f.parseOutDate("ladida 2019-01-01 fodidu 2018-02-02"), dt.date(2019, 1, 1)
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
        fasta = [("A|2019-06-01", "AAAAA"), ("B|2019-06-06", "AAAAA"), ("C|2019-06-06", "TTAAA"), ("D|2018-06-06", "AAAAA")]
        a = sorted(f.represent(fasta, max_day_sep=5, min_pident_sep=100))
        b = [("B|2019-06-06", "AAAAA"), ("C|2019-06-06", "TTAAA"), ("D|2018-06-06", "AAAAA")]
        self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
