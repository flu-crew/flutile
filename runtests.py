#!/usr/bin/env python3

import flutile.main as f
import unittest


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


if __name__ == "__main__":
    unittest.main()
