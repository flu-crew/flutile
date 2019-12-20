"""
Flu utilities

Usage:
    flutile aadiff [--make-consensus] [--use-consensus-as-reference] [<filename>]
"""

import signal
import os
from docopt import docopt
from collections import Counter
from flutile.version import __version__
import sys


def parseFasta(filehandle):
    seqList = []
    header = None
    seq = None
    for line in filehandle:
        line = line.strip()
        if not line or line[0] == "#":
            continue
        if line[0] == ">":
            if header and seq:
                seqList.append((header, seq))
            header = line[1:]
            seq = ""
        else:
            seq += line
    seqList.append((header, seq))
    return seqList


def aadiff_table(s, consensus=False, consensus_as_ref=False):
    def find_consensus(i):
        counts = Counter(s[j][1][i] for j in range(1, len(s)))
        return counts.most_common()[0][0]

    # add the consensus header column, if needed
    if consensus or consensus_as_ref:
        consensus_seq = "".join(find_consensus(i) for i in range(len(s[0][1])))
        # the consensus column goes first if it is being used as a reference
        if consensus_as_ref:
            s = [("Consensus", consensus_seq)] + s
        # otherwise it goes last
        else:
            s = s + [("Consensus", consensus_seq)]

    header = [x[0] for x in s]
    yield (["site"] + header)

    seq_ids = list(range(len(s)))
    aa_ids = list(range(len(s[0][1])))

    # for each amino acid position in the alignment
    for i in aa_ids:
        ref = s[0][1][i]
        position = str(i + 1)
        # for each sequence in the alignment
        for j in seq_ids:
            if s[j][1][i] != ref:
                row = [position, ref]
                for k in seq_ids[1:]:
                    aa = s[k][1][i]
                    if aa == ref:
                        row.append("")
                    else:
                        row.append(aa)
                yield row
                break


def main():
    if os.name is "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    args = docopt(__doc__, version=f"flutile {__version__}")

    if args["<filename>"]:
        f = open(args["<filename>"])
    else:
        f = sys.stdin

    if args["aadiff"]:
        s = parseFasta(f)
        for row in aadiff_table(
            s,
            consensus=args["--make-consensus"],
            consensus_as_ref=args["--use-consensus-as-reference"],
        ):
            print("\t".join(row))


if __name__ == "__main__":
    main()
