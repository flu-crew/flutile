"""
Flu utilities

Usage:
    flutile compare [--make-consensus] [--use-consensus-as-reference] [<alignment>]
    flutile represent --max-day-sep=<days> --min-pident-sep=<pident> [<alignment>]

Options:
    --max-day-sep=<days>       Maximum number of days separating members of a group [default: 60]
    --min-pident-sep=<pident>  Minimum percent identity difference between members of a group [default: 100]
"""

import signal
import os
from docopt import docopt
from collections import Counter
from flutile.version import __version__
import sys
import datetime as dt
import re
from collections import defaultdict


class InputError(Exception):
    pass


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


def parseOutDate(s):
    datepat = re.compile("\d\d\d\d-\d\d-\d\d")
    match = re.search(datepat, s)
    if match:
        return dt.date(*[int(x) for x in match.group().split("-")])
    else:
        return None


def pident(s1, s2):
    # assert that the sequences are of equal length
    if len(s1) != len(s2):
        raise InputError("Cannot compare sequences of different length")
    # count the number of identities (not counting gaps)
    identities = 0
    N = 0
    for (x, y) in zip(s1, s2):
        if x != "-" and y != "-":
            N += 1
            identities += x == y
    return 100 * identities / N


def represent(s, max_day_sep, min_pident_sep):
    dates = [parseOutDate(header) for (header, seq) in s]
    N = len(s)
    pairs = []
    paired = set()
    for i in range(N - 1):
        for j in range(i + 1, N):
            close_by_time = abs((dates[i] - dates[j]).days) <= max_day_sep
            close_by_seq = pident(s[i][1], s[j][1]) >= min_pident_sep
            if close_by_time and close_by_seq:
                pairs.append((i, j))
                paired.update([i,j])

    seqs = set(i for i in range(N) if i not in paired)

    for group in components(pairs):
        group = sorted(list(group), key=lambda i: dates[i])
        seqs.add(group[-1])

    return [s[i] for i in seqs]


def components(pairs):
    if len(pairs) == 0:
        return []

    groupmap = defaultdict(set)

    for (x, y) in pairs:
        groupmap[x].add(y)
        groupmap[y].add(x)

    def group(x, xs, xss):
        for y in xss[x]:
            if not y in xs:
                xs.add(y)
                xs.update(group(y, xs, xss))
        return xs

    groups = []
    while groupmap:
        (a, b) = list(groupmap.items())[0]
        xs = group(a, set(), groupmap)
        groupmap = {k: v for (k, v) in groupmap.items() if not k in xs}
        groups.append(xs)

    return groups


def main():
    if os.name is "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)

    args = docopt(__doc__, version=f"flutile {__version__}")

    if args["<filename>"]:
        f = open(args["<filename>"])
    else:
        f = sys.stdin

    if args["compare"]:
        s = parseFasta(f)
        for row in aadiff_table(
            s,
            consensus=args["--make-consensus"],
            consensus_as_ref=args["--use-consensus-as-reference"],
        ):
            print("\t".join(row))

    if args["represent"]:
        s = parseFasta(f)
        try:
          pident = float(args["--min-pident-sep"])
          if not (0.0 <= pident <= 100.0):
            print("Expected pident between 0 and 100", file=sys.stderr)
            exit(1)
        except TypeError as e:
          print("Expected pident to be a float", file=sys.stderr)
          exit(1)

        for (header, seq) in represent(
            s,
            max_day_sep=int(args["--max-day-sep"]),
            min_pident_sep=pident,
        ):
            print(">" + header)
            print(seq)


if __name__ == "__main__":
    main()
