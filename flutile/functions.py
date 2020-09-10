import collections
import sys
import datetime
import re
import tqdm
import tempfile
import smof
import os
from flutile.functions import *
from flutile.version import __version__


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
        counts = collections.Counter(s[j][1][i] for j in range(1, len(s)))
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
        return datetime.date(*[int(x) for x in match.group().split("-")])
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
    if N < 1:
        return 0
    else:
        return identities / N


USA_STATES = [
    "alaska",
    "alabama",
    "arkansas",
    "arizona",
    "california",
    "colorado",
    "connecticut",
    "district_of_columbia",
    "delaware",
    "florida",
    "georgia",
    "hawaii",
    "iowa",
    "idaho",
    "illinois",
    "indiana",
    "kansas",
    "kentucky",
    "louisiana",
    "massachusetts",
    "maryland",
    "maine",
    "michigan",
    "minnesota",
    "missouri",
    "mississippi",
    "montana",
    "north_carolina",
    "north_dakota",
    "nebraska",
    "new_hampshire",
    "new_jersey",
    "new_mexico",
    "nevada",
    "new_york",
    "ohio",
    "oklahoma",
    "oregon",
    "pennsylvania",
    "rhode_island",
    "south_carolina",
    "south_dakota",
    "tennessee",
    "texas",
    "utah",
    "virginia",
    "vermont",
    "washington",
    "wisconsin",
    "west_virginia",
    "wyoming",
]


def getUsaState(x):
    x = x.lower().replace(" ", "_")
    for state in USA_STATES:
        if state in x:
            return state
    return None


def represent(s, max_day_sep, min_pident_sep, same_state):
    if max_day_sep is not None:
        dates = [parseOutDate(header) for (header, seq) in s]
    if same_state:
        states = [getUsaState(header) for (header, seq) in s]
    N = len(s)
    pairs = []
    paired = set()
    for i in tqdm.tqdm(range(N - 1)):
        for j in range(i + 1, N):
            close_by_seq = pident(s[i][1], s[j][1]) >= min_pident_sep
            if max_day_sep is not None:
                close_by_time = abs((dates[i] - dates[j]).days) <= max_day_sep
            else:
                close_by_time = True
            if same_state:
                close_by_state = states[i] and states[j] and states[i] == states[j]
            else:
                close_by_state = True
            if close_by_time and close_by_seq and close_by_state:
                pairs.append((i, j))
                paired.update([i, j])

    seqs = set(i for i in range(N) if i not in paired)

    groups = components(pairs)

    grouped = set()
    for group in groups:
        grouped.update(group)
        if max_day_sep is not None:
            # if we are using time, then keep the most recent
            group = list(reversed(sorted(list(group), key=lambda i: dates[i])))
        else:
            # otherwise keep the first alphabetically
            group = sorted(list(group), key=lambda i: s[i])
        seqs.add(group[0])

    groups = groups + [{i} for i in range(N) if not i in grouped]

    return (groups, [s[i] for i in seqs])


def components(pairs):
    if len(pairs) == 0:
        return []

    groupmap = collections.defaultdict(set)

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


def align(seq, mafft_exe="mafft"):
    from Bio.Align.Applications import MafftCommandline

    with open(".tmp", "w+") as fasta_fh:
        smof.print_fasta(seq, out=fasta_fh)

    o, e = MafftCommandline(mafft_exe, input=".tmp")()

    aligned_fh = open(".tmp_aln", "w+")

    print(o, file=aligned_fh)

    return list(smof.open_fasta(".tmp_aln"))


def extract_motif(alignment, motif, ref=None):
    matches = smof.grep(
        alignment,
        pattern=motif,
        match_sequence=True,
        gapped=True,
        perl_regexp=True,
        gff=True,
    )

    if ref:
        print("Checking against a reference", file=sys.stderr)
        for row in matches:
            els = row.split("\t")
            if els[0] == ref:
                lower, upper = els[3:5]
                break
        else:
            print(
                "The reference does not contain the motif, this is a bug",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        print("No reference, using the most common locus", file=sys.stderr)
        bounds = [(a, b) for a, b in [x.split("\t")[3:5] for x in matches]]
        boundCount = collections.Counter(bounds)
        lower, upper = boundCount.most_common()[0][0]

    mot = smof.subseq(alignment, int(lower), int(upper))

    return mot


def extract_aa_motif_from_dna(fasta_file, motif, aa_ref_file=None, mafft_exe="mafft"):
    fna = list(smof.open_fasta(fasta_file))
    faa = smof.translate(fna, all_frames=True)

    # add reference sequences
    if aa_ref_file:
        print(f"Using ref '{aa_ref_file}'", file=sys.stderr) 
        aa_ref = list(smof.open_fasta(aa_ref_file))
        faa = list(aa_ref) + list(faa)
    else:
        print(f"No ref found", file=sys.stderr) 
        aa_ref = Fasta([])

    refs = [s.header for s in aa_ref]

    if len(refs) > 0:
        mot_ref = refs[0]
    else:
        mot_ref = None

    aln = align(faa, mafft_exe=mafft_exe)

    mot = extract_motif(aln, motif, ref=mot_ref)

    # remove any reference sequences
    for header in refs:
        mot = smof.grep(mot, pattern=header, invert_match=True)

    return mot


def extract_h1_ha1(fasta_file, motif="DT[LI]C.*QSR", mafft_exe="mafft", cds=False):
    h1_ref_file = os.path.join(os.path.dirname(__file__), "data", "h1-ha1-ref.faa")

    ha1 = extract_aa_motif_from_dna(
        fasta_file, motif=motif, aa_ref_file=h1_ref_file, mafft_exe=mafft_exe
    )

    smof.print_fasta(ha1)


def extract_h3_ha1(fasta_file, motif="QKL.*QTR", mafft_exe="mafft", cds=False):
    h3_ref_file = os.path.join(os.path.dirname(__file__), "data", "h3-ha1-ref.faa")

    ha1 = extract_aa_motif_from_dna(
        fasta_file, motif=motif, aa_ref_file=h3_ref_file, mafft_exe=mafft_exe
    )

    smof.print_fasta(ha1)
