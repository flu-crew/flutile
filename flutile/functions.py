import collections
import sys
import datetime
import re
import tqdm
import tempfile
import smof
import os
import flutile.motifs as motifs
from flutile.functions import *
from flutile.version import __version__


class InputError(Exception):
    pass


def is_aligned(fasta):
    lengths = {len(s.seq) for s in fasta}
    if len(lengths) != 1:
        err("Expected all files to be of equal length")


def err(msg):
    raise InputError("msg")


def with_aligned_pairs(alignment, f, *args, **kwargs):
    fasta_obj = smof.open_fasta(alignment)
    is_aligned(fasta_obj)
    s = [(s.header, s.seq) for s in fasta_obj]
    f(s, *args, **kwargs)


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

    print(e, file=sys.stderr)

    with open(".tmp_aln", "w+") as fh:
        print(o, file=fh)

    aln = list(smof.open_fasta(".tmp_aln"))

    return aln


def is_aligned(fasta):
    lengths = {len(s.seq) for s in fasta}
    if len(lengths) != 1:
        print(lengths, file=sys.stderr)
        err("Expected all files to be of equal length")


def ungap_indices(start, end, fasta):
    """
    Map indices in an ungapped sequence to indices in the gapped string
    """
    original_index = 0
    (a, b) = (0, 0)
    for (i, c) in enumerate(fasta):
        if c not in "._-":
            original_index += 1
            if original_index == start:
                a = i + 1
            if original_index == end:
                b = i + 1
                break
    return (a, b)


def extract_aa2aa(start, end, faa, ref, mafft_exe="mafft"):
    # add reference sequences to the input seq
    faa = list(ref) + list(faa)

    # align the reference and input protein sequences
    aln = align(faa, mafft_exe=mafft_exe)

    # find reference gapped start and stop positions
    (gapped_start, gapped_end) = ungap_indices(start, end, list(aln)[0].seq)

    # extract the subset region
    extracted = smof.subseq(aln, gapped_start, gapped_end)

    # return sans reference
    return list(extracted)[1:]


def extract_dna2dna(start, end, fna, ref, mafft_exe="mafft"):
    # translate the DNA inputs (longest uninterupted CDS)
    faa = smof.translate(fna, all_frames=True)

    # because generators are the root of all evil
    ref = list(ref)
    faa = list(faa)
    fna = list(fna)
    aln = list(align(ref + faa, mafft_exe=mafft_exe))

    # find reference gapped start and stop positions
    (a, b) = ungap_indices(start, end, list(aln)[0].seq)

    # find (start, length) bounds for each DNA entry
    for i in range(len(fna)):
        # start is 0-based
        # whereas 'a' and 'b' above are 1-based
        start, length = smof.find_max_orf(fna[i].seq, from_start=False)

        seq = aln[i + 1].seq

        # motif start position, 0-based
        aa_offset = len(seq[0 : (a - 1)].replace("-", ""))
        aa_length = len(seq[(a - 1) : b].replace("-", ""))

        dna_start = start + aa_offset * 3
        dna_end = dna_start + aa_length * 3 + 1

        fna[i].seq = fna[i].seq[dna_start:dna_end]

    return fna


def _dispatch_extract(
    start, end, subtype, fasta_file, conversion=None, *args, **kwargs
):
    # get reference for this subtype
    ref_file = os.path.join(os.path.dirname(__file__), "data", "subtype-refs.faa")
    ref_fasta = smof.open_fasta(ref_file)
    ref = list(smof.grep(ref_fasta, pattern=f"|H{subtype}N"))

    # open input sequences
    entries = smof.open_fasta(fasta_file)
    # remove gaps
    entries = list(smof.clean(entries, toseq=True))
    # dispatch by sequence type
    if conversion == "aa2aa":
        f = extract_aa2aa
    elif conversion == "dna2aa":
        entries = smof.translate(entries, all_frames=True)
        f = extract_aa2aa
    elif conversion == "dna2dna":
        f = extract_dna2dna
    else:
        err("You shouldn't have done that")

    return f(start, end, entries, ref, *args, **kwargs)


def extract_ha1(subtype, *args, **kwargs):
    """
    Get the HA1 range relative to the subtype of interest by mapping the H1 HA1
    range (18,344) range to the subtype of interest
    """
    (start, end) = map_ha_range(start=18, end=344, subtype1=1, subtype2=subtype)
    out = _dispatch_extract(start=start, end=end, subtype=subtype, *args, **kwargs)
    smof.print_fasta(out)


def gapped_indices(seq):
    def _handle_run(run, ind):
        if ind == 0:
            return list(reversed([("-" + str(i)) for i in range(1, run + 1)]))
        else:
            return [f"{str(ind)}+{str(i)}" for i in range(1, run + 1)]

    indices = []
    run = 0
    ind = 0
    for i in range(len(seq)):
        if seq[i] == "-":
            run += 1
        else:
            if run > 0:
                indices += _handle_run(run, ind)
                run = 0
            ind += 1
            indices.append(str(ind))
    if run > 0:
        indices += _handle_run(run, ind)
    return indices


def aadiff_table(
    s,
    make_consensus=False,
    consensus_as_reference=False,
    remove_unchanged=True,
    **kwargs,
):
    """
    From the input alignment s, make an aadiff table
    """

    def find_consensus(i):
        counts = collections.Counter(s[j][1][i] for j in range(1, len(s)))
        return counts.most_common()[0][0]

    # add the consensus header column, if needed
    if make_consensus or consensus_as_reference:
        consensus_seq = "".join(find_consensus(i) for i in range(len(s[0][1])))
        # the consensus column goes first if it is being used as a reference
        if consensus_as_reference:
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
            if not remove_unchanged or s[j][1][i] != ref:
                row = [position, ref]
                for k in seq_ids[1:]:
                    aa = s[k][1][i]
                    if aa == ref:
                        row.append("")
                    else:
                        row.append(aa)
                yield row
                break


def annotate_table(
    table,
    subtype=None,
    annotation_tables="",
    join_annotations=False,
    keep_signal=False,
    caton82=False,
    wiley81=False,
    **kwargs,
):
    def load_builtin_annotations(name, subtype_exp):
        if subtype_exp != subtype:
            err(f"{name} annotation is only defined for {subtype_exp}")

        if keep_signal:
            err("{name} annotation is incompatible with --keep_signal")

        return os.path.join(os.path.dirname(__file__), "data", f"{name}.txt")

    all_anntables = []

    if caton82:
        all_anntables.append(load_builtin_annotations("caton82", "H1"))

    if wiley81:
        all_anntables.append(load_builtin_annotations("wiley81", "H3"))

    annotations = {x[0]: [] for x in table[1:]}

    new_column_names = []

    if annotation_tables:
        all_anntables += [f.strip() for f in annotation_tables.split(",")]

    for annotation_file in all_anntables:

        with open(annotation_file, "r") as f:
            lines = [line.split("\t") for line in f.readlines()]
            new_annotations = {line[0]: [x.strip() for x in line[1:]] for line in lines}

            new_column_names += lines[0][1:]
            for k, vs in annotations.items():
                if k in new_annotations:
                    annotations[k] += new_annotations[k]
                else:
                    annotations[k] += [""] * len(lines[0][1:])

    if join_annotations:
        annotations = {
            k: [", ".join([v for v in vs if v])] for (k, vs) in annotations.items()
        }
        new_column_names = ["annotations"]

    for i, row in enumerate(table):
        if i == 0:
            table[0] += new_column_names
        else:
            table[i] += annotations[row[0]]

    return table


def referenced_table(
    faa,
    subtype=None,
    mafft_exe="mafft",
    keep_signal=False,
    remove_unchanged=True,
    **kwargs,
):
    if subtype:
        ref = motifs.get_ha_subtype_nterm_motif(subtype)
    else:
        ref = None

    if keep_signal:
        # each of the reference sequences starts at the signal peptide's
        # initial methionine, so to keep the signal we trim nothing
        trim = 0
    elif subtype:
        # the motif here is an exact match to the reference signal peptide, we
        # want to remove the signal peptide, so the trim length is simply the
        # peptide length
        trim = len(ref.motif)

    # open input sequences
    entries = smof.open_fasta(faa)
    # remove gaps
    entries = list(smof.clean(entries, toseq=True))

    if subtype:
        # trim off signal peptides or other gunk at the beginning of the reference
        seq = ref.sequence[trim:]

        refseq = smof.FastaEntry(header=ref.defline, seq=seq)

        entries = [refseq] + entries

    # align the reference and input protein sequences
    aln = align(entries, mafft_exe=mafft_exe)
    aln = [(s.header, s.seq) for s in aln]

    indices = gapped_indices(aln[0][1])

    if subtype:
        # remove the reference sequence
        aln = aln[1:]

    table = list(aadiff_table(aln, remove_unchanged=remove_unchanged, **kwargs))

    # set relative indices
    for i, row in enumerate(table):
        if i > 0:
            table[i][0] = indices[int(row[0]) - 1]

    return table


def referenced_aadiff_table(faa, **kwargs):
    table = referenced_table(faa, remove_unchanged=True, **kwargs)
    return annotate_table(table, **kwargs)


def map_ha_range(start, end, subtype1, subtype2):
    """
    Map AA indices between subtypes using the Burke2014 index map

    @param start integer start position from initial methionine with respect to subtype1
    @param end integer end position from initial methionine with respect to subtype1
    @param subtype1 an integer (1-18 currently) for the HA subtype
    @param subtype2 an integer (1-18 currently) for the HA subtype
    @return index range with respect to subtype2
    """
    subtype_file = os.path.join(
        os.path.dirname(__file__), "data", "burke2014-index-map.txt"
    )
    indices = []
    with open(subtype_file, "r") as f:
        i = 0
        for line in f.readlines():
            row = [x.strip() for x in line.split("\t")]
            try:
                i = int(row[subtype1 - 1])
            except:
                pass
            if i >= start:
                try:
                    indices.append(int(row[subtype2 - 1]))
                except:
                    pass
            if i == end:
                break

    try:
        return (indices[0], indices[-1])
    except:
        return (None, None)


def referenced_annotation_table(faa, subtype, **kwargs):

    table = referenced_table(faa, subtype=subtype, remove_unchanged=False, **kwargs)

    subtype_file = os.path.join(
        os.path.dirname(__file__), "data", "burke2014-index-map.txt"
    )
    subtype_number = int(subtype[1:])
    subtype_annotation = dict()
    with open(subtype_file, "r") as f:
        for line in f.readlines():
            row = [x.strip() for x in line.split("\t")]
            try:
                i = int(row[subtype_number - 1])
                subtype_annotation[str(i)] = row
            except:
                pass
    table[0] += ["H" + str(i) for i in range(1, 19)]
    for (idx, row) in enumerate(table[1:]):
        try:
            table[idx + 1] += subtype_annotation[row[0]]
        except:
            table[idx + 1] += ["\t"] * 18

    return annotate_table(table, subtype=subtype, **kwargs)
