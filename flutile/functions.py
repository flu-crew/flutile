from __future__ import annotations
import collections
import sys
import datetime
import re
import tqdm  # type: ignore
import smof  # type: ignore
import os
import flutile.motifs as motifs
from flutile.parameters import AadiffOpts, SeqOpts, MafftOpts, Conversion
from typing import List, Tuple, Set, Dict, TextIO, Optional, TypeVar, Callable, Iterable

A = TypeVar("A")


class InputError(Exception):
    pass


def is_aligned(fasta: List[smof.FastaEntry]) -> None:
    """
    Die if the input is not aligned.

    This is not very polite.
    """
    lengths = {len(s.seq) for s in fasta}
    if len(lengths) != 1:
        raise InputError(
            "FASTA file is not aligned, all entries should be of equal length"
        )


def with_aligned_pairs(alignment: TextIO, f, *args, **kwargs):
    fasta_obj = smof.open_fasta(alignment)
    is_aligned(fasta_obj)
    s = [(s.header, s.seq) for s in fasta_obj]
    f(s, *args, **kwargs)


def parseOutDate(s: str) -> datetime.date:
    datepat = re.compile("\d\d\d\d-\d\d-\d\d")
    match = re.search(datepat, s)
    if match:
        return datetime.date(*[int(x) for x in match.group().split("-")])
    else:
        raise InputError(f"No date (xxxx-xx-xx) found in header '{s}'")


def pident(s1: str, s2: str) -> float:
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
        return 0.0
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


def getUsaState(x: str) -> Optional[str]:
    x = x.lower().replace(" ", "_")
    for state in USA_STATES:
        if state in x:
            return state
    return None


def represent(
    s: List[smof.FastaEntry], max_day_sep: int, min_pident_sep: float, same_state: bool
) -> Tuple[List[Set[int]], List[smof.FastaEntry]]:
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
            group_lst = list(reversed(sorted(list(group), key=lambda i: dates[i])))
        else:
            # otherwise keep the first alphabetically
            group_lst = sorted(list(group), key=lambda i: s[i])
        seqs.add(group_lst[0])

    groups = groups + [{i} for i in range(N) if i not in grouped]

    return (groups, [s[i] for i in seqs])


def components(pairs: List[Tuple[int, int]]) -> List[Set[int]]:
    if len(pairs) == 0:
        return []

    groupmap: collections.defaultdict[int, Set[int]]
    groupmap = collections.defaultdict(set)

    for (x, y) in pairs:
        groupmap[x].add(y)
        groupmap[y].add(x)

    def group(x, xs, xss):
        for y in xss[x]:
            if y not in xs:
                xs.add(y)
                xs.update(group(y, xs, xss))
        return xs

    shrinking_groupmap = dict(groupmap)
    groups = []
    while shrinking_groupmap:
        (a, b) = list(shrinking_groupmap.items())[0]
        xs = group(a, set(), shrinking_groupmap)
        shrinking_groupmap = {
            k: v for (k, v) in shrinking_groupmap.items() if k not in xs
        }
        groups.append(xs)

    return groups


def is_ha(subtype: str) -> bool:
    return bool(re.fullmatch("H\d+", subtype))


def is_na(subtype: str) -> bool:
    return bool(re.fullmatch("N\d+", subtype))


def get_ref(subtype: str) -> List[smof.FastaEntry]:
    synonyms = {
        "M": "M1",
        "MP": "M1",
        "NS": "NS1",
    }
    if subtype in synonyms:
        subtype = synonyms[subtype]

    if is_ha(subtype):
        try:
            i = int(subtype[1:])
        except:
            raise InputError("Expected HA, but H was not followed by an integer")
        ref_file = os.path.join(os.path.dirname(__file__), "data", "subtype-refs.faa")
        ref_fasta = smof.open_fasta(ref_file)
        ref = list(smof.grep(ref_fasta, no_color=True, pattern=f"|H{i}N"))
    elif is_na(subtype):
        try:
            i = int(subtype[1:])
        except:
            raise InputError("Expected NA, but N was not followed by an integer")
        ref_file = os.path.join(
            os.path.dirname(__file__), "data", "na-subtype-refs.faa"
        )
        ref_fasta = smof.open_fasta(ref_file)
        ref = list(
            smof.grep(
                ref_fasta, perl_regexp=True, no_color=True, pattern=f"\\|H[0-9]+N{i}\\|"
            )
        )
    elif subtype in ["PB2", "PB1", "PA", "NP", "M1", "NS1"]:
        ref_file = os.path.join(os.path.dirname(__file__), "data", "internal-refs.faa")
        ref_fasta = smof.open_fasta(ref_file)
        ref = list(
            smof.grep(
                ref_fasta, perl_regexp=True, no_color=True, pattern=f"\\|{subtype}$"
            )
        )
    else:
        raise InputError("Unexpected subtype")
    return ref


def align(
    seq: List[smof.FastaEntry], mafft_opts: MafftOpts = MafftOpts()
) -> List[smof.FastaEntry]:
    from Bio.Align.Applications import MafftCommandline  # type: ignore
    from tempfile import mkstemp

    hashy = smof.md5sum(seq)
    (n, fasta_filename) = mkstemp(suffix=f"-{hashy}.fa")

    with open(fasta_filename, "w+") as fasta_fh:
        smof.print_fasta(seq, out=fasta_fh)

    o, e = MafftCommandline(mafft_opts.mafft_exe, input=fasta_filename)()

    if mafft_opts.verbose:
        print(e, file=sys.stderr)

    (n, aln_filename) = mkstemp(suffix=f"-{hashy}.aln")

    with open(aln_filename, "w+") as fh:
        print(o, file=fh)

    aln = list(smof.open_fasta(aln_filename))

    return aln


def ungap_indices(
    start: int, end: int, fasta: List[smof.FastaEntry]
) -> Tuple[int, int]:
    """
    Map an index range in an ungapped sequence to indices in the gapped sequence
    """
    original_index = 0
    (a, b) = (-1, 0)
    for (i, c) in enumerate(fasta):
        if c not in "._-":
            original_index += 1
            if original_index == start:
                a = i + 1
            if original_index == end:
                b = i + 1
                break
    # truncate to sequence length if the end index is greater than sequence length
    else:
        b = i + 1
    if a == -1:
        # The start is beyond the end of the sequence
        # For now I will provide a dumby value that is guaranteed to index to an empty list
        dumby = len(fasta) + 1
        return (dumby, dumby)
    else:
        return (a, b)


def extract_aa2aa(
    bounds: List[Tuple[int, int]],
    seq: List[smof.FastaEntry],
    ref: List[smof.FastaEntry],
    mafft_opts: MafftOpts = MafftOpts(),
) -> List[List[smof.FastaEntry]]:
    # empty inputs silently return empty results
    if len(seq) > 0:

        # there must exactly 1 referece
        if len(ref) != 1:
            raise InputError("Expected exactly one reference file")

        # add reference sequences to the input seq
        with_ref = ref + seq

        # align the reference and input protein sequences
        aln = align(with_ref, mafft_opts)

        # find reference gapped start and stop positions
        intervals = [
            ungap_indices(start, stop, list(aln)[0].seq) for (start, stop) in bounds
        ]

        # extract the subset regions
        extracts = [smof.subseq(aln, start, stop) for (start, stop) in intervals]

        return [list(extracted)[1:] for extracted in extracts]

    else:

        return []


def map_dna2dna(
    bounds: List[Tuple[int, int]],
    fna: List[smof.FastaEntry],
    aln: List[smof.FastaEntry],
) -> List[List[smof.FastaEntry]]:

    motif_sets = []

    for (gapped_start, gapped_stop) in bounds:
        # find reference gapped start and stop positions
        (a, b) = ungap_indices(gapped_start, gapped_stop, list(aln)[0].seq)

        motif = []

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
            dna_end = dna_start + aa_length * 3

            motif.append(
                smof.FastaEntry(header=fna[i].header, seq=fna[i].seq[dna_start:dna_end])
            )

        motif_sets.append(motif)

    return motif_sets


def extract_dna2dna(
    bounds: List[Tuple[int, int]],
    seq: List[smof.FastaEntry],
    ref: List[smof.FastaEntry],
    mafft_opts: MafftOpts = MafftOpts(),
) -> List[List[smof.FastaEntry]]:
    # translate the DNA inputs (longest uninterupted CDS)
    faa = smof.translate(seq, all_frames=True)

    # because generators are the root of all evil
    ref = list(ref)
    faa = list(faa)
    seq = list(seq)
    aln = list(align(ref + faa, mafft_opts))

    extracted = map_dna2dna(bounds, seq, aln)

    return extracted


def _dispatch_extract(
    fasta_file: TextIO,
    bounds: List[Tuple[int, int]],
    seq_opts: SeqOpts = SeqOpts(),
    mafft_opts: MafftOpts = MafftOpts(),
) -> List[smof.FastaEntry]:
    if seq_opts.subtype is not None:
        # get reference for this subtype
        ref = get_ref(seq_opts.subtype)
    else:
        raise InputError(
            "HA subtype MUST be defiend to extract a range relative to a reference"
        )

    # open input sequences
    entries = list(smof.open_fasta(fasta_file))
    # remove gaps
    entries = list(smof.clean(entries, toseq=True))
    # dispatch by sequence type
    if seq_opts.conversion is None:
        raise InputError(
            "Conversion is undefined, did you forget a --conversion argument?"
        )
    else:
        if seq_opts.conversion == Conversion.AA_TO_AA:
            f = extract_aa2aa
        elif seq_opts.conversion == Conversion.DNA_TO_AA:
            entries = list(smof.translate(entries, all_frames=True))
            f = extract_aa2aa
        elif seq_opts.conversion == Conversion.DNA_TO_DNA:
            f = extract_dna2dna
        else:
            raise ValueError("Invalid Conversion value - this is a bug in flutile")

    motifs = f(bounds, entries, ref, mafft_opts)

    return motifs


def extract_ha1(
    fasta_file: TextIO,
    seq_opts: SeqOpts = SeqOpts(),
    mafft_opts: MafftOpts = MafftOpts(),
) -> List[smof.FastaEntry]:
    """
    Get the HA1 range relative to the subtype of interest by mapping the H1 HA1
    range (18,344) range to the subtype of interest
    """

    if seq_opts.subtype is not None:
        (start, end) = motifs.GENBANK_HA1_REGIONS[seq_opts.subtype]
    else:
        raise InputError("Please specify a subtype")

    outs = _dispatch_extract(
        fasta_file=fasta_file,
        bounds=[(start, end)],
        seq_opts=seq_opts,
        mafft_opts=mafft_opts,
    )

    # for ha1 extract, there will be exactly one region extracted for sequence
    out = list(outs)[0]
    return out


def print_fasta(seq: List[smof.FastaEntry], out: TextIO = sys.stdout) -> None:
    """
    This is just a little helper to avoid requiring the UI import smof
    """
    smof.print_fasta(seq, out=out)


def get_signal_offset(subtype: str) -> int:
    """
    subtype should be a string from H1 to H18
    """
    try:
        return len(motifs.NTERM_MOTIFS[subtype])
    except KeyError:
        raise InputError("Expected a subtype string (e.g., H1)")


def extract_bounds(
    fasta_file: TextIO,
    bounds: List[Tuple[int, int]],
    seq_opts: SeqOpts = SeqOpts(),
    mafft_opts: MafftOpts = MafftOpts(),
) -> List[Tuple[str, List[str]]]:
    """
    Extract a motif
    """
    bounds = [(min(xs), max(xs)) for xs in bounds]

    if not seq_opts.keep_signal:
        if seq_opts.subtype is None:
            raise InputError("Please specify a subtype")
        elif is_ha(seq_opts.subtype):
            offset = get_signal_offset(seq_opts.subtype)
            bounds = [(a + offset, b + offset) for (a, b) in bounds]
        else:
            # only HA has a signal peptide to trim (right?)
            pass

    extracts = [
        list(seq)
        for seq in _dispatch_extract(
            fasta_file=fasta_file,
            bounds=bounds,
            seq_opts=seq_opts,
            mafft_opts=mafft_opts,
        )
    ]

    pairs = []
    if len(extracts) > 0:
        for n_seq in range(len(extracts[0])):
            defline = extracts[0][n_seq].header
            motif_seqs = [motif[n_seq].seq for motif in extracts]
            pairs.append((defline, motif_seqs))

    return pairs


def parse_motif(
    motif_str: str, subtype: Optional[str]
) -> Tuple[str, List[Tuple[int, int]]]:
    """
    Parse a motif string provided by the user through the `flutile trim motif` command.

    Here is the string grammar:

    motif = range
          | name "=" range
    range = number
          | number "-" number
          | range "," range

    Ranges are relative to indices on the reference sequecne, so subtype is
    ultimately NOT optional. However, this function is not responsible for
    looking up references; missing subtypes will raise an error further
    downstream.

    Examples:
      * 42
      * 42-44
      * 33,42-44
      * foofoo=33,42-44
    """
    motif_str = motif_str.strip()
    pair = motif_str.split("=")
    if len(pair) == 1:
        if subtype is None:
            name = motif_str
        else:
            name = subtype + ":" + motif_str
        ranges_str = pair[0]
    elif len(pair) == 2:
        (name, ranges_str) = pair
    else:
        raise InputError("There should be only one equal sign per motif")
    name = name.strip()
    submotif_strs = ranges_str.strip().split(",")

    bounds: List[Tuple[int, int]]
    bounds = []
    for submotif_str in [s.strip() for s in submotif_strs]:
        pair = submotif_str.split("-")
        if len(pair) == 1:
            (a, b) = (submotif_str, submotif_str)
        elif len(pair) == 2:
            (a, b) = pair
        else:
            raise InputError("Badly formed motif string")
        bounds.append((int(a), int(b)))

    return (name, bounds)


def concat(xss: Iterable[Iterable[A]]) -> List[A]:
    return [x for xs in xss for x in xs]


def unconcat(
    xs: List[str],
    widths: List[int],
    joiner: Callable[[List[str]], str] = lambda ys: "".join(ys),
) -> List[str]:
    i = 0
    collapsed = []
    for w in [w for w in widths if w > 0]:
        if i >= len(xs):
            break
        collapsed.append(joiner(xs[i : (i + w)]))
        i += w
    return collapsed


def make_motifs(
    fasta_file: TextIO,
    motif_strs: List[str],
    seq_opts: SeqOpts = SeqOpts(),
    mafft_opts: MafftOpts = MafftOpts(),
):

    motifs: Dict[str, List[Tuple[int, int]]] = dict()
    motifs = {
        name: bounds
        for (name, bounds) in [
            parse_motif(motif, subtype=seq_opts.subtype) for motif in motif_strs
        ]
    }

    bounds: List[Tuple[int, int]] = concat(motifs.values())
    pairs_flat = extract_bounds(
        fasta_file, bounds, seq_opts=seq_opts, mafft_opts=mafft_opts
    )

    widths = [len(v) for v in motifs.values()]
    pairs = []
    for (defline, segs) in pairs_flat:
        pairs.append((defline, unconcat(segs, widths)))

    return (motifs.keys(), pairs)


def gapped_indices(seq: str) -> List[str]:
    """
    Given the gapped reference sequence, create index strings for use in aadiff site columns.

    Any gaps before the first reference base translate to negative numbers
    specifying the negative offset from the first reference base. Any gaps
    after a reference base translate to strings specifying the positive offset.

    For example:
      "--AT--G-"  -->  ["-2", "-1", "1", "2", "2+1", "2+2", "3", "3+1"]
    """

    def _handle_run(run: int, ind: int) -> List[str]:
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
    msa: List[Tuple[str, str]], aadiff_opts: AadiffOpts = AadiffOpts()
) -> List[List[str]]:
    """
    From the input alignment s, make an aadiff table
    """

    def find_consensus(i):
        counts = collections.Counter(msa[j][1][i] for j in range(1, len(msa)))
        return counts.most_common()[0][0]

    # add the consensus header column, if needed
    if aadiff_opts.make_consensus or aadiff_opts.consensus_as_reference:
        consensus_seq = "".join(find_consensus(i) for i in range(len(msa[0][1])))
        # the consensus column goes first if it is being used as a reference
        if aadiff_opts.consensus_as_reference:
            msa = [("Consensus", consensus_seq)] + msa
        # otherwise it goes last
        else:
            msa = msa + [("Consensus", consensus_seq)]

    header = ["site"] + [x[0] for x in msa]

    table = [header]

    seq_ids = list(range(len(msa)))
    aa_ids = list(range(len(msa[0][1])))

    # for each amino acid position in the alignment
    for i in aa_ids:
        ref = msa[0][1][i]
        position = str(i + 1)
        # for each sequence in the alignment
        for j in seq_ids:
            if not aadiff_opts.remove_unchanged or msa[j][1][i] != ref:
                row = [position, ref]
                for k in seq_ids[1:]:
                    aa = msa[k][1][i]
                    if aa == ref:
                        row.append("")
                    else:
                        row.append(aa)
                table.append(row)
                break

    return table


def transpose(xss: List[List[A]]) -> List[List[A]]:
    """
    Transpose a list of lists. Each element in the input list of lists is
    considered to be a column. The list of rows is returned.

    The nested lists are assumed to all be of the same length. If they are not,
    the shortest list will be used as the number of output rows and other rows
    wil be lost without raising an error.
    """
    # get the number of rows in the input
    N = min([len(xs) for xs in xss])
    return [[xs[i] for xs in xss] for i in range(N)]


def annotate_table(
    table: List[List[str]],
    seq_opts: SeqOpts = SeqOpts(),
    aadiff_opts: AadiffOpts = AadiffOpts(),
) -> List[List[str]]:
    def load_builtin_annotations(name: str, subtype_exp: str) -> str:
        if subtype_exp != seq_opts.subtype:
            raise InputError(f"{name} annotation is only defined for {subtype_exp}")

        if seq_opts.keep_signal:
            raise InputError("{name} annotation is incompatible with --keep_signal")

        return os.path.join(os.path.dirname(__file__), "data", f"{name}.txt")

    all_anntables = []

    if aadiff_opts.caton82:
        all_anntables.append(load_builtin_annotations("caton82", "H1"))

    if aadiff_opts.wiley81:
        all_anntables.append(load_builtin_annotations("wiley81", "H3"))

    all_anntables += aadiff_opts.annotation_tables

    annotations: Dict[str, List[str]]
    annotations = {x[0]: [] for x in table[1:]}

    new_column_names = []

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

    if aadiff_opts.count:
        # aggregate all non-index rows into a single counts column
        for i in range(len(table)):
            ag = ", ".join(
                [
                    f"{str(v)} {k}"
                    for (k, v) in collections.Counter(table[i][2:]).items()
                    if k
                ]
            )
            table[i] = [table[i][0], table[i][1], f"({ag})"]
        new_column_names = ["site", "reference", "changes"]

    if aadiff_opts.join_annotations:
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
    faa: TextIO,
    seq_opts: SeqOpts = SeqOpts(),
    mafft_opts: MafftOpts = MafftOpts(),
    aadiff_opts: AadiffOpts = AadiffOpts(),
) -> List[List[str]]:
    # if this is an HA and we want to trim off the signal peptide
    if (
        seq_opts.subtype is not None
        and is_ha(seq_opts.subtype)
        and not seq_opts.keep_signal
    ):
        ref = motifs.get_ha_subtype_nterm_motif(seq_opts.subtype)
        # the motif here is an exact match to the reference signal peptide, we
        # want to remove the signal peptide, so the trim length is simply the
        # peptide length
        trim = len(ref.motif)

        # trim off signal peptides or other gunk at the beginning of the reference
        seq = ref.sequence[trim:]

        refseq = [smof.FastaEntry(header=ref.defline, seq=seq)]
    # if we aren't trimming, we can use the default (M-initialized) full protein references
    elif seq_opts.subtype:
        refseq = get_ref(seq_opts.subtype)
    # unless we aren't given a reference, then we do nothing special
    else:
        refseq = []

    # open input sequences
    entries = smof.open_fasta(faa)

    # remove gaps
    entries = list(smof.clean(entries, toseq=True))

    # prepend any references
    entries = refseq + entries

    # align the reference and input protein sequences
    aln_objs = align(entries, mafft_opts)
    # To simple list of string pairs, rather than list of FastaEntry objects
    aln = [(s.header, s.seq) for s in aln_objs]

    try:
        indices = gapped_indices(aln[0][1])
    except IndexError:
        raise InputError("Empty fasta sequence")

    # remove any reference sequences (1 or 0)
    aln = aln[len(refseq) :]

    table = aadiff_table(aln, aadiff_opts)

    # set relative indices
    for i, row in enumerate(table):
        if i > 0:
            table[i][0] = indices[int(row[0]) - 1]

    return table


def referenced_aadiff_table(
    faa: TextIO,
    mafft_opts: MafftOpts = MafftOpts(),
    seq_opts: SeqOpts = SeqOpts(),
    aadiff_opts: AadiffOpts = AadiffOpts(),
) -> List[List[str]]:
    table = referenced_table(
        faa, mafft_opts=mafft_opts, seq_opts=seq_opts, aadiff_opts=aadiff_opts
    )
    return annotate_table(table, aadiff_opts=aadiff_opts, seq_opts=seq_opts)


def map_ha_range(
    start: int, end: int, subtype1: int, subtype2: int
) -> Tuple[Optional[int], Optional[int]]:
    """
    Map AA indices between HA subtypes using the Burke2014 index map

    If the region maps off the edge of the subtype2, then return (None, None).

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


def referenced_annotation_table(
    faa: TextIO,
    mafft_opts: MafftOpts = MafftOpts(),
    seq_opts: SeqOpts = SeqOpts(),
    aadiff_opts: AadiffOpts = AadiffOpts(),
) -> List[List[str]]:
    """
    This is just for HA
    """

    table = referenced_table(faa, seq_opts, mafft_opts, aadiff_opts)

    subtype_file = os.path.join(
        os.path.dirname(__file__), "data", "burke2014-index-map.txt"
    )

    if seq_opts.subtype is None:
        raise InputError("Please provide a subtype")
    elif not is_ha(seq_opts.subtype):
        raise InputError("Subtype must be an HA")
    else:
        try:
            subtype_number = int(seq_opts.subtype[1:])
        except TypeError:
            raise InputError("HA subtype must have the form 'H<number>', e.g., H1, H17")

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

    return annotate_table(table, seq_opts, aadiff_opts)
