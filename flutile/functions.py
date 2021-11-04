import collections
import sys
import datetime
import re
import tqdm
import smof
import os
import flutile.motifs as motifs


class InputError(Exception):
    pass


def err(msg):
    raise InputError(msg)


def is_aligned(fasta):
    lengths = {len(s.seq) for s in fasta}
    if len(lengths) != 1:
        err("FASTA file is not aligned, all entries should be of equal length")


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

    groups = groups + [{i} for i in range(N) if i not in grouped]

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
            if y not in xs:
                xs.add(y)
                xs.update(group(y, xs, xss))
        return xs

    groups = []
    while groupmap:
        (a, b) = list(groupmap.items())[0]
        xs = group(a, set(), groupmap)
        groupmap = {k: v for (k, v) in groupmap.items() if k not in xs}
        groups.append(xs)

    return groups


def is_ha(subtype):
    return bool(re.fullmatch("H\d+", subtype))


def is_na(subtype):
    return bool(re.fullmatch("N\d+", subtype))


def get_ref(subtype):
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
            err("Expected HA, but H was not followed by an integer")
        ref_file = os.path.join(os.path.dirname(__file__), "data", "subtype-refs.faa")
        ref_fasta = smof.open_fasta(ref_file)
        ref = list(smof.grep(ref_fasta, no_color=True, pattern=f"|H{i}N"))
    elif is_na(subtype):
        try:
            i = int(subtype[1:])
        except:
            err("Expected NA, but N was not followed by an integer")
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
        err("Unexpected subtype")
    return ref


def align(seq, mafft_exe="mafft"):
    from Bio.Align.Applications import MafftCommandline
    from tempfile import mkstemp

    hashy = smof.md5sum(seq)
    (n, fasta_filename) = mkstemp(suffix=f"-{hashy}.fa")

    with open(fasta_filename, "w+") as fasta_fh:
        smof.print_fasta(seq, out=fasta_fh)

    o, e = MafftCommandline(mafft_exe, input=fasta_filename)()

    print(e, file=sys.stderr)

    (n, aln_filename) = mkstemp(suffix=f"-{hashy}.aln")

    with open(aln_filename, "w+") as fh:
        print(o, file=fh)

    aln = list(smof.open_fasta(aln_filename))

    return aln


def ungap_indices(start, end, fasta):
    """
    Map indices in an ungapped sequence to indices in the gapped string
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


def extract_aa2aa(bounds, faa, ref, mafft_exe="mafft"):
    faa = list(faa)
    ref = list(ref)

    # empty inputs silently return empty results
    if (len(faa) > 0):

      # there must exactly 1 referece
      if len(ref) != 1:
          err("Expected exactly one reference file")

      # add reference sequences to the input seq
      with_ref = ref + faa

      # align the reference and input protein sequences
      aln = align(with_ref, mafft_exe=mafft_exe)

      # find reference gapped start and stop positions
      intervals = [
          ungap_indices(start, stop, list(aln)[0].seq) for (start, stop) in bounds
      ]

      # extract the subset regions
      extracts = [smof.subseq(aln, start, stop) for (start, stop) in intervals]

      return [list(extracted)[1:] for extracted in extracts]

    else:

      return []

def map_dna2dna(bounds, fna, aln):

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


def extract_dna2dna(bounds, fna, ref, mafft_exe="mafft"):
    # translate the DNA inputs (longest uninterupted CDS)
    faa = smof.translate(fna, all_frames=True)

    # because generators are the root of all evil
    ref = list(ref)
    faa = list(faa)
    fna = list(fna)
    aln = list(align(ref + faa, mafft_exe=mafft_exe))

    extracted = map_dna2dna(bounds, fna, aln)

    return extracted


def _dispatch_extract(bounds, subtype, fasta_file, conversion=None, *args, **kwargs):
    # get reference for this subtype
    ref = get_ref(subtype)

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

    motifs = f(bounds, entries, ref, *args, **kwargs)

    return motifs


def extract_ha1(subtype, *args, **kwargs):
    """
    Get the HA1 range relative to the subtype of interest by mapping the H1 HA1
    range (18,344) range to the subtype of interest
    """

    (start, end) = motifs.GENBANK_HA1_REGIONS[subtype]

    outs = _dispatch_extract(bounds=[(start, end)], subtype=subtype, *args, **kwargs)

    # for ha1 extract, there will be exactly one region extracted for sequence
    out = list(outs)[0]
    smof.print_fasta(out)


def get_signal_offset(subtype):
    """
    subtype should be a string from H1 to H18
    """
    try:
        return len(motifs.NTERM_MOTIFS[subtype])
    except KeyError:
        err("Expected a subtype string (e.g., H1)")


def extract_bounds(bounds, keep_signal, subtype, *args, **kwargs):
    """
    Extract a motif
    """
    bounds = [(min(xs), max(xs)) for xs in bounds]

    if is_ha(subtype) and not keep_signal:
        offset = get_signal_offset(subtype)
        bounds = [(a + offset, b + offset) for (a, b) in bounds]

    extracts = [
        list(seq)
        for seq in _dispatch_extract(bounds=bounds, subtype=subtype, *args, **kwargs)
    ]

    pairs = []
    if len(extracts) > 0:
        for n_seq in range(len(extracts[0])):
            defline = extracts[0][n_seq].header
            motif_seqs = [motif[n_seq].seq for motif in extracts]
            pairs.append((defline, motif_seqs))

    return pairs


def parse_motifs(motif_strs, subtype):
    motif_set = dict()
    for (i, motif_str) in enumerate(motif_strs):
        motif_str = motif_str.strip()
        pair = motif_str.split("=")
        if len(pair) == 1:
            name = subtype + ":" + motif_str
            ranges_str = pair[0]
        elif len(pair) == 2:
            (name, ranges_str) = pair
        else:
            raise InputError("There should be only one equal sign per motif")
        name = name.strip()
        submotif_strs = ranges_str.strip().split(",")
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
        motif_set[name] = bounds
    return motif_set


def concat(xss):
    return [x for xs in xss for x in xs]


def unconcat(xs, widths, joiner=lambda ys: "".join(ys)):
    i = 0
    collapsed = []
    for w in [w for w in widths if w > 0]:
        if i >= len(xs):
            break
        collapsed.append(joiner(xs[i : (i + w)]))
        i += w
    return collapsed


def make_motifs(motif_strs, subtype, *args, **kwargs):
    motifs = parse_motifs(motif_strs=motif_strs, subtype=subtype)
    bounds = concat(motifs.values())
    pairs_flat = extract_bounds(bounds, subtype=subtype, *args, **kwargs)

    widths = [len(v) for v in motifs.values()]
    pairs = []
    for (defline, segs) in pairs_flat:
        pairs.append((defline, unconcat(segs, widths)))

    return (motifs.keys(), pairs)


def write_bounds(tabular=False, outfile=sys.stdout, *args, **kwargs):
    (names, pairs) = make_motifs(*args, **kwargs)

    if isinstance(outfile, str):
        outfile = open(outfile, "w")

    if tabular:
        print("\t" + "\t".join(names), file=outfile)
        for (defline, seqs) in pairs:
            print(defline + "\t" + "\t".join(seqs), file=outfile)
    else:
        pairs = [(header, "".join(seqs)) for (header, seqs) in pairs]
        smof.print_fasta(pairs, out=outfile)


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


def transpose(xss):
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
    table,
    subtype=None,
    annotation_tables="",
    join_annotations=False,
    keep_signal=False,
    count=False,
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

    if count:
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
    # if this is an HA and we want to trim off the signal peptide
    if is_ha(subtype) and not keep_signal:
        ref = motifs.get_ha_subtype_nterm_motif(subtype)
        # the motif here is an exact match to the reference signal peptide, we
        # want to remove the signal peptide, so the trim length is simply the
        # peptide length
        trim = len(ref.motif)

        # trim off signal peptides or other gunk at the beginning of the reference
        seq = ref.sequence[trim:]

        refseq = [smof.FastaEntry(header=ref.defline, seq=seq)]
    # if we aren't trimming, we can use the default (M-initialized) full protein references
    elif subtype:
        refseq = get_ref(subtype)
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
    aln = align(entries, mafft_exe=mafft_exe)
    aln = [(s.header, s.seq) for s in aln]

    indices = gapped_indices(aln[0][1])

    # remove any reference sequences (1 or 0)
    aln = aln[len(refseq) :]

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
    Map AA indices between HA subtypes using the Burke2014 index map

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
    """
    This is just for HA
    """

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
