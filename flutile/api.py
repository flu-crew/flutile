# This module contains functions that are exported as part of the API that
# other packages may use. Specifically octofludb. Any change to these functions
# may break packages that depend on them, so take care.

from __future__ import annotations
from typing import TextIO, List, Union
import flutile.parameters as par
import flutile.functions as fun
import sys


def write_bounds(
    motif_strs: List[str],
    subtype: str,
    fasta_file: Union[str, TextIO],
    tabular: bool = True,
    keep_signal: bool = False,
    conversion: str = "dna2aa",
    outfile: Union[str, TextIO] = sys.stdout,
    mafft_exe: str = "mafft",
    verbose: bool = False,
) -> None:
    """
    Print motif bounds. It is used by octofludb.
    """

    if isinstance(fasta_file, str):
        fasta_fh = open(fasta_file, "r")
    else:
        fasta_fh = fasta_file

    if isinstance(outfile, str):
        outfile_fh = open(outfile, "w")
    else:
        outfile_fh = outfile

    conversion_enum = par.parse_conversion(conversion)

    mafft_opts = par.MafftOpts(mafft_exe=mafft_exe, verbose=verbose)
    seq_opts = par.SeqOpts(
        keep_signal=keep_signal, conversion=conversion_enum, subtype=subtype
    )

    (names, pairs) = fun.make_motifs(
        fasta_file=fasta_fh,
        motif_strs=motif_strs,
        seq_opts=seq_opts,
        mafft_opts=mafft_opts,
    )

    if tabular:
        print("\t" + "\t".join(names), file=outfile_fh)
        for (defline, seqs) in pairs:
            print(defline + "\t" + "\t".join(seqs), file=outfile_fh)
    else:
        pairs = [(header, "".join(seqs)) for (header, seqs) in pairs]
        fun.print_fasta(pairs, out=outfile_fh)
