import signal
import os
import click
from flutile.version import __version__
import sys
import click
from flutile.functions import *

INT_SENTINEL = 1e9

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

mafft_exe_opt = click.option(
    "--mafft-exe",
    default="mafft",
    type=str,
    help="Path to MAFFT alignment tool executable",
)

@click.command(
    name="aadiff",
    help="Compare differences between sequences. The input fasta file does NOT need to be aligned. If a subtype is specified, indexing will be relative to the index reference (Burke 2014). If no subtype is specified, the numbering will be relative to the first entry in fasta file (or the consensus sequence if --consensus-as-reference is set)",
)
@click.argument("faa", default=sys.stdin, type=click.File())
@click.option(
    "--subtype",
    help="Currently HA subtypes from H1 to H18 are supported and will number relative to the start of the mature peptide, using the offsets described in (Burke 2014). If the flag --keep-signal is set, then numbering is relative to the initial methionine.",
)
@click.option("--make-consensus", is_flag=True, help="Add a sequence consensus column")
@click.option(
    "--consensus-as-reference",
    is_flag=True,
    help="Use the consensus as the reference column",
)
@click.option(
    "--keep_signal",
    is_flag=True,
    help="Number relative to the initial Methionine (do not trim the signal peptide)",
)
@mafft_exe_opt
def aadiff_cmd(*args, **kwargs):
    for row in referenced_aadiff_table(*args, **kwargs):
        print("\t".join(row))


@click.command(
    name="represent",
    help="Representative subsampling by sequence similarity. The input must be a sequence alignment.",
)
@click.argument("alignment", default=sys.stdin, type=click.File())
@click.option(
    "--max-day-sep",
    help="Maximum number of days separating members of a group",
    type=click.IntRange(min=0, max=INT_SENTINEL),
    default=INT_SENTINEL,
)
@click.option(
    "--min-pident-sep",
    help="Minimum proportion identity difference between members of a group",
    type=click.FloatRange(min=0, max=1),
    default=1,
)
@click.option(
    "--same-state", help="Group strains only if they are in the same sate", is_flag=True
)
@click.option(
    "--print-groups",
    help="Rather than subsetting the fasta, print the groups of similar strains",
    is_flag=True,
)
def represent_cmd(alignment, max_day_sep, min_pident_sep, same_state, print_groups):
    if max_day_sep == INT_SENTINEL:
        max_day_sep = None

    (groups, seqs) = with_aligned_pairs(
        alignment,
        represent,
        max_day_sep=max_day_sep,
        min_pident_sep=min_pident_sep,
        same_state=same_state,
    )

    if print_groups:
        for group in groups:
            for i in group:
                print(s[i][0])
            print("")
    else:
        for (header, seq) in seqs:
            print(">" + header)
            print(seq)


conversion_opt = click.option(
    "--conversion",
    type=click.Choice(["dna2aa", "dna2dna", "aa2aa"], case_sensitive=False),
    default="dna2aa",
    help="aa2aa: align the input and reference AA sequences and extract the HA1 regions based on the reference. dna2aa: translate the input DNA sequences, add the AA reference, align, and extract the HA1. dna2dna: translate the input DNA (keeping track of the CDS start position), add reference, align, find AA HA1 motif, map HA1 regions back to DNA.",
)


@click.command(
    name="h1-ha1",
    help="Trim H1 DNA down to the HA1 AA using Brisbane/10/2007 template.",
)
@click.argument("fasta_file", default=sys.stdin, type=click.File())
@mafft_exe_opt
@conversion_opt
def h1_ha1_cmd(fasta_file, mafft_exe, conversion):
    extract_h1_ha1(fasta_file, mafft_exe=mafft_exe, conversion=conversion)


@click.command(
    name="h3-ha1",
    help="Trim H3 DNA down to the HA1 AA using California/07/2009 template.",
)
@click.argument("fasta_file", default=sys.stdin, type=click.File())
@mafft_exe_opt
@conversion_opt
def h3_ha1_cmd(fasta_file, mafft_exe, conversion):
    extract_h3_ha1(fasta_file, mafft_exe=mafft_exe, conversion=conversion)


@click.group(
    name="trim",
    help="Trim flu sequences in various ways. The trim operations use subtype-specific templates that are stored in the flutile package. You should not need to change these.",
    context_settings=CONTEXT_SETTINGS,
)
def trim_grp():
    pass


trim_grp.add_command(h1_ha1_cmd)
trim_grp.add_command(h3_ha1_cmd)


@click.group(help="Flu-crew utilities", context_settings=CONTEXT_SETTINGS)
def cli_grp():
    pass


cli_grp.add_command(aadiff_cmd)
cli_grp.add_command(represent_cmd)
cli_grp.add_command(trim_grp)


def main():
    cli_grp()


if __name__ == "__main__":
    if os.main is "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()
