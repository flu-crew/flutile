import signal
import os
import click
from flutile.version import __version__
import sys
import click
from flutile.functions import *

INT_SENTINEL = 1e9

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.command(
    name="h1",
    help="Compare AA differences between H1 strains. The input fasta file does NOT need to be aligned. Indexing is relative to the pandemic A/California/07/2009 strain. Optional trimming removes the first 17 residues from the indexing reference. This affects only the index numbering. If you actually want to trim the amino acid sequences, see `flutile trim`.",
)
@click.argument("faa", default=sys.stdin, type=click.File())
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
def h1_aadiff_cmd(faa, keep_signal, make_consensus, consensus_as_reference):
    table = h1_aadiff_table(
        faa_file=faa,
        keep_signal=keep_signal,
        consensus=make_consensus,
        consensus_as_ref=consensus_as_reference,
    )
    for row in table:
        print("\t".join(row))

@click.command(
    name="h3",
    help="Compare AA differences between H3 strains. The input fasta file does NOT need to be aligned. Indexing is relative to A/Brisbane/10/2007.",
)
@click.argument("faa", default=sys.stdin, type=click.File())
@click.option("--make-consensus", is_flag=True, help="Add a sequence consensus column")
@click.option(
    "--consensus-as-reference",
    is_flag=True,
    help="Use the consensus as the reference column",
)
def h3_aadiff_cmd(faa, keep_signal, make_consensus, consensus_as_reference):
    table = h3_aadiff_table(
        faa_file=faa,
        keep_signal=keep_signal,
        consensus=make_consensus,
        consensus_as_ref=consensus_as_reference,
    )
    for row in table:
        print("\t".join(row))


@click.group(
    name="aadiff",
    help="Compare AA differences between H1 strains.",
    context_settings=CONTEXT_SETTINGS,
)
def aadiff_grp():
    pass

aadiff_grp.add_command(h1_aadiff_cmd)
aadiff_grp.add_command(h3_aadiff_cmd)



#      flutile represent [--max-day-sep=<days>] [--min-pident-sep=<pident>] [--same-state] [--print-groups] [<alignment>]
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


mafft_exe_opt = click.option(
    "--mafft-exe",
    default="mafft",
    type=str,
    help="Path to MAFFT alignment tool executable",
)

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
def h3_ha1_cmd(fasta_file, mafft_exe, conversion_opt):
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


cli_grp.add_command(aadiff_grp)
cli_grp.add_command(represent_cmd)
cli_grp.add_command(trim_grp)


def main():
    cli_grp()


if __name__ == "__main__":
    if os.main is "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()
