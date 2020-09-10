import signal
import os
import click
from flutile.version import __version__
import sys
import click
from flutile.functions import *

INT_SENTINEL = 1e9


@click.command(
    name="compare",
    help="Build table comparing AA differences between strains. The input must be a sequence alignment, usually a protein alignment, but DNA will work too.",
)
@click.argument("alignment", default=sys.stdin, type=click.File())
@click.option("--make-consensus", is_flag=True, help="Add a sequence consensus column")
@click.option(
    "--consensus-as-reference",
    is_flag=True,
    help="Use the consensus as the reference column",
)
def compare_cmd(alignment, make_consensus, consensus_as_reference):
    s = parseFasta(alignment)
    rows = aadiff_table(
        s, consensus=make_consensus, consensus_as_ref=consensus_as_reference
    )
    for row in rows:
        print("\t".join(row))


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
    s = parseFasta(alignment)

    if max_day_sep == INT_SENTINEL:
        max_day_sep = None

    (groups, seqs) = represent(
        s, max_day_sep=max_day_sep, min_pident_sep=min_pident_sep, same_state=same_state
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
cds_opt = click.option(
    "--cds",
    default="mafft",
    is_flag=True,
    type=str,
    help="Trim to DNA coding sequence (not AA)",
)

@click.command(name="h1-ha1", help="Trim H1 DNA down to the HA1 AA using Brisbane/10/2007 template.")
@click.argument("fasta_file", default=sys.stdin, type=click.File())
@mafft_exe_opt
@cds_opt
def h1_ha1_cmd(fasta_file, mafft_exe, cds):
    extract_h1_ha1(fasta_file, mafft_exe=mafft_exe, cds=cds)

@click.command(name="h3-ha1", help="Trim H3 DNA down to the HA1 AA using California/07/2009 template.")
@click.argument("fasta_file", default=sys.stdin, type=click.File())
@mafft_exe_opt
@cds_opt
def h3_ha1_cmd(fasta_file, mafft_exe, cds):
    extract_h3_ha1(fasta_file, mafft_exe=mafft_exe, cds=cds)


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

@click.group(name="trim", help="Trim flu sequences in various ways. The trim operations use subtype-specific templates that are stored in the flutile package. You should not need to change these.", context_settings=CONTEXT_SETTINGS)
def trim_grp():
    pass

trim_grp.add_command(h1_ha1_cmd)
trim_grp.add_command(h3_ha1_cmd)

@click.group(help="Flu-crew utilities", context_settings=CONTEXT_SETTINGS)
def cli_grp():
    pass

cli_grp.add_command(compare_cmd)
cli_grp.add_command(represent_cmd)
cli_grp.add_command(trim_grp)


def main():
    cli_grp()


if __name__ == "__main__":
    if os.main is "posix":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    main()
