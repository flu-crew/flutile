from collections import namedtuple
import os
import smof
import re
import sys

RefMotif = namedtuple("RefMotif", "defline sequence motif")


NTERM_MOTIFS = {
    "H1": "MKARLLVLLCALAATDA",
    "H2": "MAIIYLILLFTAVRG",
    "H3": "MKTIIALSYIFCLALG",  # the paper had the motif MKTIIALSYIFCLPLG, but this seems to be a mistake
    "H4": "MLSIAILFLLIAEGSS",
    "H5": "MEKIVLLFAIVSLVKS",
    "H6": "MIAIIVIATLAAAGKS",
    "H7": "MNTQILVFALVASIPTNA",
    "H8": "MEKFIAIAMLLASTNA",
    "H9": "MEAASLITILLVVTASNA",
    "H10": "MYKIVVIIALLGAVKG",
    "H11": "MEKTLLFAAIFLCVKA",
    "H12": "MEKFIILSTVLAASFA",
    "H13": "MALNVIATLTLISVCVHA",
    "H14": "MIALILVALALSHTAYS",
    "H15": "MNTQIIVILVLGLSMVRS",
    "H16": "MMIKVLYFLIIVLGRYSKA",
    "H17": "MELIILLILLNPYTFVLG",
    "H18": "MITILILVLPIVVG",
}


def get_ha_subtype_nterm_motif(ha_subtype):
    """
    ha_subtype: H1, H2, ..., H18
    """
    ref_file = "subtype-refs.faa"
    subtype_faa = os.path.join(os.path.dirname(__file__), "data", ref_file)
    ref_motif = None
    for entry in smof.open_fasta(subtype_faa):
        try:
            subtype = re.match(".*(H[0-9]+)N[0-9]+.*", entry.header).groups(1)[0]
        except AttributeError:
            print(
                "Expected the subtype to be in the '{ref_file}' defline",
                file=sys.stderr,
            )
            raise
        try:
            motif = NTERM_MOTIFS[subtype]
        except KeyError:
            print(f"Unexpected HA, found {subtype}, expected H1-H18", file=sys.stderr)
            raise

        if not motif in entry.seq:
            print(
                f"Could not find the motif '{motif}' in {subtype} refernce in '{ref_file}'"
            )
            raise

        if subtype == ha_subtype:
            ref_motif = RefMotif(defline=entry.header, sequence=entry.seq, motif=motif)
            # Technically I could stop the loop here, but the cost of checking
            # the whole file is not great enough to offset the value of
            # checking the correctness of the data.

    return ref_motif
