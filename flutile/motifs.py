from __future__ import annotations
from collections import namedtuple
import os
import smof  # type: ignore
import re
import sys
from typing import Dict, Tuple

RefMotif = namedtuple("RefMotif", "defline sequence motif")

# These bounds are calibrated to reproduce the HA1 regions that are reported in
# genbank. They are tested against random selections from genbank.
GENBANK_HA1_REGIONS: Dict[str, Tuple[int, int]] = {
    "H1": (18, 344),
    "H2": (16, 340),
    "H3": (17, 345),
    "H4": (17, 343),
    "H5": (17, 346),
    "H6": (17, 345),
    "H7": (19, 341),
    "H8": (18, 344),
    "H9": (19, 338),
    "H10": (17, 340),
    "H11": (17, 342),
    "H12": (18, 342),
    "H13": (19, 343),
    "H14": (18, 347),
    "H15": (19, 349),
    "H16": (20, 344),
    "H17": (19, 342),
    "H18": (15, 339),
}


NTERM_MOTIFS: Dict[str, str] = {
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


def get_ha_subtype_nterm_motif(ha_subtype: str) -> RefMotif:
    """
    ha_subtype: H1, H2, ..., H18
    """
    ref_file = "subtype-refs.faa"
    subtype_faa = os.path.join(os.path.dirname(__file__), "data", ref_file)
    ref_motif = None
    for entry in smof.open_fasta(subtype_faa):
        subtype_match = re.match(".*(H[0-9]+)N[0-9]+.*", entry.header)
        subtype: str
        if subtype_match is None:
            raise ValueError("Expected the subtype to be in the '{ref_file}' defline")
        else:
            subtype = subtype_match.groups(1)[0]  # type: ignore

        try:
            motif = NTERM_MOTIFS[subtype]
        except KeyError:
            print(f"Unexpected HA, found {subtype}, expected H1-H18", file=sys.stderr)
            raise

        if motif not in entry.seq:
            print(
                f"Could not find the motif '{motif}' in {subtype} refernce in '{ref_file}'"
            )
            raise

        if subtype == ha_subtype:
            ref_motif = RefMotif(defline=entry.header, sequence=entry.seq, motif=motif)
            # Technically I could stop the loop here, but the cost of checking
            # the whole file is not great enough to offset the value of
            # checking the correctness of the data.

    if ref_motif is None:
        raise ValueError("Expected an HA subtype such as 'H1' or 'H12'")
    else:
        return ref_motif
