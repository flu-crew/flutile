from typing import List, Optional
from enum import Enum


class Conversion(Enum):
    DNA_TO_DNA = 1
    AA_TO_AA = 2
    DNA_TO_AA = 3


class MafftOpts:
    """
    Options for the msa aligner
    """

    def __init__(self, mafft_exe: str = "mafft", verbose: bool = False):
        self.mafft_exe = mafft_exe
        self.verbose = verbose


class SeqOpts:
    """
    Biological options
    """

    def __init__(
        self,
        subtype: Optional[str] = None,
        keep_signal: bool = False,
        conversion: Optional[Conversion] = None,
    ):
        self.subtype = subtype
        self.keep_signal = keep_signal
        self.conversion = conversion


class AadiffOpts:
    """
    Options for generation of aadiff tables
    """

    def __init__(
        self,
        make_consensus: bool = False,
        consensus_as_reference: bool = False,
        remove_unchanged: bool = True,
        annotation_tables: List[str] = [],
        count: bool = False,
        join_annotations: bool = False,
        caton82: bool = False,
        wiley81: bool = False,
    ):
        self.make_consensus = make_consensus
        self.consensus_as_reference = consensus_as_reference
        self.remove_unchanged = remove_unchanged
        self.annotation_tables = annotation_tables
        self.count = count
        self.join_annotations = join_annotations
        self.caton82 = caton82
        self.wiley81 = wiley81
