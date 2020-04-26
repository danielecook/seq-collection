import hts
import json
import sequtils
import strutils
import strformat
import utils/helpers
import utils/vcf_header
import streams
from constants import ANN_header

proc vcf2sql*(vcf: string) =
    var v:VCF
    doAssert open(v, vcf)
