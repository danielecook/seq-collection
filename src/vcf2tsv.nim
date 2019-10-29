import tables
import sequtils
import strutils
import utils/gz
import os
import re
import sets
import streams
import zip/gzipfiles
import utils/helpers
import bitvector
import bloom
import hts

proc vcf2tsv*(vcf: string) =
    var v:VCF
    doAssert open(v, vcf)