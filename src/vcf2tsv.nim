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

proc vcf2tsv*(vcf_in: string) =
    echo vcf_in