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

proc vcf2tsv*(vcf: string, long: bool, region_list: seq[string]) =
    var v:VCF
    doAssert open(v, vcf)
    var samples = v.samples
    for rec in variants(v, region_list):
        echo rec
        