import memfiles
import streams
import zip/gzipfiles
import utils/helpers
import strformat
import strutils
import os
import hts


proc vcf_window*(vcf: string) =
    var v:VCF
    doAssert open(v, vcf)
    #echo $v.info