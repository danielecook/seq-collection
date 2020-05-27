import utils/helpers
import strformat
import algorithm
import strutils
import sequtils
import os
import hts
import math
import threadpool
import random
import algorithm

randomize()

const fa_gc_header* = ["reads",
                       "gc_content",
                       "gc_bases",
                       "n_bases",
                       "bases",
                       "basename",
                       "absolute_path"].join("\t")

proc random_site(v: VCF, prob_bins: seq[float]): Variant =
    var chrom_select = prob_bins.lowerBound(rand(1.0))
    let chrom_name = v.contigs[chrom_select].name
    let start_pos = rand(v.contigs[chrom_select].length.int)
    let end_pos = v.contigs[chrom_select].length
    let region = fmt"{chrom_name}:{start_pos}-{end_pos}"
    for i in v.query(region):
        return i
    # If nothing found, likely we are too close to the end; Try again
    return random_site(v, prob_bins)

proc print_variant(wtr: VCF, record: Variant) =
    doAssert wtr.write_variant(record)


proc sample*(vcf_fname: string, positions_in: string, n_sites: int) =
    #[
        Calculate GC at various window sizes
    ]#

    # Load all positions
    #var position_set = new_seq[Position]()
    #for pos in positions_in.iter_pos():
    #    position_set.add(pos)
    var vcf: VCF
    doAssert open(vcf, vcf_fname)
    var cum_length = sum(vcf.contigs.mapIt( it.length )).float
    var chrom_name = vcf.contigs.mapIt( it.name )

    
    # Weight chromosomes by length
    var chrom_weights = vcf.contigs.mapIt( it.length.float / cum_length )
    var prob_bins = new_seq[float](chrom_weights.len)
    prob_bins[0] = chrom_weights[0]
    for i in 1..chrom_weights.len - 1:
        prob_bins[i] = prob_bins[i-1] + chrom_weights[i]
    
    var wtr:VCF
    doAssert open(wtr, "/dev/stdout", mode="w")
    wtr.header = vcf.header
    doAssert wtr.write_header()

    # Select a chromosome and position
    for i in 0..<n_sites:
        discard wtr.write_variant(random_site(vcf, prob_bins))

    close(wtr)
    close(vcf)