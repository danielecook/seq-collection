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
import times

randomize()

proc random_site(vcf_fname: string, prob_bins: seq[float]): Variant = {.gc_safe.}:
    let now = getTime()
    var v: VCF
    doAssert open(v, vcf_fname)
    defer: close(v)
    let chrom_select = prob_bins.lowerBound(rand(1.0))
    let chrom_name = v.contigs[chrom_select].name
    let start_pos = rand(v.contigs[chrom_select].length.int)
    let end_pos = v.contigs[chrom_select].length
    let region = fmt"{chrom_name}:{start_pos}-{end_pos}"
    for i in v.query(region):
        result = i
        #close(v)
        return result
    # If nothing found, likely we are too close to the end; Try again
    return random_site(vcf_fname, prob_bins)

proc print_variant(wtr: VCF, record: Variant) =
    doAssert wtr.write_variant(record)

proc flush(wtr: VCF, variants: seq[FlowVar[Variant]]) =
    for record in variants:
        doAssert wtr.write_variant(^record)

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
    
    # Weight chromosomes by length
    var cum_length = sum(vcf.contigs.mapIt( it.length )).float
    var chrom_weights = vcf.contigs.mapIt( it.length.float / cum_length )
    var prob_bins = new_seq[float](chrom_weights.len)
    prob_bins[0] = chrom_weights[0]
    for i in 1..chrom_weights.len - 1:
        prob_bins[i] = prob_bins[i-1] + chrom_weights[i]
    
    var wtr:VCF
    doAssert open(wtr, "/dev/stdout", mode="w")
    wtr.header = vcf.header
    close(vcf)
    doAssert wtr.write_header()

    # Select a chromosome and position
    var rand_sites = newSeq[FlowVar[Variant]](n_sites)
    var last_out = 0
    const flush_cnt = 100
    setMaxPoolSize(flush_cnt)
    for idx in 0..<n_sites:
        rand_sites[idx] = spawn random_site(vcf_fname, prob_bins)
        if idx %% flush_cnt == 0 and idx > 0:
            sync()
            flush(wtr, rand_sites[idx-flush_cnt+1..<idx])
            last_out = idx
    sync()
    flush(wtr, rand_sites[last_out..<n_sites])
    close(wtr)