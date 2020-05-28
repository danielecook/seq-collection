import strformat
import algorithm
import strutils
import sequtils
import os
import hts
import math
import random
import algorithm
import genome_rand
import utils/helpers

randomize()

# proc random_site(v: VCF, prob_bins: seq[float]): Variant =
#     var chrom_select = prob_bins.lowerBound(rand(1.0))
#     let chrom_name = v.contigs[chrom_select].name
#     let start_pos = rand(v.contigs[chrom_select].length.int)
#     let end_pos = v.contigs[chrom_select].length
#     let region = fmt"{chrom_name}:{start_pos}-{end_pos}"
#     for i in v.query(region):
#         return i
#     # If nothing found, likely we are too close to the end; Try again
#     return random_site(v, prob_bins)

proc print_variant(wtr: VCF, record: Variant) =
    doAssert wtr.write_variant(record)

iterator iter_variants(v: VCF, g: genome, var_type: string, n_sites: int): Variant =
    let rng_dist = range_iter("1000")
    var i = 0
    var output: bool
    block site_iter:
        for r in g.random_site(n = 0, rng_dist):
            output = false
            for variant in v.query(r.region):
                output = 
                    (var_type == "all") or
                    (var_type == "snps" and variant.is_snp()) or
                    (var_type == "mnps" and variant.is_mnp()) or
                    (var_type == "indels" and variant.is_indel())
                
                if output == true:
                    i += 1
                    yield variant
                if i >= n_sites:
                    break site_iter

proc sample*(vcf_fname: string, positions_in: string, var_type: string, n_sites: int) =
    #[
        Calculate GC at various window sizes
    ]#

    var wtr:VCF
    doAssert open(wtr, "/dev/stdout", mode="w")
    
    var vcf: VCF
    doAssert open(vcf, vcf_fname)
    wtr.header = vcf.header
    doAssert wtr.write_header()

    var genome = get_genome(vcf)

    for variant in vcf.iter_variants(genome, var_type, n_sites):
        discard wtr.write_variant(variant)


    close(wtr)
    close(vcf)