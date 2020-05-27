# Generates a random site
import strformat
import hts
import tables
import random
import sequtils
import strformat
import algorithm
import parseutils
import utils/bed

randomize()

type
    genome = ref object
        chrom_table: Table[string, Region]
        cum_length: int
        chrom_weights: seq[float]
        chrom_bins: seq[float]
  
    site = ref object
        chrom: string
        pos: int
        one: int
  
proc `$`(s: site): string =
    return fmt"{s.chrom}:{s.pos+s.one}"

#=========#
#   FAI   #
#=========#

proc gen_chrom_table(f: Fai, bed: string): Table[string, Region] =
    # Generate a table of chrom -> chrom_len
    # For bed files this is the sum of regions on that chrom
    if bed != "":
        for region in bed.iter_bed():
            result[$region] = region
    else:
        for i in 0..<f.len:
            var reg = Region(chrom: f[i],
                                  start: 0,
                                  stop: f.chrom_len(f[i]))
            result[$reg] = reg

#=========#
#   BAM   #
#=========#

#=========#
#   VCF   #
#=========#

#==================#
#   genome table   #
#==================#

proc cum_length(chr_tbl: Table[string, Region]): int =
    # Calculate the cumulative length of all chromosomes/regions
    return toSeq(chr_tbl.values()).mapIt( it.len ).foldl( a + b )

proc chrom_weights(chr_tbl: Table[string, Region]): seq[float] =
    # Calculate weight of each chromosome/region
    let cum_length = chr_tbl.cum_length()
    return toSeq(chr_tbl.values()).mapIt( it.len / cum_length )

proc chrom_bins(chr_tbl: Table[string, Region]): seq[float] =
    # Generate a sequence of bins based on chrom lengths
    var chrom_weights = chr_tbl.chrom_weights()
    var prob_bins = new_seq[float](chrom_weights.len)
    prob_bins[0] = chrom_weights[0]
    for i in 1..chrom_weights.len - 1:
        prob_bins[i] = prob_bins[i-1] + chrom_weights[i]
    return prob_bins

proc rand_region(g: genome): Region = 
    let region_select = g.chrom_bins.lowerBound(rand(1.0))
    return toSeq(g.chrom_table.values)[region_select]

proc rand_pos(g: genome, region: Region): int =
    let r = g.chrom_table[$region]
    return random(r.len) + r.start


iterator random_site(g: genome, n: int, one: int): site = 
    for i in 0..<n:
        var region = g.rand_region()
        let pos = rand_pos(g, region)
        yield site(chrom: region.chrom,
                   pos: pos,
                   one: one)

proc get_genome(f: Fai, bed: string): genome =
    var result = genome()
    let chrom_table = f.gen_chrom_table(bed)
    result.chrom_table = chrom_table
    result.cum_length = chrom_table.cum_length()
    result.chrom_weights = chrom_table.chrom_weights()
    result.chrom_bins = chrom_table.chrom_bins()
    return result

proc genome_rand*(f: Fai, n_sites: int, bed: string, one: int) =
    var genome_ref = f.get_genome(bed)
    for i in genome_ref.random_site(n_sites, one):
        if i.pos > 0:
            echo i, "\t", f.get(i.chrom, i.pos, i.pos)
        else:
            echo i, "\t", f.get(i.chrom, i.pos, i.pos+1)[0]
