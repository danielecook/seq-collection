# Generates a random site
import strformat
import hts
import tables
import random
import sequtils
import strformat
import algorithm

randomize()

type
  genome_util = ref object
    chrom_table: Table[string, int]
    cum_length: int
    chrom_weights: seq[float]
    chrom_bins: seq[float]

  site = ref object
    chrom: string
    pos: int

proc `$`(s: site): string =
    return fmt"{s.chrom}:{s.pos}"

#=========#
#   FAI   #
#=========#

proc gen_chrom_table(f: Fai): Table[string, int] =
    # Generate a table of chrom -> chrom_len
    for i in 0..<f.len:
        result[f[i]] = f.chrom_len(f[i])

proc cum_length(f: Fai): int =
    # Calculate the cumulative length of all chromosomes
    result = 0
    for i in 0..<f.len:
        result += f.chrom_len(f[i])

proc chrom_weights(f: Fai): seq[float] =
    # Calculate weight of each chromosome
    return toSeq(0..<f.len).mapIt( f.chrom_len(f[it]) / f.cum_length() )

proc chrom_bins(f: Fai): seq[float] =
    # Generate a sequence of bins based on chrom lengths
    var chrom_weights = f.chrom_weights()
    var prob_bins = new_seq[float](chrom_weights.len)
    prob_bins[0] = chrom_weights[0]
    for i in 1..chrom_weights.len - 1:
        prob_bins[i] = prob_bins[i-1] + chrom_weights[i]
    return prob_bins

proc rand_chrom(g: genome_util): string = 
    var chrom_select = g.chrom_bins.lowerBound(rand(1.0))
    return toSeq(g.chrom_table.keys)[chrom_select]

proc rand_pos(g: genome_util, chrom: string): int =
    return random(g.chrom_table[chrom])

#=========#
#   BAM   #
#=========#

#=========#
#   VCF   #
#=========#

iterator random_site(g: genome_util, n: int): site = 
    for i in 0..<n:
        var chrom = g.rand_chrom()
        let pos = rand_pos(g, chrom)
        yield site(chrom: chrom,
                pos: pos)


proc get_genome_util(f: Fai): genome_util =
    var result = genome_util()
    result.chrom_table = f.gen_chrom_table()
    result.cum_length = f.cum_length()
    result.chrom_weights = f.chrom_weights()
    result.chrom_bins = f.chrom_bins()
    return result


proc genome_rand*(f: Fai) =
    var genome_ref = f.get_genome_util()
    for i in genome_ref.random_site(100):
        echo i

