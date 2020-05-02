import memfiles
import streams
import zip/gzipfiles
import utils/helpers
import strformat
import algorithm
import strutils
import sequtils
import os
import hts
import math
import threadpool

const fa_gc_header* = ["reads",
                          "gc_content",
                          "gc_bases",
                          "n_bases",
                          "bases",
                          "basename",
                          "absolute_path"].join("\t")

# Result container
type
  result_gc = ref object
    chr_pos: Position
    windows: seq[float]

proc calc_gc(sequence: string): float =
    return (sequence.count({'G', 'C', 'g', 'c'}) / sequence.count({'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}))

proc sub_seq(sequence: string, start_pos: int, end_pos: int): string =
    var left, right: int
    left = start_pos
    right = end_pos
    if left <= 0 or left > right:
        left = 0
    if right >= (sequence.len - 1):
        right = sequence.len - 1
    return sequence[left .. right]

proc calc_window_set(fasta: string, p: Position, windows: seq[int]): result_gc =
    # An efficient method for calculating window sizes
    # by loading complete sequence into memory and successively
    # trimming
    # Calculate the GC content for largest window first, then
    # subset, calculate, repeat, etc.
    var 
        output = new_seq[float](windows.len)
        trim_len: int

    # Convert from 1-based coordinates
    if p.pos0 > fasta.len - 1:
        warning_msg fmt"{p} is out of range"
        return

    for k, window in windows:
        output[k] = fasta.sub_seq(p.pos0 - window, p.pos0 + window).calc_gc().round(window.intToStr().len + 2)
    
    return result_gc(chr_pos: p, windows: output)


proc fa_gc*(fasta: string, positions_in: string, windows_in: seq[string]) =
    #[
        Calculate GC at various window sizes
    ]#
    var f:FAI
    doAssert open(f, fasta)

    # Parse windows
    var windows = new_seq[int](windows_in.len)
    for i in 0..<windows_in.len:
        windows[i] = (windows_in[i].sci_parse_int().float).int
        if windows[i] < 1:
            quit_error "Window lengths must be >= 1"

    # Load all positions
    var position_set = new_seq[Position]()
    for pos in positions_in.iter_pos():
        position_set.add(pos)
    
    position_set.sort(helpers.genome_cmp)

    # Print header
    echo (@["chrom", "pos"] & windows.mapIt( fmt"gc_{it * 2}" )).join("\t")

    var max_pos: int
    var curr_chrom: string
    var chrom_fasta: string
    var thread_count = 0
    var fa_gc_results = newSeq[FlowVar[result_gc]](position_set.len)
    for idx, chr_pos in position_set:
        # load chromosome
        if curr_chrom != chr_pos.chrom:
            sync()
            max_pos = position_set.filterIt( it.chrom == chr_pos.chrom ).mapIt( it.pos0 ).max()
            chrom_fasta = f.get(chr_pos.chrom, 0, max_pos + windows.max())
            curr_chrom = chr_pos.chrom
        fa_gc_results[idx] = spawn calc_window_set(chrom_fasta, chr_pos, windows)
    sync()

    for idx in 0..<fa_gc_results.len:
        var result = ^fa_gc_results[idx]
        echo [result.chr_pos.chrom,
              $result.chr_pos.pos].join("\t") & "\t" & result.windows.join("\t")