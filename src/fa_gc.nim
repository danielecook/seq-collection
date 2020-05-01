import memfiles
import streams
import zip/gzipfiles
import utils/helpers
import strformat
import algorithm
import strutils
import os
import hts
import threadpool

const fa_gc_header* = ["reads",
                          "gc_content",
                          "gc_bases",
                          "n_bases",
                          "bases",
                          "basename",
                          "absolute_path"].join("\t")

proc calc_gc(sequence: string): float =
    return (sequence.count({'G', 'C', 'g', 'c'}) / sequence.count({'A', 'T', 'C', 'G', 'a', 't', 'c', 'g'}))

proc sub_seq(sequence: string, start_pos: int, end_pos: int): string =
    var left, right: int
    left = start_pos
    right = end_pos
    if left <= 0:
        left = 0
    if right >= (sequence.len - 1):
        right = sequence.len - 1
    return sequence[left .. right]

proc calc_window_set(fasta: string, p: Position, windows: seq[int]) =
    # An efficient method for calculating window sizes
    # by loading complete sequence into memory and successively
    # trimming
    # Calculate the GC content for largest window first, then
    # subset, calculate, repeat, etc.
    var 
        output = new_seq[float](windows.len)

    for k, window in windows:
        output[k] = fasta.sub_seq(p.pos - window, p.pos + window).calc_gc()
    echo fmt"{p.chrom}:{p.pos} {output}"


proc fa_gc*(fasta: string, positions_in: string, windows_in: seq[string]) =
    #[
        Calculate GC at various window sizes
    ]#
    var f:FAI
    doAssert open(f, fasta)

    # Parse windows
    var windows = new_seq[int](windows_in.len)
    for i in 0..<windows_in.len:
        windows[i] = windows_in[i].parseInt()
        if windows[i] < 1:
            quit_error "Window lengths must be >= 1"

    # Load all positions
    var position_set = new_seq[Position]()
    for pos in positions_in.iter_pos():
        position_set.add(pos)
    
    position_set.sort(helpers.genome_cmp)

    var curr_chrom: string
    var chrom_fasta: string
    for chr_pos in position_set:
        # load chromosome
        if curr_chrom != chr_pos.chrom:
            warning_msg fmt"Loading chrom {chr_pos.chrom}"
            chrom_fasta = f.get(chr_pos.chrom)
            curr_chrom = chr_pos.chrom
        calc_window_set(chrom_fasta, chr_pos, windows)
    #sync()

    # open(b, bamfile, index=true)
    # var freq_results = newSeq[FlowVar[chrom_freqs]](b.hdr.targets().len)
    # for idx, contig in b.hdr.targets():
    #     freq_results[idx] = spawn freq_inserts(bamfile, contig.name, verbose)
    # sync()