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

proc calc_window_set(fasta: string, p: Position, windows: seq[int]) =
    # An efficient method for calculating window sizes
    # by loading complete sequence into memory and successively
    # trimming
    # Calculate the GC content for largest window first, then
    # subset, calculate, repeat, etc.
    var f:FAI
    doAssert open(f, fasta)
    var 
        output = new_seq[float](windows.len)
        max_win_len = windows[0]
        initial_sequence: string
    
    try:
        initial_sequence = f.get(p.chrom, p.pos - max_win_len, p.pos + max_win_len)
    except RangeError:
        warning fmt"Variant {p} out of range"
        return

    output[0] = initial_sequence.calc_gc()
    for k, window in windows[1 ..< windows.len]:
        # If we are near the beginning or end of a contig, just
        # fetch from beginning
        if initial_sequence.len == max_win_len * 2 + 1:
            let offset = max_win_len - window
            let sub_seq = initial_sequence[offset .. max_win_len] & initial_sequence[max_win_len + 1 ..< initial_sequence.len - offset]
            output[k + 1] = sub_seq.calc_gc()
        else:
            # If the position lies near the border, query fasta directly
            output[k + 1] = f.get(p.chrom, p.pos - window, p.pos + window).calc_gc()
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

    # Sort windows from max to min
    windows.sort()
    windows.reverse()

    # Parse positions
    for pos in positions_in.iter_pos():
        spawn calc_window_set(fasta, pos, windows)
    sync()

    # open(b, bamfile, index=true)
    # var freq_results = newSeq[FlowVar[chrom_freqs]](b.hdr.targets().len)
    # for idx, contig in b.hdr.targets():
    #     freq_results[idx] = spawn freq_inserts(bamfile, contig.name, verbose)
    # sync()