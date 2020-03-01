import hts
import streams
import zip/gzipfiles
import strutils
import deques
import strformat
import tables
import sequtils
import logging
import math
import lapper
import utils/helpers

proc between[T](val, start, stop: T): bool =
    return val >= start and val <= stop

proc in_region(chrom: string, pos: int, ls: lapper_set): bool =
    var empty:seq[region]
    if ls.hasKey(chrom):
        var l = ls[chrom]
        return l.find(pos, pos, empty)
    return false

proc accept_record(r: Record, ls: lapper_set): bool = 
    #echo r.mate_pos
    if (in_region(r.chrom, r.start.int, ls) and
        in_region(r.chrom, r.stop.int, ls) and
        in_region(r.mate_chrom, r.mate_pos.int, ls)):
        #(r.mate_chrom == r.chrom and between(r.mate_pos.uint32, bed_ln.start, bed_ln.stop)):
        # Add logic here to account for mate chrom on another contig but still within bounds
        return true
    return false


proc cmd_slice*(bamfile: string, bed_file: string) =
    #[
        Slice a bam cleanly, returning only reads found within a set of bounds
    ]#
    var 
        b: Bam
        bed_ln: region
       
    var ls = bedfile_to_lapper(bed_file)

    let stream: Stream =
        if bed_file[^3 .. ^1] == ".gz":
            newGZFileStream(bed_file)
        else:
            newFileStream(bed_file, fmRead)

    # Now iterate through file again for efficiency.
    open(b, bamfile, index=true)
    for line in lines(stream):
        bed_ln = bed_line_to_region(line)
        for record in b.query(bed_ln.chrom, bed_ln.start.int, bed_ln.stop.int):
            if record.accept_record(ls):
                echo record.tostring()
