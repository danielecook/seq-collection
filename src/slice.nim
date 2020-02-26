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
import utils/helpers

proc between[T](val, start, stop: T): bool =
    return val >= start and val <= stop

proc accept_record(r: Record, bed_ln: region_t): bool = 
    #echo r.mate_pos
    if (between(r.start.uint32, bed_ln.start, bed_ln.stop) and
        between(r.stop.uint32, bed_ln.start, bed_ln.stop)) and
        (r.mate_chrom == r.chrom and between(r.mate_pos.uint32, bed_ln.start, bed_ln.stop)):
        # Add logic here to account for mate chrom on another contig but still within bounds
        return true
    return false


proc cmd_slice*(bamfile: string, pos_file: string) =
    #[
        Slice a bam cleanly, returning only reads found within a set of bounds
    ]#
    var 
        b: Bam
        bed_ln: region_t
       
    # Parse genomic positions file
    let stream: Stream =
        if pos_file[^3 .. ^1] == ".gz":
            newGZFileStream(pos_file)
        else:
            newFileStream(pos_file, fmRead)

    open(b, bamfile, index=true)
    for line in lines(stream):
        bed_ln = bed_line_to_region(line)
        for record in b.query(bed_ln.chrom, bed_ln.start.int, bed_ln.stop.int):
            if record.accept_record(bed_ln):
                echo record.tostring()
