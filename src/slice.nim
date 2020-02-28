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


proc in_region(chrom: string, pos: uint32, region_set: Table[string, seq[region_t]]): bool =
    # TODO: speed this up with binary search (lapper library from brentp)
    if region_set.hasKey(chrom):
        for region in region_set[chrom]:
            if between(pos, region.start, region.stop):
                return true
    return false

proc accept_record(r: Record, region_set: Table[string, seq[region_t]]): bool = 
    #echo r.mate_pos
    if (in_region(r.chrom, r.start.uint32, region_set) and
        in_region(r.chrom, r.stop.uint32, region_set) and
        in_region(r.mate_chrom, r.mate_pos.uint32, region_set)):
        #(r.mate_chrom == r.chrom and between(r.mate_pos.uint32, bed_ln.start, bed_ln.stop)):
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
        region_set: Table[string, seq[region_t]]
       
    # Parse genomic positions file
    let stream: Stream =
        if pos_file[^3 .. ^1] == ".gz":
            newGZFileStream(pos_file)
        else:
            newFileStream(pos_file, fmRead)
    
    # Read bed file into array
    for line in lines(stream):
        bed_ln = bed_line_to_region(line)
        if (bed_ln.chrom in region_set) == false:
            region_set[bed_ln.chrom] = new_seq[region_t](0)
        region_set[bed_ln.chrom].add bed_ln

    stream.setPosition(0)

    open(b, bamfile, index=true)
    for line in lines(stream):
        bed_ln = bed_line_to_region(line)
        for record in b.query(bed_ln.chrom, bed_ln.start.int, bed_ln.stop.int):
            if record.accept_record(region_set):
                echo record.tostring()
