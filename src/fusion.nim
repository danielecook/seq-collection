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

var logger = newConsoleLogger()
addHandler(logger)
setLogFilter(lvlInfo)

type SplitRead = ref object
    chrom: string
    pos: int64
    strand: char
    cigar: string
    mapq: int
    edit_distance: int
    split_size: int64

type CompoundRecord = ref object
    # The compound record combines a read with its partner and split reads
    read1: Record
    read2: Record
    split1: SplitRead
    split2: SplitRead

proc strand(r: Record): char =
    return if r.flag.reverse: '-' else: '+'

proc `$`(sr: SplitRead): string =
    return fmt"{sr.chrom}:{sr.pos} -- {sr.cigar} -- SS:{sr.split_size}"

proc `$`(cr: CompoundRecord): string =
    var output = ""
    if cr.read1 != nil:
        output = output & fmt"R1={cr.read1.chrom}:{cr.read1.start}-{cr.read1.stop} ({cr.read1.strand})  "
    if cr.split1 != nil:
        output = output & fmt"S1={cr.split1.chrom}:{cr.split1.pos} ({cr.split1.strand}) {cr.split1.split_size}  "
    return output

    #return fmt"{sr.chrom}:{sr.pos} -- {sr.cigar}"

proc getSplit(record: Record): SplitRead =
    var split_read = tag[string](record, "SA")
    if split_read.isSome:
        var sr = SplitRead()
        # Only consider single splits; Multiple
        # adds too much complexity
        var splitset = split_read.get().split(";")[0].split(",")
        sr.chrom = splitset[0]
        sr.pos = splitset[1].parseInt()
        sr.strand = if splitset[2] == "-": '-' else: '+'
        sr.cigar = splitset[3]
        sr.mapq = splitset[4].parseInt()
        sr.edit_distance = splitset[5].parseInt()
        if sr.chrom == record.chrom:
            sr.split_size = min(abs(record.start - sr.pos), abs(record.stop - sr.pos))
        else:
            sr.split_size = 0 
        return sr
    return nil

proc process(rec: Record): CompoundRecord =
    var cr = CompoundRecord()
    cr.read1 = rec.copy()
    cr.split1 = cr.read1.getSplit()

    
    return nil

proc accept(r: Record): bool = 
    # Extract chimeric reads
    # Paired read with large insert size or existing on different chromosomes
    if (r.flag.pair and
        r.flag.dup == false and
        (
            # Insert size == 0: diff contigs
            (abs(r.isize()) > 10000) or
            (r.isize() == 0 and r.chrom != r.mate_chrom)
        )
        ):
        return true
    # Test if split hit read
    var sr = r.getSplit()
    if sr != nil and (sr.split_size > 10000 or sr.split_size == 0):
        return true
    return false


proc cmd_fusion*(bamfile: string) =
    #[
        Calculates insert size
    ]#
    var 
        b: Bam
       
    open(b, bamfile, index=true)
    for record in b:
        if record.accept():
            echo record.tostring()
            #if record.getSplit() != nil:
            #    echo record.getSplit()
        #var cr = record.process()
        #if $cr != "":
        #    echo $cr
        #if record.accept():
        #    echo $record
        # var split_read = tag[string](record, "SA")
        # if split_read.isSome:
        #     echo $split_read, " ", $record.isize()
        #if record.isize() == 0 or abs(record.isize) >= 10000:
        #    echo $record.tostring()
        #if record.accept():
            #echo $record, " accept"
            #echo $record.isize()