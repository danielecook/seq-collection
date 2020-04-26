import hts
import streams
import zip/gzipfiles
import strutils
import deques
import strformat
import tables
import sequtils
import logging

var logger = newConsoleLogger()
addHandler(logger)
setLogFilter(lvlInfo)

# Initialize a deque object
# for adjacent positions
type Position = ref object
    chrom: string
    pos: int

type Stats = ref object
    base_mismatch: int
    paired_read_hap: int
    single_read_hap: int

proc `$`(s: Stats): string = 
    return fmt"base_mismatch={s.base_mismatch} pe_hap:{s.paired_read_hap} se_hap:{s.single_read_hap}"

proc `$`(p: Position): string = 
    return fmt"{p.chrom}:{p.pos}"

var pos_set = initDeque[Position](2)

proc accept_record(r: Record): bool = 
    # Note this function differs from
    # the one in insert_size.nim
    if (r.flag.pair == false or
        r.flag.unmapped == true or
        r.flag.mate_unmapped or
        r.flag.secondary or
        r.flag.supplementary or
        r.flag.dup):
        return false
    return true

proc adjacent_pos(positions: Deque[Position], proximity: int): seq[Position] =
    # Collects positions on the same chromosome within 
    # proximity bp of one another.
    if positions.len == 2:
        var p1 = positions[0]
        var p2 = positions[1]
        if p1.chrom == p2.chrom and abs(p2.pos - p1.pos) <= proximity and abs(p1.pos - p2.pos) != 0:
            return @[p1, p2]
    return new_seq[Position](0)

proc get_target_bases(record: Record, targets: seq[Position], min_base_qual: uint8): seq[char] =
    var offset_1 = targets[0].pos.int - record.start.int
    var offset_2 = targets[1].pos.int - record.start.int
    var q1 = if offset_1 >= 0: record.base_quality_at(offset_1) else: 0.uint8
    var q2 = if offset_1 >= 0: record.base_quality_at(offset_2) else: 0.uint8
    var a1 = if q1 > min_base_qual: record.base_at(offset_1) else: '.'
    var a2 = if q2 > min_base_qual: record.base_at(offset_2) else: '.'
    return @[a1, a2] 

proc get_read2_target_bases(read: Record, targets: seq[Position], read_set: Table[string, Record], min_base_qual: uint8): seq[char] =
    var read2qname = read.qname & "+2"
    if read2qname in read_set:
        var r2 = read_set[read2qname]
        return r2.get_target_bases(targets, min_base_qual)


proc haplotype_complete(s: seq[char]): bool =
    return s.mapIt( (it != '.').int ).foldl(a + b) == 2

proc haplotype_mismatch(a: seq[char], b: seq[char]): bool =
    # Check two haplotypes and if there are base call differences return true
    # [A, T] == [A, T] → false
    # [., T] == [A, .] → false
    # [A, T] == [A, C] → true
    for i in 0..<2:
        if a[i] != b[i] and a[i] != '.' and b[i] != '.':
            return true
    return false

proc fill_haplotype(a: seq[char], b: seq[char]): seq[char] =
    # Combine two haplotypes
    # [A, .] + [., T] → [A, T]
    var hap: seq[char]
    newSeq(hap, 2)
    for i in 0..<2:
        hap[i] = if a[i] != '.': a[i] else: b[i]
    return hap

proc collect_reads(bam: Bam, targets: seq[Position]): Table[string, Record] =
    var chrom = targets[0].chrom
    # Collect the set of reads overlapping 
    var output: Table[string, Record]
    var read: string
    for target in targets:
        for record in bam.query(chrom, target.pos, target.pos + 1):
            if record.accept_record():
                read = if record.flag.read1: "+1" 
                       else: "+2"
                   
            output[record.qname & read] = record.copy()
    return output


iterator haplotypes(b: Bam, targets: seq[Position], stats: ptr Stats, min_base_qual: uint8): seq[char] =
    # Identify haplotypes for single and paired-end reads
    # Yield base haplotypes
    var single_read_hap_set: string
    var read_set = b.collect_reads(targets)

    for key, read in read_set.pairs:
        # In this loop, read can be either read 1 or 2
        #
        # It will only consider #3 (below), when it is the first read in a pair.
        # (#3) only occurs when haplotype is not filled
        # and only when starting with the first read.
        #
        # If a read has both variants, it's pair will be checked for mismatches
        # and discarded if any are present. 
        #
        # Target arragenements
        # 
        # (1) All in first read
        # ----x----x--->   <-------------
        #
        # (2) All in second read
        # ------------->   <---x---x-----
        #                   
        # (3) Across both reads
        # ----x-------->   <------x------
        #
        # (4) Mix
        # ----x-------x->  
        #      <------x-----
        #
        # (5) Both → Ignored as this double counts
        # ----x-------x->  
        # <---x-------x--
        #
        # See https://github.com/ncbi/ncbi-vdb/blob/master/interfaces/insdc/insdc.h#L74
        var bases, read_bases, r2_bases: seq[char]
        newSeq(bases, 2)

        if read.qname in single_read_hap_set:
            # (5): If this is read 2, and read 1 had a complete haplotype
            # we can keep going, as this has already been counted
            debug(fmt"Prev. seen read: {read.get_target_bases(targets, min_base_qual)} -- {read.get_read2_target_bases(targets, read_set, min_base_qual)}")
            continue

        read_bases = read.get_target_bases(targets, min_base_qual)
        
        if read.flag.read1:
            # (1, 4, 5): Check second read for mismatches
            if read_bases.haplotype_complete():
                r2_bases = read.get_read2_target_bases(targets, read_set, min_base_qual)
                if r2_bases.len == 2:
                    if read_bases.haplotype_mismatch(r2_bases):
                        debug(fmt"SE mismatch {read_bases} {r2_bases}")
                        stats.base_mismatch += 1
                        continue
                # Single read haplotype found
                debug(fmt"SE R1-haplo {read_bases} {r2_bases}")
                stats.single_read_hap += 1
                single_read_hap_set.add read.qname
                bases = read_bases
                yield bases

            # (3): Check second read 
            elif read_bases.haplotype_complete() == false:
                r2_bases = read.get_read2_target_bases(targets, read_set, min_base_qual)
                if r2_bases.len == 2:
                    if read_bases.haplotype_mismatch(r2_bases):
                        stats.base_mismatch += 1
                        debug(fmt"PE mismatch {read_bases} {r2_bases}")
                        continue

                    # Fill bases with calls
                    bases = read_bases.fill_haplotype(r2_bases)
                    if bases.haplotype_complete():
                        # Paired read haplotype found
                        stats.paired_read_hap += 1
                        debug(fmt"PE haplo    {read_bases} {r2_bases}")
                        yield bases
        else:
            # (2): If both targets are in read 2, we can skip checking read 1.
            if read_bases.haplotype_complete():
                stats.single_read_hap += 1
                debug(fmt"SE R2-haplo {read_bases} {r2_bases}")
                bases = read_bases
                yield bases


proc cmd_contamination*(bamfile: string, pos_file: string) =
    #[
        Calculates insert size
    ]#
    var 
        b: Bam
        chrom: string
        pos: int
        line_in: seq[string]
        min_base_qual = 10.uint8
        hap_counter: CountTable[seq[char]]
        hap_stats = Stats()
       

    # Parse genomic positions file
    let stream: Stream =
        if pos_file[^3 .. ^1] == ".gz":
            newGZFileStream(pos_file)
        else:
            newFileStream(pos_file, fmRead)

    open(b, bamfile, index=true)
    for line in lines(stream):
        if line.startsWith("#"):
            continue
        line_in = line.split("\t")
        chrom = line_in[0]
        pos = line_in[1].parseInt
        pos_set.addLast Position(chrom: chrom, pos: pos)

        if len(pos_set) > 2:
            pos_set.popFirst()

        # adjacent positions checks that 2 positions are available.
        var targets = adjacent_pos(pos_set, 1000)
        if targets.len < 2:
            continue
        
        hap_counter.clear()
        for haplotype in b.haplotypes(targets, hap_stats.addr, min_base_qual):
            hap_counter.inc haplotype
        if hap_counter.len > 2:
            echo fmt"{targets} {hap_counter} {hap_stats}"