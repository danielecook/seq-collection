import hts
import streams
import zip/gzipfiles
import strutils
import deques
import strformat
import tables

# Initialize a deque object
# for adjacent positions
type Position = ref object
    chrom: string
    pos: int

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

proc adjacent_pos(positions: Deque[Position], limit: int): seq[Position] =
    var p1 = positions[0]
    var p2 = positions[1]
    if p1.chrom == p2.chrom and abs(p2.pos - p1.pos) <= limit and  abs(p1.pos - p2.pos) != 0:
        return @[p1, p2]
    return new_seq[Position](0)

proc haplotype_clean(hap: seq[char]): bool =
    for i in hap:
        if i in ['A', 'T', 'C', 'G'] == false:
            return false
    return true

proc collect_reads(bam: Bam, chrom: string, targets: seq[int]): Table[string, Record] =
    # Collect the set of reads overlapping 
    var output: Table[string, Record]
    var read: string
    for pos in targets:
        for record in bam.query(chrom, pos, pos + 1):
            if record.accept_record():
                read = if record.flag.read1: "+1" 
                       else: "+2"
                   
            output[record.qname & read] = record.copy()
    return output



proc cmd_contamination*(bamfile: string, pos_file: string) =
    #[
        Calculates insert size
    ]#
    var 
        b: Bam
        chrom: string
        pos: int
        offset_1: int
        offset_2: int
        line_in: seq[string]
        target_1: Position
        target_2: Position
        haplotype: seq[char]
        
    

    # Parse genomic positions file
    let stream: Stream =
        if pos_file[^3 .. ^1] == ".gz":
            newGZFileStream(pos_file)
        else:
            newFileStream(pos_file, fmRead)

    open(b, bamfile, index=true)
    for line in lines(stream):
        var haps: CountTable[string]
        if line.startsWith("#"):
            continue
        line_in = line.split("\t")
        chrom = line_in[0]
        pos = line_in[1].parseInt - 1
        var pos_in = Position(chrom: chrom, pos: pos)
        pos_set.addFirst pos_in

        if len(pos_set) > 2:
            pos_set.popLast()

        if len(pos_set) == 2:
            var adj = adjacent_pos(pos_set, 1000)
            if adj.len == 2:
                target_1 = adj[1]
                target_2 = adj[0]
                # Iterate both targets, collecting the complete set of reads for each position
                var read_set = b.collect_reads(target_1.chrom, @[target_1.pos, target_2.pos])
                # Restructure read set into hash table
                for key, record in read_set.pairs:
                    # See https://github.com/ncbi/ncbi-vdb/blob/master/interfaces/insdc/insdc.h#L74
                    # =ACM ...; '=' is missing.
                    #echo target_1.pos
                    if record.flag.read1:
                        offset_1 = target_1.pos.int - record.start.int
                        offset_2 = target_2.pos.int - record.start.int
                        if offset_1 > 0:
                            var b1 = record.base_at(offset_1)
                            var b2 = record.base_at(offset_2)
                            var v1 = record.base_quality_at(offset_1)
                            var v2 = record.base_quality_at(offset_2)
                            haplotype = @[record.base_at(offset_1), record.base_at(offset_2)]
                            echo record.chrom, ":", target_1.pos, " ", record.start, " ", offset_1, " ", offset_2, " ", haplotype, " ", @[v1, v2], " ", record.qname, "-->", record.flag.read1
                            var s: string
                            record.sequence(s)
                            echo s
                    # var s: string
                    # echo record.sequence(s)
                    # echo s
                    #echo record.qname
                    #for i in haplotype:
                        #if i in ['A', 'T', 'C', 'G'] == false:
                            # Something to get mate here?
                    # if haplotype.haplotype_clean():
                    #     haps.inc haplotype.join()
        # if haps.len > 0:
        #     echo target_2.pos - target_1.pos
        #     echo haps, ":", $target_1.pos, " + ", $target_2.pos