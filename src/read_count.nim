import hts
import tables
import strformat

# How to org:
# Use an object to repr. data
# Output obj. as string function
# Set a max depth
# collect records in a seq
# Send to function to tabulate into object (thread)
# fi/fo queue
# Per library output option also...?

type 
    Position* = ref object
        chromosome*: string
        pos*: int64
        base: CountTable[char]
        base_quality: seq[uint8]
        mapping_quality: seq[uint8]
        strand_count: CountTable[bool]
        

proc `$`*(p:Position): string {.inline.} =
    echo fmt"{p.chromosome}:{p.pos} [{p.base}] [{p.strand_count}]"


proc cmd_read_count*(bamfile: string, positions: string) =
    #[
        Calculates insert size
    ]#
    var 
        b: Bam
        offset: int
        result: Position

    open(b, bamfile, index=true)
    for target in 999915..1000000:
        var read_count  = Position()
        for record in b.query("I", target, target + 1):
                offset = target - record.start.int
                var base_quality = record.base_quality_at(offset)
                read_count.chromosome = record.chrom
                read_count.pos = record.start
                read_count.base.inc record.base_at(offset)
                read_count.base_quality.add base_quality
                read_count.mapping_quality.add record.mapping_quality
                read_count.strand_count.inc (record.flag().reverse() == false)
                echo record.cigar
                #echo record.flag().has_flag(BAM_FREVERSE)
                #echo record.has_flag(16.uint16)
                #echo target, ": ", record.base_at(offset), " → ", record.base_quality_at(offset), " → ", offset
        echo read_count

# avg_se_mapping_quality → mean single ended mapping quality
# [x] num_plus_strand → number of reads on the plus/forward strand
# [x] num_minus_strand → number of reads on the minus/reverse strand
# avg_pos_as_fraction → average position on the read as a fraction (calculated with respect to the length after clipping). This value is normalized to the center of the read (bases occurring strictly at the center of the read have a value of 1, those occurring strictly at the ends should approach a value of 0)
# avg_num_mismatches_as_fraction → average number of mismatches on these reads per base
# avg_sum_mismatch_qualities → average sum of the base qualities of mismatches in the reads
# num_q2_containing_reads → number of reads with q2 runs at the 3’ end
# avg_distance_to_q2_start_in_q2_reads → average distance of position (as fraction of unclipped read length) to the start of the q2 run
# avg_clipped_length → average clipped read length of reads
# avg_distance_to_effective_3p_end → average distance to the 3’ prime end of the read (as fraction of unclipped read length)