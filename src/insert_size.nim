import tables
import hts
import sequtils
import strutils
import math
import stats
import os
import re
import threadpool
#import ggplotnim

const header* = ["median",
                 "mean",
                 "min",
                 "max_all",
                 "n_reads",
                 "n_accept",
                 "n_use",
                 "sample",
                 "filename"].join("\t")

const distribution_header = ["insert_size", "count", "bam"].join("\t")

proc accept_record(r: Record): bool = 
    if (r.flag.pair == false or
        r.flag.unmapped == true or
        r.flag.mate_unmapped or
        r.flag.read1 == true or
        r.flag.secondary or
        r.flag.supplementary or
        r.flag.dup or
        r.isize == 0):
        return false
    return true

proc mean_freq(arr: seq[int64]): float64 =
    var total = 0i64
    for idx, val in arr:
        total += (idx+1)*val
    return total.float / sum(arr).float

proc median_freq(arr: seq[int64]): int64 = 
    # Calculate the median index of a frequency distribution
    var total = sum(arr)
    var running_total = 0i64
    for idx, val in arr:
        running_total += val
        if running_total.float / total.float >= 0.5:
            return idx

proc bam_sample(b: Bam): string =
    # Fetches the sample name from the header
    for line in splitLines($(b.hdr)):
        if line.startsWith("@RG"):
            for field in line.split("\t"):
                if field.startswith("SM:"):
                    return field.replace("SM:", "")
    return ""

proc freq_inserts(bam_name: string, contig: string, contig_length: uint32): int =
    var b: Bam
    var n = 0
    open(b, bam_name, index=true)
    for record in b.query(contig, start=1, stop=contig_length.int):
        if record.accept_record():
            n += 1
    return n

proc cmd_insert_size*(bamfile: string, distfile: string, threads: int8) =
    # [ ] TODO: Reimplement this with a count table?
    #  create procs for MAD (median abs. dev)
    #  and use to calculate width...?
    #  proc for getting median also?
    #  Proc for freq proc.
    #  https://nimble.directory/pkg/rbtree
    # https://nimble.directory/pkg/binaryheap

    # https://leetcode.com/problems/find-median-from-data-stream/solution/
    var b: Bam
    let fname = bam_file.lastPathPart().changeFileExt("")
    let ins_arr = 10000
    var inserts = new_seq[int64](ins_arr)
    var inserts_trimmed = new_seq[int64](ins_arr)
    var overflow: seq[int64]
    var n_reads = 0
    var n_accept = 0

    # Option 1: 0m2.951s parallelization
    open(b, bamfile, threads=threads, index=true)
    # var freq_results = newSeq[FlowVar[int]](b.hdr.targets().len)
    # for idx, contig in b.hdr.targets():
    #     echo contig.name
    #     freq_results[idx] = spawn freq_inserts(bamfile, contig.name, contig.length)
    #     echo contig.name, " ", contig.length
    # echo "G"
    # sync()
    # for idx, contig in b.hdr.targets():
    #     echo ^freq_results[idx], " - ", contig
    #     echo "---"

    # Can potentially parallelize here across contigs
    for record in b:
        n_reads += 1
        if record.accept_record():
            n_accept += 1
            var insert_val = abs(record.isize)
            if insert_val <= inserts.len:
                inserts[insert_val-1] += 1
            else:
                overflow.add(insert_val)

    # Now calculate the median insert size
    # Add an option to include 'overflow length'
    var total_length = sum(inserts) + overflow.len
    
    # Identify midpoint index
    var running_total = 0.int64
    var cumulative_sum: seq[int64]
    
    # Trim last 0.5%
    for idx, val in inserts:
        running_total += val
        cumulative_sum.add(running_total)
        if running_total.float / total_length.float <= 0.995:
            inserts_trimmed[idx] = val
    
    #echo inserts_trimmed
    let median_insert_size = median_freq(inserts_trimmed)
    let mean_insert_size = mean_freq(inserts_trimmed)
    let min_insert_size = inserts_trimmed.filterIt(it > 0)[0]
    
    # # Plotting
    # let x = toSeq(1..inserts_trimmed.len)
    # #echo x
    # let y = inserts_trimmed
    # #echo inserts_trimmed
    # var df = seqsToDf(x, y).filter( f{ "y" > 0 } )
    # echo df
    # ggplot(df, aes(x = "x", y = "y")) +
    #     geom_point() +
    #     ggtitle("ggplotnim - or I Suck At Naming Things™") +
    #     ggsave("scatterColor.pdf")

    # Output distribution
    if distfile != "":
        var f = open(distfile, fmWrite)
        f.writeLine(distribution_header)
        for idx, val in inserts_trimmed.filterIt(it > 0):
            f.writeLine([$idx, $val, fname].join("\t"))
        f.close()

    var result = [$median_insert_size,
                  $math.round(mean_insert_size, 3),
                  $min_insert_size,
                  $max(overflow),
                  $n_reads,
                  $n_accept,
                  $(sum(inserts_trimmed)),
                  bam_sample(b),
                  fname]
    echo result.join("\t")
    # var mad = new_seq[int64](cumulative_sum.len)
    # for idx, val in inserts:
    #     echo inserts[median_insert_size], "<- mid freq"
    #     mad[idx] = abs(val - inserts[median_insert_size])
    # var mad_val = median_freq(mad)*deviations

    # # set inserts to 0 if outside mad
    # for idx, val in inserts:
    #     if val <= (mad_val + median_insert_size):
    #         inserts[idx] = 0
    # echo inserts
    # echo (sum(inserts).float/inserts.len.float)
    # # Now get MAD

    # var percentile = new_seq[float](ins_arr)
    # var max_insert_size = 0
    # for record in b:
    #    if record.accept_record():
    #         var idx = abs(record.isize)
    #         if idx <= inserts.len:
    #             inserts[idx-1] += 1
    #         else:
    #             if idx > max_insert_size:
    #                 max_insert_size = idx
    
    # var total_length = 0.int64
    # var trimmed_count = 0.int64
    # var trimmed_length = 0.int64
    # for idx, val in inserts:
    #     total_length += val.int64*(idx+1)

    # # Trim 1 pct from end.
    # var running_total = 0i64
    # var max_idx = 0
    # for idx, val in inserts:
    #     var length_add = val.int64*(idx+1)
    #     running_total += length_add
    #     var pct = running_total.float / total_length.float
    #     if pct <= 0.995:
    #         percentile[idx] = pct
    #         trimmed_count += val.int64
    #         trimmed_length += length_add
    #         insert_lengths[idx] = length_add
    #         max_idx = idx
    #     else:
    #         percentile[idx] = 1.0

    # var median_insert_size: int64
    # for idx, val in percentile:
    #     if val <= 0.50:
    #         median_insert_size = idx + 1
    #     else:
    #         break
    # echo median_insert_size, "<- median"
    # # Stats
    # let mode_insert_size = inserts.find(max(inserts))+1       
    # let min_insert_size = inserts.find(inserts.filterIt( it > 0 )[0]) + 1
    # let mean_insert_size = (trimmed_length.float / trimmed_count.float)

    # # Standard dev
    # let n = sum(inserts[1..max_idx])
    # var m = 0.int64
    # for idx, val in inserts[1..max_idx]:
    #     m += val * (idx + 1)^2
    # let variance = (m.float - (n.float*(mean_insert_size^2))) / (n - 1).float
    # let std_dev = pow(variance, 0.5)

    # var result = [
    #     $min_insert_size,
    #     $mean_insert_size,
    #     $median_insert_size,
    #     $mode_insert_size,
    #     $std_dev,
    #     $max_insert_size,
    #     bamfile
    # ]
    
    # echo $result


