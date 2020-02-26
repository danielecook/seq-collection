import hts
import strutils
import strformat
import math
import os
import colorize
import algorithm
import threadpool
import ggplotnim
import sequtils
import utils/helpers

const INS_ARR = 10000

const insert_size_header* = ["median",
                 "mean",
                 "std_dev",
                 "min",
                 "percentile_99.5",
                 "max_all",
                 "n_reads",
                 "n_accept",
                 "n_use",
                 "sample"].join("\t")

const distribution_header = ["insert_size", "count", "sample"].join("\t")

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

type
  chrom_freqs = ref object
    inserts, overflow: seq[int64]
    n_reads, n_accept: int

proc freq_inserts(bamfile: string, contig: string, verbose: bool): chrom_freqs =
    var b: Bam
    var n_reads = 0
    var n_accept = 0
    var inserts = new_seq[int64](INS_ARR)
    var overflow: seq[int64]
    open(b, bamfile, index=true)
    for record in b.query(contig):
        n_reads += 1
        if record.accept_record():
            n_accept += 1
            var insert_val = abs(record.isize)
            if insert_val <= inserts.len:
                inserts[insert_val-1] += 1
            else:
                overflow.add(insert_val)
    close(b)
    if verbose:
        stderr.write_line fmt"{contig} complete".fgBlue
    return chrom_freqs(inserts: inserts, 
                       overflow: overflow,
                       n_reads: n_reads,
                       n_accept: n_accept)

proc cmd_insert_size*(bamfile: string, distfile: string, verbose: bool, basename: bool, absolute: bool) =
    #[
        Calculates insert size
    ]#
    let fname = bam_file.lastPathPart()
    var 
        b: Bam
        inserts = new_seq[int64](INS_ARR)
        inserts_trimmed = new_seq[int64](INS_ARR)
        overflow: seq[int64]
        n_reads = 0
        n_accept = 0
        max_insert = 0i64
        p99 = 0i64
    
    open(b, bamfile, index=true)
    var freq_results = newSeq[FlowVar[chrom_freqs]](b.hdr.targets().len)
    for idx, contig in b.hdr.targets():
        freq_results[idx] = spawn freq_inserts(bamfile, contig.name, verbose)
    sync()

    # Compile the results
    for idx, contig in b.hdr.targets():
        var result = ^freq_results[idx]
        for i, v in result.inserts:
            inserts[i] += v
        for i, v in result.overflow:
            overflow.add(v)
        n_reads += result.n_reads
        n_accept += result.n_accept
    close(b)

    # Calculate the median insert size
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
            p99 = idx + 1
    

    #echo inserts_trimmed
    let median_insert_size = median_freq(inserts_trimmed)
    let mean_insert_size = mean_freq(inserts_trimmed)
    let min_insert_size = inserts_trimmed.find(inserts_trimmed.filterIt(it > 0)[0]) + 1

    # Output distribution
    # var x: seq[int]
    # var y: seq[int]
    
    if distfile != "":
        var f = open(distfile, fmWrite)
        f.writeLine(output_header(distribution_header, basename, absolute))
        for idx, val in inserts_trimmed.filterIt(it > 0):
            # x.add idx
            # y.add val.int
            f.writeLine(output_w_fnames([$idx, $val, bam_sample(b), fname].join("\t"), bamfile, basename, absolute))
        f.close()

    
    # var df = seqsToDf(x, y).filter( f{ "y" > 0 } )
    # echo df
    # ggplot(df, aes(x = "x1", y = "y1")) +
    # geom_point() +
    # ggsave("scatterColor.pdf")

    # Calc max
    if overflow.len > 0:
        max_insert = max(overflow.mapIt(it.int))
    else:
        max_insert = inserts.len.int64 - inserts.reversed().find(inserts.filterIt(it > 0)[^1])
    
    # Calc sd
    # Standard dev
    let n = sum(inserts[1..p99 - 1])
    var m = 0.int64
    for idx, val in inserts[0..p99 - 1]:
        m += val * (idx + 1)^2
    let variance = (m.float - (n.float*(mean_insert_size^2))) / (n - 1).float
    let std_dev = pow(variance, 0.5)

    var header_out = [$median_insert_size,
                  fmt"{mean_insert_size:0.3f}",
                  fmt"{std_dev:0.3f}",
                  $min_insert_size,
                  $p99,
                  $max_insert,
                  $n_reads,
                  $n_accept,
                  $(sum(inserts_trimmed)),
                  bam_sample(b)]
    echo output_w_fnames(header_out.join("\t"), bamfile, basename, absolute)

