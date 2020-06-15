import tables
import hts
import sequtils
import strutils
import strformat
import math
import stats
import os
import colorize
import regex
import algorithm
import threadpool
import utils/helpers

const INS_ARR = 10000

# const insert_size_header* = ["median",
#                  "mean",
#                  "std_dev",
#                  "min",
#                  "percentile_99.5",
#                  "max_all",
#                  "n_reads",
#                  "n_accept",
#                  "n_use",
#                  "sample"].join("\t")

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

proc covariance(a: seq[uint], b: seq[uint]): float =
    var covar: float
    var len_arr = a.len.float
    var mean_a = sum(a).float / len_arr
    var mean_b = sum(b).float / len_arr
    for idx in 0..(len_arr.int - 1):
        covar += ((a[idx].float - mean_a) * (b[idx].float - mean_b)).float / len_arr
    return covar

proc stddev(a: seq[uint]): float = 
    var
        len_arr = a.len.float
        mean_a = sum(a).float / len_arr
        s: float
    for i in a:
        s += (i.float - mean_a)^2
    return pow(s / len_arr, 0.5)

proc pcorr(a: seq[uint], b: seq[uint]): float =
    # Pearson correlation coefficient
    return covariance(a, b) / (stddev(a) * stddev(b))

iterator pos(length: uint, amt = 10000.uint): uint =
    var i = 0.uint
    while i <= length:
        i += amt
        yield i

proc read_groups(b: Bam): seq[string] =
    var rg: seq[string]
    for line in splitLines($(b.hdr)):
        if line.startsWith("@RG"):
            for field in line.split("\t"):
                if field.startswith("ID:"):
                    rg.add(field.replace("ID:", ""))
    return rg

proc depth(b: Bam, contig: string, pos: int, rgs: seq[string]): seq[uint] =
    var d = new_seq[uint](rgs.len)
    for read in b.query(contig, pos, pos+1):
        if read.start == pos:
            let rg = tag[string](read, "RG").get()
            d[rgs.find(rg)] += 1.uint
    return d

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

proc library_id*(bamfile: string, verbose: bool, basename: bool, absolute: bool) =
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
    #var freq_results = newSeq[FlowVar[chrom_freqs]](b.hdr.targets().len)
    var rgs = read_groups(b)
    var d: seq[uint]
    var depth_profile = new_seq[seq[uint]](rgs.len)
    for idx, contig in b.hdr.targets():
        for i in pos(contig.length):
            d = depth(b, contig.name, i.int, rgs)
            if sum(d) > 0.uint:
                for idx, val in d:
                    depth_profile[idx].add(val)
    for i, x in depth_profile:
        for j, y in depth_profile:
            if i != j:
                echo pcorr(x, y), " : ", rgs[i], " → ", rgs[j]
    #echo covariance(depth_profile[0], depth_profile[1])
        #freq_results[idx] = spawn freq_inserts(bamfile, contig.name, verbose)
    #sync()

    # Compile the results
    # for idx, contig in b.hdr.targets():
    #     var result = ^freq_results[idx]
    #     for i, v in result.inserts:
    #         inserts[i] += v
    #     for i, v in result.overflow:
    #         overflow.add(v)
    #     n_reads += result.n_reads
    #     n_accept += result.n_accept
    # close(b)

    # var header_out = [$median_insert_size,
    #               $math.round(mean_insert_size, 3),
    #               $round(std_dev, 3),
    #               $min_insert_size,
    #               $p99,
    #               $max_insert,
    #               $n_reads,
    #               $n_accept,
    #               $(sum(inserts_trimmed)),
    #               bam_sample(b)]
    # output_w_fnames(header_out.join("\t"), bamfile, basename, absolute)

