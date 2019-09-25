import memfiles
import streams
import zip/gzipfiles
import utils/helpers
import strformat
import strutils
import os

const header* = ["reads",
                "gc_content",
                "gc_bases",
                "bases",
                "fname"].join("\t")



proc fq_count*(fastq: string) =
    #[
        Deduplicates reads by ID in FASTQs

        Based on Brent Pedersens implementation:
        https://gist.github.com/brentp/640806
    ]#

    var 
        i = 0
        gc_cnt: int64
        n_cnt: int64
        total_len: int64
        n_reads: int
        #n_dups: int

    var basename = lastPathPart(fastq)
    var absolute_path = absolutePath(fastq)

    let stream: Stream =
        if fastq[^3 .. ^1] == ".gz":
            newGZFileStream(fastq)
        else:
            newFileStream(fastq, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fastq, 2)

    for line in lines(stream):
        i.inc()
        if (i mod 4) == 1:
            n_reads.inc()
        if (i mod 4) == 2:
            gc_cnt.inc(line.count("G") + line.count("C"))
            n_cnt.inc(line.count("N"))
            total_len.inc(line.len)
    
    var output = [$n_reads,
                $(gc_cnt.float / (total_len - n_cnt).float),
                $gc_cnt,
                $n_cnt,
                $total_len,
                basename,
                absolute_path]

    echo output.join("\t")
    
    # https://github.com/nim-lang/Nim/issues/9026#issuecomment-423632254
    # var mf = memfiles.open(fn)
    # var cs: cstring
    # var linec, wordc, bytec: int
    # var inWord: bool
    # var s: string
    # for slice in memSlices(mf):
    #   inc(linec)
    #   cs = cast[cstring](slice.data)
    #   let length = slice.size
    #   inc(bytec, length)
    #   var j = -1
    #   for i in 0..length-1:
    #     j = i
    #     if cs[i] in WhiteSpace:
    #       if inWord == true:
    #         inc(wordc)
    #         inWord = false
    #     else:
    #       inWord = true
    #   if j >= 0:
    #     inc(wordc)
    # result.linec = linec
    # result.wordc = wordc
    # result.bytec = bytec + linec
