import sequtils
import strutils
import gz
import os
import re
import streams
import zip/gzipfiles
import utils

const qual = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
const header* = ["machine",
                 "run",
                 "lane",
                 "flowcell",
                 "index1",
                 "index2",
                 "sequencer",
                 "prob_sequencer",
                 "qual_format",
                 "qual_phred",
                 "qual_multiple",
                 "min_qual",
                 "max_qual",
                 "n_lines",
                 "basename",
                 "absolute_path"].join("\t")

import tables

type 
    Fastq* = object
        name*: string
        phred*: string
        minimum*, maximum*: int


var fastq_types* = @[Fastq(name: "Sanger", phred : "Phred+33", minimum: 0, maximum: 40),
                     Fastq(name: "Solexa", phred : "Solexa+64", minimum: 59, maximum: 104),
                     Fastq(name: "Illumina 1.3+", phred: "Phred+64", minimum: 64, maximum: 104),
                     Fastq(name: "Illumina 1.5+", phred: "Phred+64", minimum: 64, maximum: 104),
                     Fastq(name: "Illumina 1.8+", phred: "Phred+33", minimum: 0, maximum: 42)]

type
    Instrument* = object
        pattern*: Regex
        sequencer*: seq[string]


let InstrumentIDs = @[Instrument(pattern: re"HWI-M[0-9]{4}$", sequencer: @["MiSeq"]),
                        Instrument(pattern: re"HWUSI", sequencer: @["Genome Analyzer IIx"]),
                        Instrument(pattern: re"M[0-9]{5}$", sequencer: @["MiSeq"]),
                        Instrument(pattern: re"HWI-C[0-9]{5}$", sequencer: @["HiSeq 1500"]),
                        Instrument(pattern: re"C[0-9]{5}$", sequencer: @["HiSeq 1500"]),
                        Instrument(pattern: re"HWI-D[0-9]{5}$", sequencer: @["HiSeq 2500"]),
                        Instrument(pattern: re"D[0-9]{5}$", sequencer: @["HiSeq 2500"]),
                        Instrument(pattern: re"J[0-9]{5}$", sequencer: @["HiSeq 3000"]),
                        Instrument(pattern: re"K[0-9]{5}$", sequencer: @["HiSeq 3000","HiSeq 4000"]),
                        Instrument(pattern: re"E[0-9]{5}$", sequencer: @["HiSeq X"]),
                        Instrument(pattern: re"NB[0-9]{6}$", sequencer: @["NextSeq"]),
                        Instrument(pattern: re"NS[0-9]{6}$", sequencer: @["NextSeq"]),
                        Instrument(pattern: re"MN[0-9]{5}$", sequencer: @["MiniSeq"])]

type
    Flowcell* = object
        pattern*: Regex
        sequencer*: seq[string]
        description*: string

let FCIDs = @[Flowcell(pattern: re"C[A-Z,0-9]{4}ANXX$", sequencer: @["HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], description: "High Output (8-lane) v4 flow cell"),
              Flowcell(pattern: re"C[A-Z,0-9]{4}ACXX$", sequencer: @["HiSeq 1000", "HiSeq 1500", "HiSeq 2000", "HiSeq 2500"], description: "High Output (8-lane) v3 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}ADXX$", sequencer: @["HiSeq 1500", "HiSeq 2500"], description: "Rapid Run (2-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BCXX$", sequencer: @["HiSeq 1500", "HiSeq 2500"], description: "Rapid Run (2-lane) v2 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BCXY$", sequencer: @["HiSeq 1500", "HiSeq 2500"], description: "Rapid Run (2-lane) v2 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BBXX$", sequencer: @["HiSeq 4000"], description: "(8-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BBXY$", sequencer: @["HiSeq 4000"], description: "(8-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}CCXX$", sequencer: @["HiSeq X"], description: "(8-lane) flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}CCXY$", sequencer: @["HiSeq X"], description: "(8-lane) flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}ALXX$", sequencer: @["HiSeq X"], description: "(8-lane) flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGXX$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGXY$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGX2$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}AFXX$", sequencer: @["NextSeq"], description: "Mid output flow cell"),
              Flowcell(pattern: re"A[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq flow cell"),
              Flowcell(pattern: re"B[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq flow cell"),
              Flowcell(pattern: re"D[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq nano flow cell"),
              Flowcell(pattern: re"G[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq micro flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}DMXX$", sequencer: @["NovaSeq"], description: "S2 flow cell")]


proc qual_to_int(q_score: char): int =
    return qual.find(q_score)

proc qual_min_max*(quality_scores: string, prev_min: int, prev_max: int): (int, int) =
    # Calculate the running min/max
    var qual_scores =  map(toSeq(quality_scores.items), qual_to_int)
    if prev_min >= 0:
        qual_scores = qual_scores.concat(@[prev_min, prev_max])
    return (qual_scores.min(), qual_scores.max())

proc union(a: seq[string], b: seq[string]): seq[string] =
    return a.concat(b).deduplicate()

proc intersect(a: seq[string], b: seq[string]): seq[string] =
    var output: seq[string]
    for i in a:
        for j in b:
            if i == j:
                output.add(i)
    return output.deduplicate()


proc detect_sequencer(machine: string, flowcell: string): (seq[string], string) =
    # Identify sequencer based on machine or flowcell
    # https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
    var seq_by_iid: seq[string]
    var seq_by_fcid: seq[string]
    for k in InstrumentIDs:
        if re.match(machine, k.pattern):
            for s in k.sequencer:
                seq_by_iid.add(s)

    for k in FCIDs:
        if re.match(flowcell, k.pattern):
            for s in k.sequencer:
                seq_by_fcid.add(s)

    # Both empty
    if concat(seq_by_iid, seq_by_fcid).len == 0:
        return (@[], "")
    
    if seq_by_iid.len == 0:
        return (seq_by_fcid, "likely:flowcell")

    if seq_by_fcid.len == 0:
        return (seq_by_iid, "likely:machine")

    var sequencers = intersect(seq_by_iid, seq_by_fcid)
    if sequencers.len > 0:
        return (sequencers, "high:machine+flowcell")
    else:
        return (union(seq_by_iid, seq_by_fcid), "uncertain")

    return (@[], "")



proc fq_meta*(fastq: string, sample_n = 20) =

    var basename = lastPathPart(fastq)
    var absolute_path = absolutePath(fastq)

    let stream: Stream =
        if fastq[^3 .. ^1] == ".gz":
            newGZFileStream(fastq)
        else:
            newFileStream(fastq, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fastq, 2)
    
    var
        machine: string
        run: string
        lane: string
        flowcell: string
        qual_min = -1
        qual_max = -1
        line: string
        sequencer_list: seq[string]
        sequencer: string
        sequencer_prob: string
        barcodes = newSeq[string](sample_n)
        i = 0
    
    while not stream.atEnd() and i < sample_n * 4:
        line = stream.readLine()
        if i %% 4 == 0:
            try:
                if i == 0:
                    let qual_line = line.split(":")
                    machine = qual_line[0].strip(chars = {'@'})
                    run = qual_line[1]
                    flowcell = qual_line[2]
                    if '_' in flowcell:
                        flowcell = flowcell.split("_")[^1]
                    lane = qual_line[3]
                barcodes[i.div(4)] = line.split(":")[^1]
            except IndexError:
                discard
            
        # Quality scores
        if i %% 4 == 3:
            (qual_min, qual_max) = qual_min_max(line, qual_min, qual_max)

        i.inc()
    stream.close()

    if machine != "" and flowcell != "":
        (sequencer_list, sequencer_prob) = detect_sequencer(machine, flowcell)
        if sequencer_list.len >= 1:
            sequencer = sequencer_list[^1]

    var fastq_scores = fastq_types.filterIt(qual_min >= it.minimum and qual_max <= it.maximum)
    var most_comm_barcode = barcodes.newCountTable().largest()[0]
    var fastq_scores_name = fastq_scores.mapIt(it.name).join(";")
    var fastq_scores_phred = fastq_scores.mapIt(it.phred).join(";")

    echo [machine,
          $run,
          $lane,
          flowcell,
          most_comm_barcode,
          "",
          sequencer,
          sequencer_prob,
          fastq_scores_name,
          fastq_scores_phred,
          $(fastq_scores.mapIt(it.name).len > 1),
          (if qual_min >= 0: $qual_min else: ""),
          (if qual_max >= 0: $qual_max else: ""),
          $(i/4).int,
          basename,
          absolute_path].join("\t")