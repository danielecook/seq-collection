import sequtils
import strutils
import utils/gz
import os
import re
import sets
import streams
import zip/gzipfiles
import utils/helpers

const qual = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
const fq_meta_header* = ["machine",
                         "sequencer",
                         "prob_sequencer",
                         "flowcell",
                         "flowcell_description",
                         "run",
                         "lane",
                         "sequence_id",
                         "index1",
                         "index2",
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
                      Instrument(pattern: re"HWUSI", sequencer: @["GenomeAnalyzerIIx"]),
                      Instrument(pattern: re"M[0-9]{5}$", sequencer: @["MiSeq"]),
                      Instrument(pattern: re"A[0-9]{5}$", sequencer: @["NovaSeq"]),
                      Instrument(pattern: re"HWI-C[0-9]{5}$", sequencer: @["HiSeq1500"]),
                      Instrument(pattern: re"C[0-9]{5}$", sequencer: @["HiSeq1500"]),
                      Instrument(pattern: re"HWI-D[0-9]{5}$", sequencer: @["HiSeq2500"]),
                      Instrument(pattern: re"D[0-9]{5}$", sequencer: @["HiSeq2500"]),
                      Instrument(pattern: re"J[0-9]{5}$", sequencer: @["HiSeq3000"]),
                      Instrument(pattern: re"K[0-9]{5}$", sequencer: @["HiSeq3000","HiSeq4000"]),
                      Instrument(pattern: re"E[0-9]{5}$", sequencer: @["HiSeqX"]),
                      Instrument(pattern: re"NB[0-9]{6}$", sequencer: @["NextSeq"]),
                      Instrument(pattern: re"NS[0-9]{6}$", sequencer: @["NextSeq"]),
                      Instrument(pattern: re"MN[0-9]{5}$", sequencer: @["MiniSeq"])]

type
    Flowcell* = object
        pattern*: Regex
        sequencer*: seq[string]
        description*: string

# source 0: https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py (everything not labeled)
# source 1: https://github.com/CFSAN-Biostatistics/snp-pipeline/blob/f9bd23caaf0f84deea1a6593dadb5d07b00e75e4/snppipeline/fastq.py#L60
let FCIDs = @[Flowcell(pattern: re"AAXX$", sequencer: @["GenomeAnalyzer"], description: ""), # <-- Source 1
              Flowcell(pattern: re"C[A-Z,0-9]{4}ANXX$", sequencer: @["HiSeq1500", "HiSeq2000", "HiSeq2500"], description: "High Output (8-lane) v4 flow cell"),
              Flowcell(pattern: re"C[A-Z,0-9]{4}ACXX$", sequencer: @["HiSeq1000", "HiSeq1500", "HiSeq2000", "HiSeq2500"], description: "High Output (8-lane) v3 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}ADXX$", sequencer: @["HiSeq1500", "HiSeq2500"], description: "Rapid Run (2-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BCXX$", sequencer: @["HiSeq1500", "HiSeq2500"], description: "Rapid Run (2-lane) v2 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BCXY$", sequencer: @["HiSeq1500", "HiSeq2500"], description: "Rapid Run (2-lane) v2 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BBXX$", sequencer: @["HiSeq4000"], description: "(8-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BBXY$", sequencer: @["HiSeq4000"], description: "(8-lane) v1 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}CCXX$", sequencer: @["HiSeqX"], description: "(8-lane) flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}CCXY$", sequencer: @["HiSeqX"], description: "(8-lane) flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}ALXX$", sequencer: @["HiSeqX"], description: "(8-lane) flow cell"),
              # source 1
              Flowcell(pattern: re"H[A-Z,0-9]{4}AGXX$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGXX$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGXY$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}BGX2$", sequencer: @["NextSeq"], description: "High output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}AFXX$", sequencer: @["NextSeq"], description: "Mid output flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}DMXX$", sequencer: @["NovaSeq"], description: "S2 flow cell"),
              Flowcell(pattern: re"H[A-Z,0-9]{4}DSXX$", sequencer: @["NovaSeq"], description: "S2 flow cell"),
              Flowcell(pattern: re"^A[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq flow cell"),
              Flowcell(pattern: re"^B[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq flow cell"),
              Flowcell(pattern: re"^D[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq nano flow cell"),
              Flowcell(pattern: re"^G[A-Z,0-9]{4}$", sequencer: @["MiSeq"], description: "MiSeq micro flow cell")]

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


proc detect_sequencer(machine: string, flowcell: string): (seq[string], string, string) =
    # Identify sequencer based on machine or flowcell
    # https://github.com/10XGenomics/supernova/blob/master/tenkit/lib/python/tenkit/illumina_instrument.py
    var seq_by_iid: seq[string]
    var seq_by_fcid: seq[string]
    var flowcell_description: string
    for k in InstrumentIDs:
        if re.find(machine, k.pattern) > -1:
            for s in k.sequencer:
                seq_by_iid.add(s)

    for k in FCIDs:
        if re.find(flowcell, k.pattern) > -1:
            flowcell_description = k.description
            for s in k.sequencer:
                seq_by_fcid.add(s)
    
    # Both empty
    if concat(seq_by_iid, seq_by_fcid).len == 0:
        return (@[], "", "")

    if seq_by_iid.len == 0:
        return (seq_by_fcid, "likely:flowcell", flowcell_description)

    if seq_by_fcid.len == 0:
        return (seq_by_iid, "likely:machine", flowcell_description)

    var sequencers = intersect(seq_by_iid, seq_by_fcid)
    if sequencers.len > 0:
        return (sequencers, "high:machine+flowcell", flowcell_description)
    else:
        return (union(seq_by_iid, seq_by_fcid), "uncertain", "")


proc extract_read_info*(line: string): (string, string, string, string, string) =
    # Parses FASTQ read lines and send back output
    var
        qual_line = line.split({':', '/', '#'})
        sequence_id: string
        machine: string
        run: string
        lane: string
        flowcell: string
    
    if qual_line.len == 1:
        sequence_id = qual_line[0].strip(chars = {'@'})
    elif qual_line.len > 1:
        machine = qual_line[0].strip(chars = {'@'})
        if '/' in line:
            # @HWUSI-EAS100R:6:73:941:1973#ATGGGC/1
            # machine:lane:tile:x:y#index/read
            lane = qual_line[1]
        else:
            # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG 
            # @EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:1
            # @D00446:1:140101_HWI-D00446_0001_C8HN4ANXX:8:2210:1238:2018 1:Y:0:GCTCGGTA
            run = qual_line[1]
            flowcell = qual_line[2]
            if '_' in flowcell:
                flowcell = flowcell.split("_")[^1]
            lane = qual_line[3]
    return (sequence_id, machine, run, lane, flowcell)

proc get_sequencer_name(sequencers: seq[string]): string =
    const set1 = ["HiSeq2000", "HiSeq2500"]
    const set2 = ["HiSeq1500", "HiSeq2500"]
    const set3 = ["HiSeq3000", "HiSeq4000"]
    #const set4 = ["HiSeq1000", "HiSeq1500"]
    #const set5 = ["HiSeq1000", "HiSeq2000"]
    
    if (set1.anyIt(it in sequencers) == true):
        return "HiSeq2000/2500"
    elif set2.anyIt(it in sequencers): 
        return "HiSeq1500/2500"
    elif set3.anyIt(it in sequencers): 
        return "HiSeq3000/4000"
    elif sequencers.len > 0:
        return sequencers[^1]

proc fq_meta*(fastq_in: string, sample_n = 20, follow_symlinks: bool) =

    var
        fastq: string
        sequence_id: string
        machine: string
        run: string
        lane: string
        flowcell: string
        barcode: string
        barcode_re = re"[ATCGN\+\-]{3,12}+"
        qual_min = -1
        qual_max = -1
        line: string
        most_comm_barcode: string
        sequencer_list: seq[string]
        sequencer: string
        sequencer_prob: string
        flowcell_description: string
        barcodes = newSeq[string](sample_n)
        i = 0

    if follow_symlinks and symlinkExists(fastq_in):
        fastq = expandSymlink(fastq_in)
    else:
        fastq = fastq_in

    var basename = lastPathPart(fastq)
    var absolute_path = absolutePath(fastq)


    let stream: Stream =
        if fastq[^3 .. ^1] == ".gz":
            newGZFileStream(fastq)
        else:
            newFileStream(fastq, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fastq, 2)

    
    while not stream.atEnd() and i < sample_n * 4:
        line = stream.readLine()
        if i %% 4 == 0:
            try:
                if i == 0:
                    (sequence_id, machine, run, lane, flowcell) = extract_read_info(line)
                var qual_line = line.split({':', '/', '#'})
                if qual_line.len > 2:
                    if '/' in line:
                        barcode = qual_line[^2]
                    else:
                        barcode = qual_line[^1]
                    if re.match(barcode, barcode_re):
                        barcodes[i.div(4)] = barcode
            except IndexError:
                raise
        
        # Quality scores
        if i %% 4 == 3:
            (qual_min, qual_max) = qual_min_max(line, qual_min, qual_max)

        i.inc()
    stream.close()

    if machine != "" or flowcell != "":
        (sequencer_list, sequencer_prob, flowcell_description) = detect_sequencer(machine, flowcell)
        sequencer = get_sequencer_name(sequencer_list)

    var fastq_scores = fastq_types.filterIt(qual_min >= it.minimum and qual_max <= it.maximum)
    var barcode_set = barcodes.filterIt(it != "")
    if barcode_set.len > 0:
        most_comm_barcode = barcode_set.newCountTable().largest()[0]
    var fastq_scores_name = fastq_scores.mapIt(it.name).join(";")
    let fastq_scores_phred = fastq_scores.mapIt(it.phred).deduplicate().join(";")

    echo [machine,
          sequencer,
          sequencer_prob,
          flowcell,
          flowcell_description,
          $run,
          $lane,
          sequence_id,
          most_comm_barcode,
          "",
          fastq_scores_name,
          fastq_scores_phred,
          $(fastq_scores.mapIt(it.name).len > 1),
          (if qual_min >= 0: $qual_min else: ""),
          (if qual_max >= 0: $qual_max else: ""),
          $(i/4).int,
          basename,
          absolute_path].join("\t")
