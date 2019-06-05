import sequtils

const qual = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""

type 
    Fastq* = object
        name*: string
        phred*: string
        minimum*, maximum*: int

#var fastq_types: seq[Fastq]
var fastq_types* = @[Fastq(name : "Sanger", phred : "Phred+33", minimum : 0, maximum : 40),
                    Fastq(name : "Solexa", phred : "Solexa+64", minimum : 59, maximum : 104),
                    Fastq(name : "Illumina 1.3+", phred: "Phred+64", minimum : 64, maximum : 104),
                    Fastq(name : "Illumina 1.5+", phred: "Phred+64", minimum : 64, maximum : 104),
                    Fastq(name : "Illumina 1.8+", phred: "Phred+33", minimum : 0, maximum : 41)]

proc qual_to_int(q_score: char): int =
    return qual.find(q_score)

proc qual_min_max*(quality_scores: string, prev_min: int, prev_max: int): (int, int) =
    # Calculate the running min/max
    var a = toSeq(quality_scores.items)
    echo fastq_types
    echo prev_min, " ", prev_max
    var qual_scores =  map(toSeq(quality_scores.items), qual_to_int).concat(@[prev_min, prev_max])
    return (qual_scores.min(), qual_scores.max())
