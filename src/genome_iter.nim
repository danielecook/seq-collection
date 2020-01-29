# Outputs ranges of given width
# e.g.
# V:18600000-18699999
# V:18700000-18799999
import strformat
import hts


proc genome_iter*(vcf: VCF, l: int) =
    if l == 0:
        for contig in vcf.contigs:
            echo contig.name
    else:
        var pos = 1
        for contig in vcf.contigs:
            while pos < contig.length.int:
                if pos+l-1 > contig.length.int:
                    echo fmt"{contig.name}:{pos}-{contig.length}"
                else:
                    echo fmt"{contig.name}:{pos}-{pos+l-1}"
                pos += l
            pos = 1


proc genome_iter*(bam: BAM, l: int) =
    if l == 0:
        for contig in bam.hdr.targets:
            echo contig.name
    else:
        var pos: int
        for contig in bam.hdr.targets:
            while pos < contig.length.int:
                if pos+l-1 > contig.length.int:
                    echo fmt"{contig.name}:{pos}-{contig.length-1}"
                else:
                    echo fmt"{contig.name}:{pos}-{pos+l-1}"
                pos += l
            pos = 0