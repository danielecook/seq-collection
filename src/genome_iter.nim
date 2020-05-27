# Outputs ranges of given width
# e.g.
# V:18600000-18699999
# V:18700000-18799999
import strformat
import hts


proc genome_iter*(f: Fai, l: int) =
    if l == 0:
        for i in 0..<f.len:
            echo f[i]
    else:
        var pos = 1
        for i in 0..<f.len:
            var contig_len = f.chrom_len(f[i])
            var contig_name = f[i]
            while pos < contig_len:
                if pos+l-1 > contig_len:
                    echo fmt"{contig_name}:{pos}-{contig_len}"
                else:
                    echo fmt"{contig_name}:{pos}-{pos+l-1}"
                pos += l
            pos = 1


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