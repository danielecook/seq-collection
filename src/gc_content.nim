# Author: Daniel E. Cook
# Calculates GC content for a given window or set of windows
import strformat
import hts



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