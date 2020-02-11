import hts
import tables
import sequtils
import strutils
import strformat
import algorithm
import sequtils
import math

# Largely ported from https://github.com/edgardomortiz/vcf2phylip

# Dictionary of IUPAC ambiguities for nucleotides
# '*' means deletion for GATK (and other software?)
# Deletions are ignored when making the consensus
# Dictionary to translate IUPAC ambiguities, lowercase letters are used when "*" or "N" were present for a position,
# however, software like Genious for example are case insensitive and will imply ignore capitalization
const AMB = {".": '.', "*":'-', "A":'A', "C":'C', "G":'G',"N":'N',"T":'T',
             "*A":'a',"*C":'c',"*G":'g',"*N":'n',"*T":'t',"AC":'M',
             "AG":'R',"AN":'a',"AT":'W',"CG":'S',"CN":'c',"CT":'Y',
             "GN":'g',"GT":'K',"NT":'t',"*AC":'m',"*AG":'r',"*AN":'a',
             "*AT":'w',"*CG":'s',"*CN":'c',"*CT":'y',"*GN":'g',"*GT":'k',
             "*NT":'t',"ACG":'V',"ACN":'m',"ACT":'H',"AGN":'r',"AGT":'D',
             "ANT":'w',"CGN":'s',"CGT":'B',"CNT":'y',"GNT":'k',"*ACG":'v',
             "*ACN":'m',"*ACT":'h',"*AGN":'r',"*AGT":'d',"*ANT":'w',"*CGN":'s',
             "*CGT":'b',"*CNT":'y',"*GNT":'k',"ACGN":'v',"ACGT":'N',"ACNT":'h',
             "AGNT":'d',"CGNT":'b',"*ACGN":'v',"*ACGT":'N',"*ACNT":'h',"*AGNT":'d',
             "*CGNT":'b',"*ACGNT":'N'}.toTable

proc is_snp(rec: Variant): bool = 
    for i in @[rec.REF].concat(rec.ALT):
        if i.len != 1:
            return false
    return true

proc missing_to_n(x: char): char = 
    if x == '.':
        return 'N' 
    else:
        return x

proc vcf2phylo*(vcf: string, region: seq[string]) =
    var v:VCF
    doAssert open(v, vcf)

    var snp_num = 0
    var snp_accepted = 0
    var snp_shallow = 0
    var snp_multinuc = 0
    var snp_biallelic = 0

    var n_rec = 0
    var gts = new_seq[int32](5)
    var gt_out: char
    var gt_set: seq[char]
    var transpose = new_seq[new_seq[char](0)](v.samples.len)

    for rec in v:
        discard rec.format.get("GT", gts)

        # Check if this is a SNP
        if rec.is_snp == false:
            continue

        # Check minimum number of samples
        var alleles = @[rec.REF].concat(rec.ALT)

        n_rec += 1
        if n_rec mod 1000 == 0:
            stderr.writeLine(fmt"Processed {n_rec} variants")

        gt_set.set_len(0)
        for idx, g in rec.format.genotypes(gts).toSeq():
            var gt: seq[char]
            for a in g:
                var agt = if a.value() >= 0: alleles[a.value()] else: "."
                if agt.len > 0:
                    gt.add(agt)
            try:
                gt.sort()
                gt_set.add(AMB[$gt.deduplicate().join("")])
            except KeyError:
                # Variant malformatted somehow
                continue

        #var n_missing = gt_set.count('.')
        #if n_missing <= 4:
        #    continue
        
        # Once checked for missing, transpose genotypes
        gt_set.applyIt(if it == '.': 'N' else: it)
        for idx, gt_in in gt_set:
            # Convert missing to 'N'
            transpose[idx].add(gt_in)

    var len_longest_sample_name = max(v.samples.mapIt(it.len))
    
    # Transpose and output
    echo &"{v.samples.len} {transpose[0].len}"
    for idx, gt_set in transpose:
        var nt_seq = gt_set.join("")
        echo &"""{alignLeft(v.samples[idx], len_longest_sample_name+3)}{nt_seq}"""


    

    # Write phylip output

            #echo gt.split("|")
            #tgt_set.add(gt)
        #echo tgt_set

        #var alleles = @[rec.REF].concat(rec.ALT)
        #site_tmp = ''.join([amb[''.join(sorted(set([nuc[broken[i][j]] for j in gt_idx])))] for i in range(9, index_last_sample)])
