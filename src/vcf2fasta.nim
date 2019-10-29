import hts
import json
import sequtils
import strutils
import strformat
import streams
import utils/helpers
import constants


proc to_fasta*(vcf: string, region_list: seq[string], sample_set: string, force: bool) =
    # vcf - VCF input
    # region_list
    # sample_set
    # force - Force even if is not phased
    var v:VCF
    doAssert open(v, vcf)
    if sample_set != "ALL":
        let samples_keep = filterIt(sample_set.split({',', ' '}), it.len > 0)
        set_samples(v, samples_keep)

    var gts = new_seq[int32](10)

    var allele_set: seq[string]
    var sample: string
    # sample → chromosome_n → genotype

    var n_samples = v.samples.len
    # Initialize streams
    echo n_samples
    var sequences = new_seq[seq[FileStream]](n_samples)
    var chrom_set: seq[FileStream]
    for n_sample in 0..<n_samples:
        sample = v.samples[n_sample]
        chrom_set = @[]
        for n_chrom in 0..<2:
            echo fmt"{n_sample} {n_chrom} - {v.samples[n_sample]}"
            chrom_set.add(newFileStream(fmt"{sample}_{n_chrom}.fa", fmWrite))
        sequences[n_sample] = chrom_set

    var n_chrom = 0
    var n_sample = 0
    var allele_out: string
    for rec in variants(v, region_list):
        allele_set = sequtils.concat(@[rec.REF], rec.ALT)
        # Record ploidy
        n_sample = 0
        for g in rec.format.genotypes(gts):
            n_chrom = 0
            for a in g:
                if a.phased == false and force == false:
                    quit_error("Genotypes are not phased", 99)
                else:
                    if a.value() != -1:
                        allele_out = allele_set[a.value()]
                    else:
                        allele_out = "N"
                    echo fmt"{allele_out} → {n_sample} → {n_chrom}"
                    sequences[n_sample][n_chrom].write(allele_out)
                n_chrom.inc()
            n_sample.inc()
