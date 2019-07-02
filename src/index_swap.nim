import os
import hts
import math
import sequtils
import sugar
import utils/sites
import utils/helpers
import colorize
import argparse
import strutils
import parseutils
import algorithm
import strformat
import threadpool
import stats
import progress
import tables

type alt_error* = object
    alt*: RunningStat
    dp*: RunningStat # depth of genotyped sites
    af*: RunningStat

type alt_contamination* = object
    alt*: RunningStat # apparent alt calls in ref case
    dp*: RunningStat # depth of genotyped sites in ref case value
    af*: RunningStat # AF; alt/dp
    pool_alt_depth*: RunningStat # Alt bases found in other samples that are apparent het/hom alt.

proc index_swaps*(bams: seq[string], sites_path: string, fasta: string, threads: int) =
    
    var
        threads = 2
        n_samples = bams.len
        group: string
        sample_names = newSeqOfCap[string](n_samples)
        sample_flowcells = newSeqofCap[string](n_samples)
        results = newSeq[seq[int8]](n_samples)
        alt_depth = newSeq[seq[int]](n_samples)
        alt_alleles = newSeq[seq[int]](n_samples)
        error_stat = newSeq[alt_error](n_samples)
        error_stat_aggregate: alt_error
        contamination_stat = newSeq[alt_contamination](n_samples)
        depth = newSeq[seq[int]](n_samples)
        stats = newSeq[Stat4](n_samples)
        responses = newSeq[FlowVarBase](n_samples)
    
    if check_file_list(bams):
        quit()
    
    if threads < 2:
        threads = 1
  
    if threads > n_samples:
        threads = n_samples
  
    # Read sites
    var fai: Fai
    if fasta.len > 0:
        discard open(fai, fasta)
    var sitelist = readSites(sites_path, fai)

    # Get sample names and flowcells
    for bam in bams:
        sample_names.add(get_sample_names(bam))
        sample_flowcells.add(get_sample_flowcells(bam))

    # Create a progress bar
    var progress_bar = newProgressBar(total = sitelist.len * sample_names.len, output=stderr)
    progress_bar.start()

    for j in 0..<bams.len:
        # Work on results
        results[j] = newSeq[int8](sitelist.len)
        alt_depth[j] = newSeq[int](sitelist.len)
        alt_alleles[j] = newSeq[int](sitelist.len)
        depth[j] = newSeq[int](sitelist.len)
        responses[j] = spawn get_bam_alts(bams[j],
                                          fasta,
                                          sitelist,
                                          results[j].addr,
                                          alt_depth[j].addr,
                                          alt_alleles[j].addr,
                                          depth[j].addr,
                                          stats[j].addr,
                                          10,
                                          progress_bar.addr)
        
    for index, fv in responses:
        blockUntil(fv)
    
    var header_line = @["sample",
                        "site",
                        "ref",
                        "alt",
                        "sample_af",
                        "sample_alt_reads",
                        "sample_depth",
                        "pool_af",
                        "pool_alt_depth",
                        "pool_total_depth",
                        "pool_ref",
                        "pool_ref_contaminated",
                        "pool_het",
                        "pool_alt",
                        "pool_fail",
                        "group",
                        "flowcell"].join("\t")
    echo header_line
    for site_n in 0..sitelist.high:

        # Pull out 'reference-like' sites
        # Infer that ALT

        ## -1 = Below min-depth threshhold
        ## 0 = HOM REF
        ## 1 = HET
        ## 2 = HOM ALT
        ## 3 = CONTAMINATED REF 

        ## Population 

        # if all(lc[(results[i][site_n] in [0,3]) | (i <- 0..<n_samples), bool], proc(x: bool): bool = return x):
        #     for i in 0..<n_samples:
        #         error_stat[i].dp.push(depth[i][site_n])
        #         error_stat[i].alt.push(alt_alleles[i][site_n])
        #         error_stat_aggregate.dp.push(depth[i][site_n])
        #         error_stat_aggregate.alt.push(alt_alleles[i][site_n])
        var x = 1
        for i in 0..<n_samples:
            if results[i][site_n] in [0, 3]:

                # Contamination
                contamination_stat[i].dp.push(depth[i][site_n])
                contamination_stat[i].alt.push(alt_alleles[i][site_n])
                var af = float(alt_alleles[i][site_n]) / float(depth[i][site_n])
                contamination_stat[i].af.push(af)

                var pool_alt_depth = 0
                var pool_total_depth = 0
                var flowcell_gts: seq[int]
                for j in 0..<n_samples:
                    # Don't include the sample of interest
                    # Only compare samples sequenced on the same flowcell
                    if i != j and sample_flowcells[i] == sample_flowcells[j]:
                        var gt = results[j][site_n]
                        flowcell_gts.add gt
                        if gt in [1,2]:
                            pool_alt_depth += alt_depth[j][site_n]
                        pool_total_depth += depth[j][site_n]
                
                var gt_counts = newCountTable(flowcell_gts)

                if pool_alt_depth > 0 and af > 0:
                    group = "switch"
                elif pool_alt_depth == 0 and af > 0:
                    group = "technical"
                else:
                    group = "NA"
                
                contamination_stat[i].pool_alt_depth.push(pool_alt_depth)
                var out_line = @[sample_names[i],
                                 fmt"{sitelist[site_n].chrom}:{sitelist[site_n].position}",
                                 $sitelist[site_n].ref_allele,
                                 $sitelist[site_n].alt_allele,
                                 $af,
                                 $alt_alleles[i][site_n],
                                 $depth[i][site_n],
                                 $(pool_alt_depth.float / pool_total_depth.float),
                                 $pool_alt_depth,
                                 $pool_total_depth,
                                 $(gt_counts[0]),
                                 $(gt_counts[3]),
                                 $(gt_counts[1]),
                                 $(gt_counts[2]),
                                 $(gt_counts[-1]),
                                 $group,
                                 sample_flowcells[i]]
                echo out_line.join("\t")
    
    stderr.write_line fmt"Analysis complete {sites_path}".fgGreen
  
    shallow(results)
