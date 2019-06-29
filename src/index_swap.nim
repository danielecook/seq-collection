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

type alt_error* = object
    alt*: RunningStat
    dp*: RunningStat # depth of genotyped sites
    af*: RunningStat

type alt_contamination* = object
    alt*: RunningStat # apparent alt calls in ref case
    dp*: RunningStat # depth of genotyped sites in ref case value
    af*: RunningStat # AF; alt/dp
    alt_load*: RunningStat # Alt bases found in other samples that are apparent het/hom alt.

proc index_swaps*(bams: seq[string], sites_path: string, fasta: string, debug: bool, threads: int) =
    
    var
        threads = 2
        n_samples = bams.len
        group: string
        sample_names = newSeqOfCap[string](n_samples)
        sample_flowcells = newSeqOfCap[string](n_samples)
        results = newSeq[seq[int8]](n_samples)
        alt_depth = newSeq[seq[int]](n_samples)
        alt_alleles = newSeq[seq[int]](n_samples)
        error_stat = newSeq[alt_error](n_samples)
        error_stat_aggregate: alt_error
        contamination_stat = newSeq[alt_contamination](n_samples)
        contamination_aggregate: alt_contamination
        depth = newSeq[seq[int]](n_samples)
        exit = false
        stats = newSeq[Stat4](n_samples)
        responses = newSeq[FlowVarBase](n_samples)
    
    if bams.len == 0:
        print_error(fmt"No Bams Specified")
        quit()

    for i in bams:
        # Check that file exists
        if not os.existsFile(i):
            print_error(fmt": File not found: {i}")
            exit = true
    if exit:
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
  
    if debug and sitelist.len > 2000:
        stderr.write_line "Debug".bgWhite.fgGreen & " Truncating site list".fgGreen
        sitelist = sitelist[0..<2000]

    for j in 0..<bams.len:
        # Get sample name
        sample_names.add(get_sample_names(bams[j]))

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
                                          10)
        
    echo sample_names, "<-- SAMPLE NAMES"
    
    for index, fv in responses:
        blockUntil(fv)
    
    var error_alt_bases: RunningStat
    var error_af_freq: RunningStat
    var error_dp_freq: RunningStat
    var contamination_freq: RunningStat

    var last_pos = 101
    var header_line = @["sample",
                        "site",
                        "af",
                        "gt",
                        "depth",
                        "alt_load",
                        "group"]
    echo header_line.join("\t")
    for rowi in 0..sitelist.high:

        # Pull out 'reference-like' sites
        # Infer that ALT
        if all(lc[(results[i][rowi] in [0,3]) | (i <- 0..<n_samples), bool], proc(x: bool): bool = return x):
            for i in 0..<n_samples:
                error_stat[i].dp.push(depth[i][rowi])
                error_stat[i].alt.push(alt_alleles[i][rowi])
                error_stat_aggregate.dp.push(depth[i][rowi])
                error_stat_aggregate.alt.push(alt_alleles[i][rowi])

        for i in 0..<n_samples:
            if results[i][rowi] in [0, 3]:

                # Contamination
                contamination_stat[i].dp.push(depth[i][rowi])
                contamination_stat[i].alt.push(alt_alleles[i][rowi])
                var af = float(alt_alleles[i][rowi]) / float(depth[i][rowi])
                contamination_stat[i].af.push(af)

                var alt_load = 0
                for j in 0..<n_samples:
                    if results[j][rowi] in [1,2]:
                        alt_load += alt_depth[j][rowi]

                if alt_load > 0 and af > 0:
                    group = "switch"
                elif alt_load == 0 and af > 0:
                    group = "technical"
                else:
                    group = "NA"
                
                contamination_stat[i].alt_load.push(alt_load)
                var out_line = @[sample_names[i],
                                 fmt"{sitelist[i].chrom}:{sitelist[i].position}",
                                 $af,
                                 $alt_alleles[i][rowi],
                                 $depth[i][rowi],
                                 $alt_load,
                                 $group]
                echo out_line.join("\t")
    
    stderr.write_line fmt"Analysis complete {sites_path}".fgGreen
  
    shallow(results)
