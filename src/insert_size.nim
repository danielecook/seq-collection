import tables
import hts
import sequtils
import math
import stats

proc accept_record(r: Record): bool = 
    if (r.flag.pair == false or
        r.flag.unmapped == true or
        r.flag.mate_unmapped or
        r.flag.read1 == true or
        r.flag.secondary or
        r.flag.supplementary or
        r.flag.dup or
        r.isize == 0):
        return false
    return true

proc cmd_insert_size*(bamfile: string) =
    # [ ] TODO: Reimplement this with a count table?
    #  create procs for MAD (median abs. dev)
    #  and use to calculate width...?
    #  proc for getting median also?
    #  Proc for freq proc.
    #  https://nimble.directory/pkg/rbtree
    # https://nimble.directory/pkg/binaryheap

    # https://leetcode.com/problems/find-median-from-data-stream/solution/
    var b: Bam
    let ins_arr = 100000
    const deviations = 10
    var inserts = new_seq[int64](ins_arr)
    var insert_lengths = new_seq[int64](ins_arr)
    var percentile = new_seq[float](ins_arr)
    var max_insert_size = 0
    open(b, bamfile, index=true)
    for record in b:
       if record.accept_record():
            var idx = abs(record.isize)
            if idx <= inserts.len:
                inserts[idx-1] += 1
            else:
                if idx > max_insert_size:
                    max_insert_size = idx
    
    var total_length = 0.int64
    var trimmed_count = 0.int64
    var trimmed_length = 0.int64
    for idx, val in inserts:
        total_length += val.int64*(idx+1)

    # Trim 1 pct from end.
    var running_total = 0i64
    var max_idx = 0
    for idx, val in inserts:
        var length_add = val.int64*(idx+1)
        running_total += length_add
        var pct = running_total.float / total_length.float
        if pct <= 0.995:
            percentile[idx] = pct
            trimmed_count += val.int64
            trimmed_length += length_add
            insert_lengths[idx] = length_add
            max_idx = idx
        else:
            percentile[idx] = 1.0

    var median_insert_size: int64
    for idx, val in percentile:
        if val <= 0.50:
            median_insert_size = idx + 1
        else:
            break
    echo median_insert_size, "<- median"
    # Stats
    let mode_insert_size = inserts.find(max(inserts))+1       
    let min_insert_size = inserts.find(inserts.filterIt( it > 0 )[0]) + 1
    let mean_insert_size = (trimmed_length.float / trimmed_count.float)

    # Standard dev
    let n = sum(inserts[1..max_idx])
    var m = 0.int64
    for idx, val in inserts[1..max_idx]:
        m += val * (idx + 1)^2
    let variance = (m.float - (n.float*(mean_insert_size^2))) / (n - 1).float
    let std_dev = pow(variance, 0.5)

    var result = [
        $min_insert_size,
        $mean_insert_size,
        $median_insert_size,
        $mode_insert_size,
        $std_dev,
        $max_insert_size,
        bamfile
    ]
    
    echo $result


