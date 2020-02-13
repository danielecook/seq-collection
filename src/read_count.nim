import hts

# How to org:
# Use an object to repr. data
# Output obj. as string function
# Set a max depth
# collect records in a seq
# Send to function to tabulate into object (thread)
# fi/fo queue
# Per library output option also...?


proc cmd_read_count*(bamfile: string, positions: string) =
    #[
        Calculates insert size
    ]#
    var 
        b: Bam
        offset: int

    open(b, bamfile, index=true)
    for target in 999915..999930:
        for record in b.query("I", target, target + 1):
                offset = target - record.start.int
                echo target, ": ", record.base_at(offset), " → ", record.base_quality_at(offset), " → ", offset
