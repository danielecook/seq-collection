#!/usr/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

PARENT_DIR=$(git rev-parse --show-toplevel)
export PATH="${PATH}:${PARENT_DIR}"

tear_down() {
    rm test.fq.gz
}

#===================#
# Utility Functions #
#===================#
printf_stdin() { local stdin; read -d '' -u 0 stdin; printf "$@" "$stdin"; }
export -f printf_stdin

echo "STDIN"
cat tests/fastq/dup.fq | ${PARENT_DIR}/sc --debug fq-dedup &
echo "END"

echo "AS FILE"
${PARENT_DIR}/sc --debug fq-dedup ${PARENT_DIR}/tests/fastq/dup.fq
echo "END2"

#======#
# json #
#======#

set -o nounset
run test_json sc json "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal \"X\" "$(jq '.CHROM' "${STDOUT_FILE}")"
assert_equal 17276844 "$(jq '.POS' "${STDOUT_FILE}")"
assert_equal "\"PASS\"" "$(jq '.FILTER[0]' "${STDOUT_FILE}" )"
assert_equal 999 "$(jq '.QUAL' "${STDOUT_FILE}")"
assert_equal \"T\" "$(jq '.REF' "${STDOUT_FILE}")"

run test_json_pretty sc json --pretty "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal 13 "$(wc -l "${STDOUT_FILE}" | awk '{ print $1 }' )"

# INFO
run single_info_item sc json --info="DP" "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal 9836 "$(jq '.INFO.DP' "${STDOUT_FILE}")"

run multi_info_item sc json --info="DP,MQ,DP4,HOB,INDEL" "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal 9836 "$(jq '.INFO.DP' "${STDOUT_FILE}")"
assert_equal 60 "$(jq '.INFO.MQ' "${STDOUT_FILE}")"
assert_equal 92 "$(jq '.INFO.DP4[3]' "${STDOUT_FILE}")"
assert_equal 0.5 "$(jq '.INFO.HOB' "${STDOUT_FILE}")"

# FORMAT
run format_dp sc json --format="DP" "${PARENT_DIR}/tests/vcf/test.vcf.gz" I:41947-41947 | jq '.FORMAT.DP|add'
assert_equal 2094 "$(jq '.FORMAT.DP|add' "${STDOUT_FILE}")"

# GENOTYPE
run gt sc json -f GT "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal 0 "$(jq '.FORMAT.GT[0][0]' "${STDOUT_FILE}")"
run gt_all sc json -f ALL "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844
assert_equal 0 "$(jq '.FORMAT.GT[0][0]' "${STDOUT_FILE}")"

# CONVERT MISSING VALUES TO NULL
run missing_float_to_null sc json -f PL "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844 
assert_equal "null" "$(jq -c '.FORMAT.PL[0]' "${STDOUT_FILE}")"
run missing_floats_to_null sc json -f PL "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17276844-17276844 
assert_equal "[null,null]" "$(jq -c '.FORMAT.PL[0:2]' "${STDOUT_FILE}")"

#BCSQ
run bcsq_gene sc json -i BCSQ -n "${PARENT_DIR}/tests/vcf/test.bcsq.vcf.gz" chr22:40679539-40679539 
assert_equal \"MCHR1\" "$(jq -c '.INFO.BCSQ[0].gene' "${STDOUT_FILE}")"

# ARRAY
run array_trailing_comma sc json -a "${PARENT_DIR}/tests/vcf/test.vcf.gz" X:17261695-17276844
assert_equal , "$(tail -n 3 "${STDOUT_FILE}" | head -n1 | awk '{print substr($0,length,1)}')"
assert_equal \} "$(tail -n 2 "${STDOUT_FILE}" | head -n1 | awk '{print substr($0,length,1)}')"


#==========#
# fq-dedup #
#==========#

run fq_dedup sc fq-dedup tests/fastq/dup.fq
assert_equal 4 "$(grep -c '@' "${STDOUT_FILE}")"
assert_equal 4 "$(grep -c '@' "${STDOUT_FILE}")"

run fq_dedup_gz sc fq-dedup tests/fastq/dup.fq.gz
assert_equal 4 "$(grep -c '@' "${STDOUT_FILE}")"
assert_equal 4 "$(grep -c '@' "${STDOUT_FILE}")"

#=========#
# fq-meta #
#=========#
# Download parse script if not exists
if ! test -s "${PARENT_DIR}/illumina_instrument.py"; then
    curl https://raw.githubusercontent.com/10XGenomics/supernova/master/tenkit/lib/python/tenkit/illumina_instrument.py | \
            grep -v 'import martian' | \
            sed 's/martian.exit/exit/g' | \
            sed 's/f.readline()/f.readline().decode("utf-8")/g' | \
            sed 's/exit/print/g' > "${PARENT_DIR}/illumina_instrument.py"

    # Fix python script
    echo "Fixing ${PARENT_DIR}/illumina_instrument.py"
    2to3 --write --nobackups "${PARENT_DIR}/illumina_instrument.py"
fi;

function test_fq() {
    gzip -c $1 > test.fq.gz
    python -c "import illumina_instrument as ia; print(ia.sequencer_detection_message(['test.fq.gz']))"
}

# Genome Analyzer II
run fq_meta_genome_analyzer1 sc fq-meta tests/fastq/illumina_1.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" \
             "$(test_fq tests/fastq/illumina_1.fq | grep -E -o "Genome Analyzer IIx" | uniq | sed 's/ //g')"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_1.fq | grep -E -o "likely" | uniq)"

run fq_meta_genome_analyzer2 sc fq-meta tests/fastq/illumina_2.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" \
             "$(test_fq tests/fastq/illumina_2.fq | grep -E -o "Genome Analyzer IIx" | uniq | sed 's/ //g')"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_2.fq | grep -E -o "likely" | uniq)"

# Unresolvable; This FASTQ is from wikipedia...
run fq_meta3 sc fq-meta tests/fastq/illumina_3.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" ""
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" ""
assert_equal "$(test_fq tests/fastq/illumina_3.fq | grep -E -o 'could not detect')" "could not detect"

run fq_meta4 sc fq-meta tests/fastq/illumina_4.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}" )" ""
assert_equal "$(cut -f 3  "${STDOUT_FILE}"| cut -f 1 -d ':')" ""
assert_equal "$(test_fq tests/fastq/illumina_4.fq | grep -E -o 'could not detect')" "could not detect"

run fq_meta_hiseq_2000_2500 sc fq-meta tests/fastq/illumina_2000_2500.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" "HiSeq2000/2500"
assert_equal "$(cut -f 3 "${STDOUT_FILE}")" "high:machine+flowcell"
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" \
             "$(test_fq tests/fastq/illumina_2000_2500.fq | grep -E -o "HiSeq2000/2500" | uniq)"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_2000_2500.fq | grep -E -o 'high' )"

# HiSeq 4000
run fq_meta_hiseq_4000 sc fq-meta tests/fastq/illumina_3000_4000.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}" )" "HiSeq3000/4000"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" )" "high:machine+flowcell"
assert_equal "$(cut -f 2 "${STDOUT_FILE}" )" \
             "$(test_fq tests/fastq/illumina_3000_4000.fq | grep -E -o "HiSeq3000/4000" | uniq)"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_3000_4000.fq | grep -E -o 'high')"

# HiSeq X
run fq_meta_hiseq_x sc fq-meta tests/fastq/illumina_hiseq_x.fq
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" "HiSeqX"
assert_equal "$(cut -f 3 "${STDOUT_FILE}")" "high:machine+flowcell"
assert_equal "$(cut -f 3 "${STDOUT_FILE}" | cut -f 1 -d ':')" \
             "$(test_fq tests/fastq/illumina_hiseq_x.fq | grep -E -o 'high')"

# Novaseq
run fq_meta_novaseq sc fq-meta tests/fastq/novaseq.fq
assert_equal "$(cut -f 3 "${STDOUT_FILE}")" "high:machine+flowcell"
assert_equal "$(cut -f 2 "${STDOUT_FILE}")" "NovaSeq"

#=============#
# Insert size #
#=============#

run insert_size sc insert-size tests/bam/test.bam
assert_equal "$(cut -f 1 "${STDOUT_FILE}" | tail -n 1)" "179"
assert_equal "$(cut -f 2 "${STDOUT_FILE}" | printf_stdin "%4.1f" )" "176.5"
assert_equal "$(cut -f 4 "${STDOUT_FILE}" | tail -n 1)" "38"
assert_equal "$(cut -f 5 "${STDOUT_FILE}" | tail -n 1)" "358"
assert_equal "$(cut -f 6 "${STDOUT_FILE}" | tail -n 1)" "359"

#======#
# iter #
#======#

run sci_format sc iter tests/vcf/test.vcf.gz 1e6
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I:1-1000000"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA:1-13794"

run comma_iter sc iter tests/vcf/test.vcf.gz 100,000
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I:1-100000"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA:1-13794"

run vcf_iter sc iter tests/vcf/test.vcf.gz 100000
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I:1-100000"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA:1-13794"

run bam_iter sc iter tests/bam/test.bam 1000000
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I:0-999999"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA:0-13793"

run iter_chrom_vcf sc iter tests/vcf/test.vcf.gz 0
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA"

run iter_chrom_bam sc iter tests/bam/elegans.bam 0
assert_equal "$(head -n 1 "${STDOUT_FILE}")" "I"
assert_equal "$(tail -n 1 "${STDOUT_FILE}")" "MtDNA"
