# This benchmarking script is use to benchmark trained nanopore models
# bash BenchmarkModel.sh <model_path> <test_data_path> <reference_file> <groundtruth_file>
dorado basecaller $1 --trim all --kit-name SQK-NBD114-24 $2 > basecalls.bam
samtools fastq basecalls.bam > basecalls.fq
# length filter
fastplong --disable_adapter_trimming --disable_quality_filtering --length_required 600 --length_limit 800 --thread 32 --in basecalls.fq --out basecalls_600_800.fq
# quality filter
fastplong --disable_adapter_trimming -m 15 --disable_length_filtering --thread 32 --in basecalls_600_800.fq --out basecalls_q15.fq
# 
