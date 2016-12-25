Duplex Sequence Pyspark Processing Pipeline

This takes two paired end FASTQ files and analyses them for Duplex Sequencing amplicons then matches valied Duplex 
Sequencing Duplex Consensus Sequences against CODIS sequences.

1) python3 streamFastqTransform.py -s1 <seq1.fastq> -s2 <seq2.fastq> -o <sample>_spark_input.txt

2) spark-submit dsSeqModel.py -i CODIS_spark_input.txt -o CODIS_output -s CODIS

Creates three output files:
1) CODIS_output.txt - A list of all of the DCS and the counts of sequences that support it with the CODIS match
2) CODIS_output.vcf - A VCF format list of all DCS and counts.
3) CODIS_output_sum.txt - A list of all clans that support a CODIS allele.

Note: D21S11 is too long for a 156 bp read, so it requires a match count (mc) setting of 2 to match CODIS allele 
and it will output an incorrect allele count because the read cannot span the full allele sequence.

Requires Pyspark >= 2.x

Workflow
1) use streamFastqTransform.py in parallel to upload FASTQ files to a HDFS
2) use dsSeqModel.py to analyze files

A Docker repository will be available soon.
