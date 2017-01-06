# Duplex Sequence Pyspark Processing Pipeline

This takes two paired end FASTQ files and analyses them for Duplex Sequencing amplicons then matches valied Duplex 
Sequencing Duplex Consensus Sequences against CODIS sequences.

1) python3 streamFastqTransform.py -s1 seq1.fastq -s2 seq2.fastq -o sample_spark_input.txt

2) spark-submit dsSeqModel.py -i sample_spark_input.txt -o output_prefix -s sample_name

Creates three output files:
1) output_prefix.txt - A list of all of the DCS and the counts of sequences that support it with the CODIS match
2) output_prefix.vcf - A VCF format list of all DCS and counts.
3) output_prefix_sum.txt - A list of all clans that support a CODIS allele (or a matching sequence in the default
CODIS_surrounding.tsv or a similar format file containing other sequences supplied with the -c flag).

Example
The repository provides example CODIS_spark_input.txt and CODIS_surrounding.tsv that can be used to run an 
example analysis with:

spark-submit dsSeqModel.py -i examples/CODIS_spark_input.txt -o CODIS_output -s CODIS -c examples/CODIS_surrounding.tsv

Note: D21S11 is too long for a 156 bp read, so it requires a match count (mc) setting of 2 to match CODIS allele 
and it will output an incorrect allele count because the read cannot span the full allele sequence.

Requires Pyspark >= 2.x, TRAL, and TRF

# Installation

Install Spark
https://spark.apache.org/downloads.html

Install TRAL
http://elkeschaper.github.io/tral/
> pip install TRAL

Install TRF
http://tandem.bu.edu/trf/trf.html



# Workflow
1) use streamFastqTransform.py in parallel to upload FASTQ files to a HDFS

2) use dsSeqModel.py to analyze files

A Docker repository will be available soon.
