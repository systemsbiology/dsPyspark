import collections
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import motifs
from contextlib import contextmanager
from pyspark import SparkContext, SparkConf
from tral.sequence import sequence
import re
import pandas as pd
import os
import argparse
import datetime
import warnings
import sys

SPARK_MASTER='local[*]'
SPARK_APP_NAME='dsSeqModel'
VERSION=0.1

# RUN with:  spark-submit dsSeqModel.py -i <input> -o <output prefix> -s <sample name for vcf>
# this is from streamFastqTransform.py and expects that to output:
# barcode seq1 seq1_quality seq2 seq2_quality read_type read_name
# where barcode combined from seq1 and seq2

# This script uses direct matches - this should be the advantage of using duplex sequencing
# If there is a high error rate for some reason, it will produce few results
# If the person has a mutation in the 12 bp surrounding each CODIS loci, it will fail if you use the reference genome
# It is always best to create a CODIS_surrounding.tsv for an individual person based on sequencing their genome

parser = argparse.ArgumentParser(description="Create alleles from Duplex Sequencing CODIS FASTQ files")
parser.add_argument('-i', '-input',
                    help="Input tab delimited TXT file of Duplex Sequencing information", required=True)
parser.add_argument('-o', '-output',
                    help="Output prefix for TXT and VCF files", required=True)
parser.add_argument('-s', '-sample', help="Sample Name for VCF file", required=True)
parser.add_argument('-nr', '-num_reads', help="Number of reads to count as a Duplex Consensus Sequence",
                    type=int, default=3)
parser.add_argument('-nd', '-num_duplex', help="Number of Duplex Consensus Sequences required to clan",
                    type=int, default=1)
parser.add_argument('-na', '-num_duplex_alleles', help="Number of Duplex Consensus Sequences required to support a CODIS allele report",
                    type=int, default=1)
parser.add_argument('-mc', '-match_count', help="Number of matches required for a sequence to match a CODIS allele",
                    type=int, default=4)
parser.add_argument('-c', '-codis', help="CODIS information for matching loci to genome",
                    default="CODIS_surrounding.tsv")
parser.add_argument('-mp', '-min_partitions', help="Minimum number of partitions to split the text file",
                    type=int, default=20)
parser.add_argument('-t', '-tmp', help="Tmp directory to store files in", default="/tmp")
args = vars(parser.parse_args())
if not os.path.isfile(args['c']):
    raise IOError("Codis file not found: {}".format(args['c']))
if (args['nr'] <= 1):
    warnings.warn("WARNING: Keeping all reads with args['nr'] set to {}.".format(args['nr']), Warning)
if (args['nd'] < 1):
    warnings.warn("WARNING: Invalid option for number of duplex sequences required.", Warning)
    sys.exit()
# test args= {'i': 'CODIS1_spark_input.txt', 's': 'CODIS1', 't': '/40TB_2/workspace/tmp', 'c': 'CODIS_surrounding.tsv', 'o': 'CODIS1_output', 'nd': 1, 'nr': 3, 'mp': 20}

@contextmanager
def spark_manager():
    conf = SparkConf().setMaster(SPARK_MASTER).setAppName(SPARK_APP_NAME)
    conf.set("spark.local.dir", args['t'])
    #conf.set("spark.local.dir", "/40TB_2/workspace/tmp")
    conf.set("spark.driver.memory", "5g")
    spark_context = SparkContext(conf=conf)

    try:
        yield spark_context
    finally:
        spark_context.stop()

with spark_manager() as context:
    dnaRDD = context.textFile(args['i'], minPartitions=args['mp'])
    #dnaArray is
    #['NNNNAAATGCNAGGCGGCCA', 'TTTGGGGGCATCTCTTATACTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTCACCAAGATTG', 'E/EEAAEEEEEEEEEEEEAEEAEEAA<EEAAE<AAE<AEEAAAEEAAAAA<AAAA<A<EAEEE6AAE/<6<E<<AAAAA<EA6<E<<A<A<6AEAEEEAEA<6AAEEE</6AEAE<E/E//<EEAE/AEEE/E<A/6/EEEE//', 'ba', '1', 'NS500773:56:H7LK3AFXX:1:11101:15840:1018']
    # 0 is barcode
    # 1 is sequence
    # 2 is sequence quality
    # 3 is read type (ab or ba indicating whether the barcode sequence was in seq1/seq2 or seq2/seq1)
    # 4 is sequence file (this sequence came from seq1 or seq2)
    # 5 is read name
    dnaArray = dnaRDD.map(lambda line: line.split('\t'))

    def merge_two_dicts(x,y):
        z = x.copy()
        z.update(y)
        return z

    # seqBarcodes is
    #  seqBarcodes.take(1)[0][0] is the barcode and read orientation
    #  seqBarcodes.take(1)[0][1] is the sequence
    seqBarcodes = dnaArray.map(lambda x: ((x[0], x[3]+'-'+x[4]), [x[1]])).reduceByKey(lambda x,y: x+y)
    # seqCombinedBarcodes is the sequenced by barcode and read orientation
    # seqCombinedBarcodes[0] = the barcode
    # seqCombinedBarcodes[1] = a dictionary of sequences by read orientation ex: 'ab-1': [list of seqs]
    seqCombinedBarcodes = seqBarcodes.map(lambda x: [x[0][0], {x[0][1]: x[1]}]).reduceByKey(lambda x,y: merge_two_dicts(x,y))
    # seqCombinedBarcodesFiltered only keeps barcodes that have ab-1, ab-2, ba-1, and ba-2 entries above the cutoff for number of reads (args['nr'])
    seqCombinedBarcodesFiltered = seqCombinedBarcodes.filter(lambda x: len(x[1].get("ab-1", [])) > args['nr'] and len(x[1].get("ab-2",[])) > args['nr']  and len(x[1].get("ba-1",[])) > args['nr'] and len(x[1].get("ba-2",[])) > args['nr'])

    # create consensus sequences for seq using BioPython Seq module
    def create_consensus(seqArray):
        seqDNA   = [Seq(xs, IUPAC.ambiguous_dna) for xs in seqArray]
        seqMotif = motifs.create(seqDNA)
        seqCons  = str(seqMotif.consensus)
        return seqCons

    # combine ab:1 and ba:2 together;  combine ab:2 and ba:1 together
    def overall_consensus(combined):
        barcode, groupedSeqs = combined
        forward = []
        reverse = []
        forward.extend(groupedSeqs["ab-1"])
        reverse.extend(groupedSeqs["ab-2"])
        if len(forward) != len(reverse):
            warnings.warn("WARNING: Nonmatching read counts in overall_consensus for barcode {}.".format(barcode), Warning)
        forward.extend(groupedSeqs["ba-2"])
        reverse.extend(groupedSeqs["ba-1"])
        if len(forward) != len(reverse):
            warnings.warn("WARNING: Nonmatching read counts in overall_consensus for barcode {}.".format(barcode), Warning)
        for_cons = create_consensus(forward)
        if (for_cons != forward[0]):
            print("ERROR: forward doesn't match consensus")
            print(for_cons)
            print(forward[0])
        rev_cons = create_consensus(reverse)
        if (rev_cons != reverse[0]):
            print("ERROR: reverse doesn't match consensus")
            print(rev_cons)
            print(reverse[0])
        return [barcode,  for_cons, rev_cons, len(forward)]

    barcodeConsensus = seqCombinedBarcodesFiltered.map(overall_consensus)
    # break the parts of the sequence up so that we can combine them all together by consensus
    consensusSeqFor = barcodeConsensus.map(lambda x: [x[1], [[x[0], x[3]]]])
    consensusSeqRev = barcodeConsensus.map(lambda x: [x[2], [[x[0], x[3]]]])
    consensusSeqComb = consensusSeqFor.union(consensusSeqRev)

    ##############################################################################################
    # COMBINE BY CONSENSUS RATHER THAN BARCODE
    # consensusSeqComb
    # [['AATATTGGTAATTAAATGTTTACTATAGACTATTTAGTGAGATAAAAAAAAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTTAGTTCGTTCTAAACTATGACAAGTGTTCTATCATACCCTT', ['TACGTTTCAGAATTTGCGAT', 310]]]
    # create
    #[('ACCAGCCTGGGGAAAATAGCAGGACTCCATCTCTACGAAAAATTTGAAAATTAGCCCGGCATGGTGGTACATGCCTGTAGTCCTAGCTATTCAGGAGGCTGAGGCAGGAGGACTGCTTGAGCCCAGGAGTTCGAGGCTGCAGTG', ['TGTAAAACCAAAACTGATAT', 10])]
    consBySequence = consensusSeqComb.reduceByKey(lambda x,y: x+y)

    ###############################################################################################
    # SUM READ COUNTS
    # now that we have all of the sequences matched together, create summary counts
    # consBySequence is [sequence, [list of [barcode, count]]]
    def count_reads(group):
        dcs_count = 0
        dcs_max_reads = 0
        dcs_sum_reads = 0
        # a group is [barcode, count]
        for g in group:
            if (g[1] > dcs_max_reads):
                dcs_max_reads = g[1]
            dcs_count += 1
            dcs_sum_reads += g[1]
        return [dcs_count, dcs_max_reads, dcs_sum_reads]

    consCounted = consBySequence.map(lambda x: [x[0], count_reads(x[1])])

    # Filter for number of DCS required
    consFiltered = consCounted.filter(lambda x: x[1][0] >= args['nd'])

    ###############################################################################################
    # MATCH VERSUS CODIS LOCATIONS

    # load codis information
    codisPickleFile = "codis_regex.p"
    if os.path.isfile(codisPickleFile):
        codis = pd.read_pickle(codisPickleFile)
    else:
        codis = pd.read_csv(args['c'], sep="\t")
        regexer = lambda x: re.compile(x)
        rev_comp = lambda x: str(Seq(x).reverse_complement())
        codis["FULL LEFT COMP"] = codis["FULL LEFT"].apply(rev_comp).apply(regexer)
        codis["FULL RIGHT COMP"] = codis["FULL RIGHT"].apply(rev_comp).apply(regexer)
        codis["ALT LEFT COMP"] = codis["ALT LEFT"].apply(rev_comp).apply(regexer)
        codis["ALT RIGHT COMP"] = codis["ALT RIGHT"].apply(rev_comp).apply(regexer)
        codis["FULL LEFT"] = codis["FULL LEFT"].apply(regexer)
        codis["FULL RIGHT"] = codis["FULL RIGHT"].apply(regexer)
        codis["ALT LEFT"] = codis["ALT LEFT"].apply(regexer)
        codis["ALT RIGHT"] = codis["ALT RIGHT"].apply(regexer)
        codis.to_pickle("codis_regex.p")

    # match vs codis
    def match_codis(codis,seq):
        matches = []
        for i,row in codis.iterrows():
            match_count = 0
            for reg in row[["FULL LEFT","FULL RIGHT", "ALT LEFT","ALT RIGHT","FULL LEFT COMP","FULL RIGHT COMP", "ALT LEFT COMP", "ALT RIGHT COMP"]]:
                print("matching CODIS {} reg {} against seq {}".format(row["ID"], reg, seq))
                if reg.search(seq):
                    match_count = match_count +1
            print("match count {}".format(match_count))
            if match_count >= args['mc']:
                matches.append(row["ID"])
        if len(matches) == 0:
            return ''
        elif len(matches) == 1:
            return matches[0]
        else:
            return 'multiple'

    # codis_matched is the CODIS ID, the second element is the array of tuples of sequence results
    codisCompared = consFiltered.map(lambda x: [match_codis(codis,x[0]), [x]]).reduceByKey(lambda x,y: x+y)
    #codisUnmatched = codisCompared.filter(lambda x: x[0] == '')
    #codisMultiple = codisCompared.filter(lambda x: x[0] == 'multiple')
    codisMatched = codisCompared.filter(lambda x: x[0] != '' and x[0] != 'multiple')
    if (codisMatched.count() == 0):
        print("No results matching CODIS with the num_reads restriction {} and the num_dcs restriction {}".format(args['nr'], args['nd']))
        context.stop()
        sys.exit(1)
    # codis matched is the CODIS loci it matched via sequence primer comparison
    # the second entry is an array of sequence information
    # each sequence is: [sequence, [ dcs_count, max # reads in all dcs, sum of reads in all dcs]]


    ###############################################################################################
    # DETECT REPEATS
    # use tral = uses TRF
    def get_repeats(cr):
        match = cr[0]
        motif = codis.loc[codis["ID"] == match,"motif"].item()
        rc_motif = str(Seq(motif,IUPAC.ambiguous_dna).reverse_complement())
        group = cr[1]
        repeats = []
        for g in group:
            info = [(i.n, i.repeat_region_length, tuple(i.msaD)) for i in sequence.Sequence(g[0], g[0]).detect(denovo=True,detection={"detectors":["TRF"]}).repeats]
            if len(info) > 0:
                g.append(info)
                repeats.append(g)
                # currently if we match against the motif, TRF identifies several repeats with - in the center A-TAG for CSF1PO
                # CSF1PO is detected buggy...
        return repeats

    codisRepeats = codisMatched.map(lambda x: [x[0], get_repeats(x)])
    crDir = 'codis_repeats_'+datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    codisRepeats.saveAsPickleFile(crDir) # we want to save this so that we can avoid re-running the repeat identification software
    codisRepeats = context.pickleFile(crDir)
    # codisRepeats is [codis allele, [list of matches [sequence, [num dcs, max dcs_read, sum_dcs_read], [repeat count, repeat length, repeat seq]]]]

    ###############################################################################################
    # COLLAPSE ALLELES
    # make a key out of the number of repeats to combine alleles together
    # ['CTGGGCTCTTCGTCTTCCGAGTGTTTCTATTTTTAGACCGTTTGGTGTTTGGATAGATAGATAGATAGATAGATATATAAACAAATACTGTTTTGTCTTTCAATGATATCTATCTATCTATCTATCTATCTATCTATCTATATA', [1, 3302, 4037], [(9, 35, ('TATC', 'TATC', 'TATC', 'TATC', 'TATC', 'TATC', 'TATC', 'TATC', 'TAT-'))]]
    def allele_key(codis):
        groups = []
        for g in codis[1]:
            # create an individual entry for the codis and the sequence entry
            c = [codis[0]] + g
            alleles =  ''
            # if it matches > 1 allele we want to discard it
            if len(g[2]) == 1:
                alleles = [(codis[0], g[2][0]), [c]]
            else:
                alleles = [(codis[0], ()), [c]]
            groups.append(alleles)
        return groups
#'codis_repeats_20161217-170655'
    # codisAlleles is [(CODIS ALLELE, (repeat count, repeat length, repeat sequence)), [list of [sequence, [dcs num, dcs read count max, dcs read count sum], [(repeat count, repeat length, repeat sequence)]]]]
    # (('TH01', (8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))), [
    #    ['TH01', 'GGTATCTGGGCTCTGGGGTGATTCCCATTGGCCTGTTCCTCCCTTATTTCCCTCATTCATTCATTCATTCATTCATTCATTCACCATGGAGTCTGTGTTCCCTGTGACCTGCACTCGGAAGCCCTGTGTACAGGGGACTGTGTG', [1, 4977, 4977], [(8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))]],
    #    ['TH01', 'CTTATTTCCCTCATTCATTCATTCATTCATTCATTCATTCACCATGGAGTCTGTGTTCCCTGTGACCTGCACTCGGAAGCCCTGTGTACAGGGGACTGTGTGGGCCAGGCTGGATAATCGGGAGCTTTTCAGCCCACAGGAGGG', [1, 1237, 1237], [(8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))]],
    #    ['TH01', 'TCTGGGGTGATTCCCATTGGCCTGTTCCTCCCTTATTTCCCTCATTCATTCATTCATTCATTCATTCATTCACCATGGAGTCTGTGTTCCCTGTGACCTGCACTCGGAAGCCCTGTGTACAGGGGACTGTGTGGGCCAGGCTGG', [7, 11432, 26282], [(8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))]],
    #    ['TH01', 'CCAGGCTCTAGCAGCAGCTCATGGTGGGGGGTCCTGGGCAAATAGGGGGCAAAATTCAAAGGGTATCTGGGCTCTGGGGTGATTCCCATTGGCCTGTTCCTCCCTTATTTCCCTCATTCATTCATTCATTCATTCATTCATTCA', [1, 133, 133], [(8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))]],
    #    ['TH01', 'TGGGGTGATTCCCATTGGCCTGTTCCTCCCTTATTTCCCTCATTCATTCATTCATTCATTCATTCATTCACCATGGAGTCTGTGTTCCCTGTGACCTGCACTCGGAAGCCCTGTGTACAGGGGACTGTGTGGGCCAGGCTGGAT', [1, 334, 334], [(8, 31, ('TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCAT', 'TCA-'))]]
    # ])
    codisAlleles = codisRepeats.map(lambda x: allele_key(x)).flatMap(lambda x: x).reduceByKey(lambda x,y: x+y)
    # remove the ones that have () as the 'repeat' because they were filtered since they had more than one allele identified in allele_key
    codisAllelesRepeats = codisAlleles.filter(lambda x: len(x[0][1]) == 3)


    ###############################################################################################
    # PRINT OUTPUT
    # output tab file
    # codis_alleles
    # [(('CSF1PO', (2, 31, ('CTATCTATGAATATA', 'CTATCTATCAATATA'))), [['CSF1PO', 'GTGGGCCTCTCTTTTGCTCATGAAGTTTTCAGTGAAATGCATGTGCCTTAATCTCAGCCTTTTTGAAGATCTATCTATGAATATAACTATCTATCAATATATTAATACAAGACTGGGGAACACAGCCAGACACCATGTCTTAAA', [1, 481, 770], [(2, 31, ('CTATCTATGAATATA', 'CTATCTATCAATATA'))]]])]
    # ((codisid, allele), [array of results])
    # array of results is
    # [codis id, sequence, [num dcs, max dcs reads, sum dcs reads], [(repeat num, repeat length, (repeats))]

    # to make the alleles print out near each other, we want to sortByKey
    codisOrdered = codisAllelesRepeats.sortByKey().cache()
    if (codisOrdered.count() == 0):
        print("No results matching CODIS with the num_reads restriction {} and the num_dcs restriction {}".format(args['nr'], args['nd']))
        context.stop()
        sys.exit(1)

    codisOrderedFiltered = codisOrdered.filter(lambda x: len(x[1]) >= args['na'])

    # output the combined list
    def codis_sum(codis_match):
        codis_name = codis_match[0][0]
        codis_alleles = codis_match[0][1] # this is a tuple
        line_fields = ''
        als = [i for e in codis_alleles for i in ([":".join(e)] if isinstance(e,tuple) else [e])]
        for seq in codis_match[1]:
            se = [i for e in seq for i in (e if isinstance(e,list) else [e])][2:-1]
            if line_fields:
                # we want to add 0,2,3 and max 1
                max_ab = max(line_fields[1], se[1])
                lf = [sum(x) for x in zip(line_fields, se)]
                lf[1] = max_ab
                line_fields = lf
            else:
                line_fields = se
        line = [i for e in [codis_name, als, len(codis_match[1]), line_fields] for i in (e if isinstance(e,list) else [e])]
        return line

    codis_sum_out = codisOrderedFiltered.map(lambda x: codis_sum(x))
    sum_filename = args['o'] + '_sum.txt'
    sum_out = open(sum_filename, 'a')
    sum_out.write("#input	" + args['i']+"\n")
    sum_out.write("#output	"+ sum_filename+"\n")
    sum_out.write("#date	"+str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))+"\n")
    pd.DataFrame(codis_sum_out.collect()).to_csv(sum_out, sep="\t", header=["CODIS LOCI", "ALLELE COUNT", "ALLELE LENGTH", "ALLELE", "NUM CLANS", "ALL DCS COUNT", "ALL DCS MAX READS", "ALL DCS SUM READS"], index=False)
    sum_out.close()


    # there is no information on whether this sequence was in seq1 or seq2
    # we're making an output file of:
    # codis_id  allele_count, allele_length, allele_pairs[: delimited], sequence, dcs_count, read_max, read_sum, sscs_count
    # the first line is a key
    # convert each read in the group to a single line and output
    def codis_lines(codis_match):
        codis_name = codis_match[0][0]
        codis_alleles = codis_match[0][1] # this is a tuple
        lines = []
        als = [i for e in codis_alleles for i in ([":".join(e)] if isinstance(e,tuple) else [e])]
        for info in codis_match[1]:
            line = [codis_name, als, info[1], info[2]]
            # remove the inner lists
            l = [i for e in line for i in (e if isinstance(e,list) else [e])]
            lines.append(l)
        return lines

    # the key [0] is a combination of the codis and the allele, but the value [1] is ( codis id, seq, ab array, ba array, [allele arrays]).  So I want to take allele from [0][1] and the rest from value [1]
    codis_out = codisOrderedFiltered.map(lambda x: codis_lines(x)).flatMap(lambda x: x)
    filename = args['o'] + '.txt'
    out = open(filename, 'a')
    out.write("#input	" + args['i']+"\n")
    out.write("#output	"+ filename+"\n")
    out.write("#date	"+str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))+"\n")
    pd.DataFrame(codis_out.collect()).to_csv(out, sep="\t", header=["CODIS LOCI", "ALLELE COUNT", "ALLELE LENGTH", "ALLELE", "SEQ", "DCS COUNT", "DCS MAX READS", "DCS SUM READS"], index=False)
    out.close()

    # VCF
    # must have a header of:
    ##fileformat=VCFv4.2
    ##fileDate=YYYYMMDD
    ##source=DSWF_spark
    ##reference=build37.fa
    ##INFO=<ID=CA,Number=A,Type=String,Description="Codis Allele Sequence">
    ##INFO=<ID=AC,Number=A,Type=String,Description="Codis Allele Count">
    ##INFO=<ID=AL,Number=A,Type=String,Description="Codis Allele Length">
    ##INFO=<ID=DC,Number=1,Type=Integer,Description="DCS Count">
    ##INFO=<ID=DMRC,Number=1,Type=Integer,Description="DCS Max Read Count">
    ##INFO=<ID=DSRC,Number=1,Type=Integer,Description="DCS Sum Read Count">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    # each line is
    # CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  <GENOME ID>
    def vcf_lines(cl):
        codis_id = cl[0]
        codis_chr   = codis[codis["ID"] == codis_id]["chromosome"].item()
        codis_start = codis[codis["ID"] == codis_id]["start"].item()
        codis_motif = codis[codis["ID"] == codis_id]["motif"].item()
        return [codis_chr, codis_start, codis_id, codis_motif, cl[3], cl[5], 'PASS', ";".join(['AC='+str(cl[1]), 'AL='+str(cl[2]), 'CA='+cl[3], 'DC='+str(cl[5]), 'DMRC='+str(cl[6]), 'DSRC='+str(cl[7])]), 'GT', '1']

    vcf_out = codis_out.map(lambda x: vcf_lines(x))

    vcf_filename = args['o'] + '.vcf'
    vout = open(vcf_filename,'a')
    vout.write("##fileformat=VCFv4.2\n")
    vout.write("#fileDate="+str(datetime.datetime.now().strftime("%Y%m%d"))+"\n")
    vout.write("##source=DSWF Pyspark " + args['i']+"\n")
    vout.write("##reference=build37.fa\n")
    vout.write("##INFO=<ID=CA,Number=A,Type=String,Description=\"Codis Allele Sequence\">\n")
    vout.write("##INFO=<ID=AC,Number=A,Type=String,Description=\"Codis Allele Count\">\n")
    vout.write("##INFO=<ID=AL,Number=A,Type=String,Description=\"Codis Allele Length\">\n")
    vout.write("##INFO=<ID=DC,Number=1,Type=Integer,Description=\"DCS Count\">\n")
    vout.write("##INFO=<ID=DMRC,Number=1,Type=Integer,Description=\"DCS Max Read Count\">\n")
    vout.write("##INFO=<ID=DSRC,Number=1,Type=Integer,Description=\"DCS Sum Read Count\">\n")
    vout.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Codis Genotype\"\n")
    vout.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{}\n".format(args['s']))
    pd.DataFrame(vcf_out.collect()).to_csv(vout, sep="\t", header=False, index=False)
    vout.close()

    print("Finished Building consensus")
