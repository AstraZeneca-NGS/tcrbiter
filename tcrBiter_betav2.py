#! python

# TCRbiter: Detecting T cell rearrangements in silico from non-targeted DNA
# Sequencing (WGS/WES)
# Authors: Tristan Lubinski, Lara McGrath
# Date: Sept 12, 2016
# Purpose: Detect and count reads supporting T cell rearrangements

# Built for MiXCR v1.7
# For download, use, and license information for MiXCR visit
# https://github.com/milaboratory/mixcr

import sys, copy, csv, os, glob, subprocess, argparse
import gzip
from operator import itemgetter
import logging
import itertools
from io import BufferedReader
import re

import pandas as pd
import numpy as np

__version__ = "0.0.1"

# setup main progress logger
FORMAT = '%(asctime)s - %(levelname)s - %(message)s'
formatter = logging.Formatter(FORMAT)
logger = logging.getLogger("TCRbiter")
streamhandler = logging.StreamHandler(sys.stderr)
streamhandler.setLevel(logging.INFO)
streamhandler.setFormatter(formatter)
logger.addHandler(streamhandler)
logger.setLevel(logging.INFO)

def which(program):
    """ returns the path to an executable or None if it can't be found"""
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def is_mixcr2(cmd):
    cmd += " --version"
    child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True)
    pout, perr = child.communicate()
    rc = child.returncode
    if "MiXCR v2" in perr:
        return True
    else:
        return False

######DEFINED FUNCTIONS##################
def convertToCsv(Filename, outfolder):
    #Formatting of MiXCR export alignments output
    #Purpose to modify the input txt file into a csv with necessary formatting
    header = ["readId", "read1", "read2", "sequence", "targets", "CDR3",
              "CDR3sequence", "bestVhit", "allVhits", "Vstart", "Vend",
              "Vgenelength", "Vquerystart", "Vqueryend", "Vmut", "Vscore",
              "bestJhit", "allJhits", "Jstart", "Jend", "Jgenelength",
              "Jquerystart", "Jqueryend", "Jmut", "Jscore"]
    fh = open(Filename)
    fhw = open(os.path.splitext(Filename)[0] + ".csv", "w")
    fhw.write(",".join(header) + "\n")
    fh.readline() #remove header
    commathings=[7,8,10,11]
    for line in fh:
        line2 = line.split('\t')
        line3= line.split()
        if len(line2) == 11:
            line2.insert(5,' ')
        for place in commathings:
            line2[place]= line2[place].replace('"','').replace(',',' ')
        if int(line2[4])>1:
            allseq = line2[3].replace('"','').replace(',','')
            line2[3] = allseq
            allscore = line2[9].replace('"','').split(',')
            if len(allscore)>1: #have a line to split
                maxscore = 0.00
                maxset =''
                for item in allscore:
                    if "|" in item:
                        newscore = float(item.split('|')[-1])
                        if newscore > maxscore:
                            maxscore = newscore
                            maxset = item
                line2[9]=maxset #put the new seq in there
            allscore = line2[-1].replace('"','').split(',')
            if len(allscore)>1: #have a line to split
                maxscore = 0.00
                maxset =''
                for item in allscore:
                    if "|" in item:
                        newscore = float(item.strip().split('|')[-1])
                        if newscore > maxscore:
                            maxscore = newscore
                            maxset = item.strip()
                line2[-1]=maxset+"\n" #put the new seq in there
        finalline = '\t'.join(line2).replace(",","").replace("|",",").replace("\t",",")
        fhw.write(finalline)
    fh.close()
    fhw.close()

def formatPE(filename):
    df = pd.read_csv(filename)
    df["alignmentlength"] = df["sequence"].apply(len)
    df["Jcliplength"] = df["Jend"] - df["Jstart"]
    df["Vcliplength"] = df["Vend"] - df["Vstart"]
    df["Jcoverage"] = df["Jcliplength"]/df["Jgenelength"]
    df["Vcoverage"] = df["Vcliplength"]/df["Vgenelength"]
    df["Vdistancefromend"] = df["Vgenelength"] - df["Vend"]
    df["Jdistancealignmentend"] = abs(df["Jqueryend"]-df["alignmentlength"])
    order = ["readId", "read1", "read2", "sequence", "alignmentlength",
             "targets", "CDR3", "bestVhit", "allVhits", "Vstart", "Vend",
             "Vcliplength", "Vgenelength", "Vquerystart", "Vqueryend",
             "Vcoverage", "Vdistancefromend", "Vmut", "Vscore", "bestJhit",
             "allJhits", "Jstart", "Jend", "Jcliplength", "Jgenelength",
             "Jquerystart", "Jqueryend", "Jcoverage", "Jdistancealignmentend",
             "Jmut", "Jscore"]
    return(df[order])

def filterPE(df, CDR3high, CDR3low, Vdistancefromend, Jstart, Vquerystart):
    keepcdr3 = ((df["CDR3"] < CDR3high) &
                (df["CDR3"] > CDR3low))
    keepvj = ((df["Vdistancefromend"] < Vdistancefromend) &
              (df["Jstart"] < Jstart) &
              (df["Vquerystart"] < Vquerystart))
    return(df[keepcdr3 | keepvj])

def filtered_read_statistics(df):
    numerics = df.select_dtypes(include=['float64', 'int64'])
    dfo = pd.DataFrame.from_dict({
        "mean": numerics.apply(np.nanmean),
        "min": numerics.apply(np.min),
        "max": numerics.apply(np.max),
        "median": numerics.apply(np.nanmedian),
        "var": numerics.apply(lambda x: np.nanvar(x, ddof=1)),
        "sd": numerics.apply(lambda x: np.nanstd(x, ddof=1))
    })
    dfo = dfo.transpose()
    dfo["category"] = dfo.index
    keep = ["category", "alignmentlength", "CDR3", "Vstart", "Vend",
            "Vcliplength", "Vgenelength", "Vquerystart", "Vqueryend",
            "Vcoverage", "Vdistancefromend", "Jstart", "Jcliplength",
            "Jgenelength", "Jquerystart", "Jqueryend", "Jcoverage"]
    dfo = dfo[keep]
    dfo = dfo.reindex(["mean", "min", "max", "median", "var", "sd"])
    return(dfo)

def raw_vs_filtered(df, df_filt):
    dfo = pd.DataFrame.from_dict({
        "aligned": [len(df.index)],
        "filtered": [len(df_filt.index)],
        "CDR3": [df_filt["CDR3"].count()],
        "meanCDR3": [np.mean(df_filt["CDR3"])],
        "minCDR3": [np.min(df_filt["CDR3"])],
        "maxCDR3": [np.max(df_filt["CDR3"])]
    })
    dfo["CDR3isNA"] = dfo["filtered"] - dfo["CDR3"]
    return(dfo)

def blast2bed(blastfile,bedout):
	#takes in the cleaned blast file and makes a bed file out of it
	fh=open(blastfile)
	fhw=open(bedout,'w')
	for line in fh:
		line2 = line.split('\t')
		lstart = int(line2[8])
		lstop = int(line2[9])
		fhw.write("%s\t%s\t%s\t%s\n"%(line2[1],min(lstart,lstop),max(lstart,lstop),line2[0]))
	fh.close()
	fhw.close()

def stream_fastq(file_handler):
    """
    Generator which gives all four lines if a fastq read as one string
    borrowed from https://github.com/vals/umis
    """
    next_element = ''
    for i, line in enumerate(file_handler):
        next_element += line
        if i % 4 == 3:
            yield next_element
            next_element = ''

def read_fastq(filename):
    """
    return a stream of FASTQ entries, handling gzipped and empty files
    borrowed from https://github.com/vals/umis
    """
    if not filename:
        return itertools.cycle((None,))
    if filename.endswith('gz'):
        filename_fh = BufferedReader(gzip.open(filename, mode='rt'))
    else:
        filename_fh = open(filename)
    return stream_fastq(filename_fh)

def create_filtered_set(filename):
    with open(filename) as in_handle:
        return {x.split()[0] for x in in_handle}

def load_BLAST(filename):
    cnames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    return(pd.read_csv(filename, sep="\t", names=cnames))

def load_intersect_bed(filename):
    cnames = ["targetchr", "targetstart", "targetend", "gene", "querychr",
              "querystart", "queryend", "qseqid", "alignmentlength"]
    return(pd.read_csv(filename, sep="\t", names=cnames))

def detect_subtype(symbol):
    if symbol.startswith("TRBV"):
        return "V"
    elif symbol.startswith("TRBJ"):
        return "J"
    else:
        return None

###k###MAIN CODE##################
# This script can be run from anywhere but must be in the same folder as
# the files: myFields.alignmentExport.txt, mixcrFiltering.R, TRBsequences.bed,
# and intersectBlastmerging.R

# Additional Requirements: Python, mixcr, BLAST, R, and bedtools must be
# installed and in the path correctly. A BLAST database of hg38 or hg19
# must be created with makeblastdb.
usage = ("python /path/to/%(prog)s.py --r1 path/to/some_R1.fastq.gz --r2 "
         "path/to/some_r2.fastq.gz --blastdb /path/to/blast/db/for/hg38")
blastdbhelp = ("This is the BLAST database made for hg38 (if you don't have "
               "one please run makeblastdb and create it first")


parser = argparse.ArgumentParser(prog='TCRbiter', usage=usage)
parser.add_argument("--Jstart", default=15,
                   help=("Drop J reads which align to a J target gene, but after "
                         "this base in that gene."))
parser.add_argument("--Vquerystart", default=20,
                   help=("Drop V reads which align to a V target gene, but after "
                         "this base in the query sequence (the sequencing read)."))
parser.add_argument("--Vdistancefromend", default=20,
                   help=("Drop V reads which align to a V target gene, but the "
                         "end of the alignment is further than this to the end "
                         "of the V target gene."))
parser.add_argument("--CDR3high", default=60,
                   help=("Drop reads which have a longer CDR3 than this."))
parser.add_argument("--CDR3low", default=20,
                   help=("Drop reads which have a shorter CDR3 than this."))
parser.add_argument('--version', help="TCRbiter version", action='version',
                    version = __version__)
parser.add_argument('--build', help='human build to use', default='hg38')
parser.add_argument('--r1', help='Read 1 of gzipped Fastq pair', required=True)
parser.add_argument('--r2', help='Read 2 of gzipped Fastq pair', required=True)
parser.add_argument('--blastdb', help=blastdbhelp, required=True)
args = parser.parse_args()

# figure out if the requirements are actually installed
if not which("mixcr"):
    print "mixcr not found, please install it or put it in your path."
    sys.exit(1)
if not which("blastn"):
    print "blastn not found, please install it or put it in your path."
    sys.exit(1)
if not which("bedtools"):
    print "bedtools not found, please install it or put it in your path."
    sys.exit(1)

rp1 = args.r1
rp2 = args.r2
blastdb = args.blastdb
if "_R1" in rp1:
	readpairkey = rp1.split('/')[-1].split('_R1')[0]
else:
	readpairkey = rp1.split('/')[-1].split('_1')[0]

# setup analysis directories
outdir = os.getcwd()
tcr_outdir = os.path.join(outdir, readpairkey, "tcrOutput")
tcr_errordir = os.path.join(outdir, readpairkey, "tcrError")
if not os.path.exists(tcr_outdir):
	os.makedirs(tcr_outdir)
if not os.path.exists(tcr_errordir):
	os.makedirs(tcr_errordir)

# setup command logger
commandlogfile = os.path.join(tcr_outdir, "tcrCommands.txt")
cmdlogger = logging.getLogger("TCRbiter-commands")
cmdloggerhandler = logging.FileHandler(commandlogfile)
cmdloggerhandler.setFormatter(logging.Formatter("%(message)s"))
cmdloggerhandler.setLevel(logging.INFO)
cmdlogger.addHandler(cmdloggerhandler)
cmdlogger.setLevel(logging.INFO)

scriptfolder = sys.path[0]
bedfile = os.path.join(scriptfolder, "bed", "TRB-" + args.build + ".bed")
prereqs = ["mixcrFiltering.R", "myFields.alignmentExport.txt", bedfile]
missing_prereq_msg = ("Missing {} (this should have been included in the git "
                      "package, please try cloning the package again if you "
                      "continue to get this error).")
for fn in prereqs:
    prereq = os.path.join(scriptfolder, fn)
    if not os.path.exists(prereq):
        logger.error(missing_prereq_msg.format(prereq))
        sys.exit(1)

logger.info("Starting analysis.")

# mixcr alignments to T cell receptor beta sequences
logger.info("Starting mixcr align for %s." % readpairkey)
alignment_report = os.path.join(tcr_outdir, "%s.alignmentReport.log" % readpairkey)
chaincmd = " -c TRB " if is_mixcr2("mixcr") else " -l TRB "
outfile = os.path.join(tcr_outdir, readpairkey + ".vdjca")
cline = ("mixcr align %s -OvParameters.geneFeatureToAlign=Vgene -s hsa "
         "--save-description --report %s %s %s %s"
         %(chaincmd, alignment_report, rp1, rp2, outfile))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "MixcrAlign.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "MixcrAlign.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("Completed mixcr align for %s." % readpairkey)

# export mixcr alignments
logger.info("Starting mixcr exportAlignments for %s." % readpairkey)
vdjca = os.path.join(tcr_outdir, readpairkey + ".vdjca")
results = os.path.join(tcr_outdir, readpairkey + ".results.txt")
cline = ("mixcr exportAlignments -s -pf %s/myFields.alignmentExport.txt %s %s"
         %(scriptfolder, vdjca, results))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "MixcrExport.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "MixcrExport.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("Completed mixcr exportAlignments for %s." % readpairkey)

# convert to CSV process the file
logger.info("Starting conversion of mixcr alignments to CSV for %s." % readpairkey)
csvout = os.path.join(tcr_outdir, readpairkey + ".results.txt")
convertToCsv(csvout, tcr_outdir)
cmdlogger.info("convertToCsv (internal to script)")
logger.info("Through conversion to CSV for %s." % readpairkey)

logger.info("Starting filtering of mixcr results for %s." % readpairkey)
stem = os.path.join(tcr_outdir, readpairkey)
results = os.path.join(tcr_outdir, readpairkey + ".results.csv")
filteredresults = os.path.join(tcr_outdir, readpairkey + ".filteredResults.csv")
description = os.path.join(tcr_outdir, readpairkey + ".description.txt")
mixcr = formatPE(results)
mixcr_filt = filterPE(mixcr, args.CDR3high, args.CDR3low, args.Vdistancefromend,
                      args.Jstart, args.Vquerystart)
stats = filtered_read_statistics(mixcr_filt)
descriptive = raw_vs_filtered(mixcr, mixcr_filt)

mixcr_filt.to_csv(filteredresults, index=False)
stats.to_csv(stem + ".statsFiltered.csv", index=False)
descriptive.to_csv(stem + ".description.txt", index=False, sep="\t")
mixcr_filt["read1"].to_csv(stem + ".filteredread1.txt", index=False)
mixcr_filt["read2"].to_csv(stem + ".filteredread2.txt", index=False)
mixcr.to_csv(stem + ".results.curated.csv", index=False)
logger.info("Completed filtering Rscript for %s." % readpairkey)

# now we're going to parse the output into sets and filter the fastq files
# read1
logger.info("Starting filtering of the FASTQ files for %s." % readpairkey)
filteredreadfile = os.path.join(tcr_outdir, readpairkey + ".filteredread1.txt")
fhfastafile = os.path.join(tcr_outdir, readpairkey + ".read1.filtered.fasta")
filtered_set = create_filtered_set(filteredreadfile)
parser_re = re.compile('@(?P<name>.*) .*\\n(?P<seq>.*)\\n.*\\n.*\\n')
with open(fhfastafile, "w") as fasta_out:
    for read in read_fastq(rp1):
        match = parser_re.match(read)
        name = match.group('name').split()[0]
        seq = match.group('seq')
        if name not in filtered_set:
            continue
        fasta_out.write(">" + name + "\n")
        fasta_out.write(seq + "\n")
logger.info("Through fastq for Read 1 of the pair for %s." % readpairkey)
# read2
logger.info("Starting filtering of the FASTQ files for %s." % readpairkey)
filteredreadfile = os.path.join(tcr_outdir, readpairkey + ".filteredread2.txt")
fhfastafile = os.path.join(tcr_outdir, readpairkey + ".read2.filtered.fasta")
filtered_set = create_filtered_set(filteredreadfile)
parser_re = re.compile('@(?P<name>.*) .*\\n(?P<seq>.*)\\n.*\\n.*\\n')
with open(fhfastafile, "w") as fasta_out:
    for read in read_fastq(rp2):
        match = parser_re.match(read)
        name = match.group('name').split()[0]
        seq = match.group('seq')
        if name not in filtered_set:
            continue
        fasta_out.write(">" + name + "\n")
        fasta_out.write(seq + "\n")
logger.info("Through fastq for Read 2 of the pair. Fastq keys have been output "
            "and fasta files created for each of the filtered sets.")

logger.info("Starting BLAST confirmation for %s." % readpairkey)
read1fasta = os.path.join(tcr_outdir, readpairkey + ".read1.filtered.fasta")
blastoutput = os.path.join(tcr_outdir, readpairkey + ".read1.blast.txt")
cline = ('blastn -db %s -query %s -outfmt 6 -num_threads 6 > %s'
         %(blastdb, read1fasta, blastoutput))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "blastread1.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "blastread1.stdout.txt"), "w") as oh:
    oh.write(pout)
read2fasta = os.path.join(tcr_outdir, readpairkey + ".read2.filtered.fasta")
blastoutput = os.path.join(tcr_outdir, readpairkey + ".read2.blast.txt")
cline = ('blastn -db %s -query %s -outfmt 6 -num_threads 6 > %s'
         %(blastdb, read2fasta, blastoutput))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "blastread2.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "blastread2.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("BLAST confirmation completed for %s." % readpairkey)

# Remove all hits with HLA as it prohibits uploading in R later
logger.info("Filtering hits with HLA for %s." % readpairkey)
read1blast = os.path.join(tcr_outdir, readpairkey + ".read1.blast.txt")
read1filteredblast = os.path.join(tcr_outdir, readpairkey + ".read1.blastClean.txt")
cline = 'grep -vwE "HLA" %s > %s' %(read1blast, read1filteredblast)
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "blastcleanread1.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "blastcleanread1.stdout.txt"), "w") as oh:
    oh.write(pout)

read2blast = os.path.join(tcr_outdir, readpairkey + ".read2.blast.txt")
read2filteredblast = os.path.join(tcr_outdir, readpairkey + ".read2.blastClean.txt")
cline = 'grep -vwE "HLA" %s > %s' %(read2blast, read2filteredblast)
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "blastcleanread2.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "blastcleanread2.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("Completed filtering BLAST results for %s." % readpairkey)
#now modify this to a bed file

logger.info("Starting conversion of BLAST hits to BED format for %s." % readpairkey)
read1filteredblast = os.path.join(tcr_outdir,
                                  readpairkey + ".read1.blastClean.txt")
read1bed = os.path.join(tcr_outdir, readpairkey + ".read1.bed")
blast2bed(read1filteredblast, read1bed)
read2filteredblast = os.path.join(tcr_outdir,
                                  readpairkey + ".read2.blastClean.txt")
read2bed = os.path.join(tcr_outdir, readpairkey + ".read2.bed")
blast2bed(read2filteredblast, read2bed)
logger.info("Completed conversion of BLAST hits to BED for %s." % readpairkey)

read1bed = os.path.join(tcr_outdir, readpairkey + ".read1.bed")
read1bedintersect = os.path.join(tcr_outdir,
                                 readpairkey + ".read1.intersectBed.txt")
logger.info("Starting intersecting BED files for %s." % readpairkey)
cline = ('bedtools intersect -a %s -b %s -wo > %s'
         %(bedfile, read1bed, read1bedintersect))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "bedtoolsintersect.read1.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "bedtoolsintersect.read1.stdout.txt"), "w") as oh:
    oh.write(pout)
read2bed = os.path.join(tcr_outdir, readpairkey + ".read2.bed")
read2bedintersect = os.path.join(tcr_outdir,
                                 readpairkey + ".read2.intersectBed.txt")
logger.info("Starting intersecting BED files for %s." % readpairkey)
cline = ('bedtools intersect -a %s -b %s -wo > %s'
         %(bedfile, read2bed, read2bedintersect))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "bedtoolsintersect.read2.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "bedtoolsintersect.read2.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("Completed BEDfile intersection for %s." % readpairkey)

logger.info("Starting merging intersected BLAST results for %s." % readpairkey)
read1bedintersect = os.path.join(tcr_outdir,
                                 readpairkey + ".read1.intersectBed.txt")
read2bedintersect = os.path.join(tcr_outdir,
                                 readpairkey + ".read2.intersectBed.txt")
read1blastclean = os.path.join(tcr_outdir,
                               readpairkey + ".read1.blastClean.txt")
read2blastclean = os.path.join(tcr_outdir,
                               readpairkey + ".read2.blastClean.txt")
r1_intersect_bed = load_intersect_bed(read1bedintersect)
r2_intersect_bed = load_intersect_bed(read2bedintersect)
r1_blastclean = load_BLAST(read1blastclean)
r2_blastclean = load_BLAST(read2blastclean)
sseqid_per_read1 = r1_blastclean[["qseqid", "sseqid"]].drop_duplicates().groupby(["qseqid"]).count()
sseqid_per_read1.columns = ["blastchrcount"]
sseqid_per_read1["read"] = "1"
sseqid_per_read1["qseqid"] = sseqid_per_read1.index
r1join = pd.merge(r1_intersect_bed, sseqid_per_read1, on="qseqid")
r1 = pd.merge(r1join, r1_blastclean)

sseqid_per_read2 = r2_blastclean[["qseqid", "sseqid"]].drop_duplicates().groupby(["qseqid"]).count()
sseqid_per_read2.columns = ["blastchrcount"]
sseqid_per_read2["read"] = "2"
sseqid_per_read2["qseqid"] = sseqid_per_read2.index
r2join = pd.merge(r2_intersect_bed, sseqid_per_read2, on="qseqid")
r2 = pd.merge(r2join, r2_blastclean)

results = r1.append(r2)
results["type"] = results["gene"].apply(detect_subtype)
results = results.sort_values(["alignmentlength", "evalue"], ascending=[False, True]).drop_duplicates(["qseqid", "type", "read"])
results.to_csv(os.path.join(tcr_outdir,
                            readpairkey + ".intersectVDJ.txt"), sep="\t", index=False)

stats = results.groupby("qseqid").type.agg(lambda x: "".join(sorted(set(x))))
rstats = results.groupby("qseqid").read.agg(lambda x: len(set(x)))

stats = {"sample": readpairkey,
         "total": stats.shape[0],
         "VDJ": sum(stats == "JV"),
         "VV": sum(stats == "VV"),
         "JJ": sum(stats == "JJ"),
         "VNA": sum(stats == "V"),
         "JNA": sum(stats == "J"),
         "acrossreads": sum(rstats == 2),
         "singleread": sum(rstats == 1)}
pd.DataFrame(stats, index=[1]).to_csv(os.path.join(tcr_outdir,
                                                   readpairkey + ".intersectStats.txt"),
                                      sep="\t", index=False)
logger.info("Finished merging intersected BLAST results for %s." % readpairkey)

logger.info("Started creating final report for %s." % readpairkey)
# estimate genome size
with open(blastdb + ".fai") as ih:
    # don't doublte-count the alts
    chromlines = [x for x in ih if "_alt" not in x]
    chromsizes = map(int, [x.split()[1].strip() for x in chromlines if x.isdigit()])
genomesize = sum(chromsizes)

# estimate total alignments
with open(alignment_report) as ih:
    for line in ih:
        if line.startswith("Total sequencing"):
            totalreadsline = line
        if line.startswith("Input file"):
            inputfileline = line
totalreads = int(totalreadsline.split(":")[1].strip())
inputfile = inputfileline.split(":")[1].strip()

# output report
reportfile = os.path.join(tcr_outdir, readpairkey + ".report.txt")
with open(reportfile, "w") as oh:
    oh.write("# samplename: %s\n" % stats["sample"])
    oh.write("# inputfile: %s\n" % inputfile)
    oh.write("# totalreads: %s\n" % totalreads)
    oh.write("# TCRreads: %s\n" % stats["VDJ"])
    oh.write("# falsehits: %s\n" % (stats["VNA"] + stats["JNA"]))

    with open(reportfile, "a") as oh:
        oh.write("readid gene1 gene2\n")
        resreport = results.groupby("qseqid").gene.agg(lambda x: " ".join(set(x)))
        for qname, genes in zip(resreport.index, resreport):
            oh.write(" ".join([qname, genes]) + "\n")

logger.info("Finished creating final report for %s." % readpairkey)
logger.info("Processing complete for %s" % readpairkey)
