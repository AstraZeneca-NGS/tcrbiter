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
	fh = open(Filename)
	fhw = open(os.path.splitext(Filename)[0] + ".csv", "w")
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

######MAIN CODE##################
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
prereqs = ["intersectBlastmerging.R", "mixcrFiltering.R",
           "myFields.alignmentExport.txt", bedfile]
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
results = os.path.join(tcr_outdir, readpairkey + ".results.csv")
filteredresults = os.path.join(tcr_outdir, readpairkey + ".filteredResults.csv")
cline = ("Rscript --vanilla %s/mixcrFiltering.R %s %s"
         %(scriptfolder, results, filteredresults))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "MixcrFiltering.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "MixcrFiltering.stdout.txt"), "w") as oh:
    oh.write(pout)
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
cline = ('Rscript --vanilla %s/intersectBlastmerging.R %s %s %s %s'
         %(scriptfolder, read1bedintersect, read2bedintersect,
           read1blastclean, read2blastclean))
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
with open(os.path.join(tcr_errordir, "intersectBlastMerge.error.txt"), "w") as oh:
    oh.write(perr)
with open(os.path.join(tcr_errordir, "intersectBlastMerge.stdout.txt"), "w") as oh:
    oh.write(pout)
logger.info("Finished merging intersected BLAST results for %s." % readpairkey)

logger.info("Started creating final report for %s." % readpairkey)

logger.info("Finished creating final report for %s." % readpairkey)

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

# estimate T-cells found
intersect_stats = os.path.join(tcr_outdir, readpairkey + ".intersectStats.txt")
with open(intersect_stats) as ih:
    ih.next()
    statsline = ih.next()
    tokens = [x.strip() for x in statsline.split()]
    samplename = os.path.basename(tokens[0])
    vdjhits = int(tokens[2])
    falsehits = int(tokens[3]) + int(tokens[4])

# output report
reportfile = os.path.join(tcr_outdir, readpairkey + ".report.txt")
with open(reportfile, "w") as oh:
    oh.write("# samplename: %s\n" % samplename)
    oh.write("# inputfile: %s\n" % inputfile)
    oh.write("# totalreads: %s\n" % totalreads)
    oh.write("# TCRreads: %s\n" % vdjhits)
    oh.write("# falsehits: %s\n" % falsehits)

# VDJ intersections
vdj_intersection = os.path.join(tcr_outdir, readpairkey + ".intersectVDJ.txt")
with open(vdj_intersection) as ih, open(reportfile, "a") as oh:
    ih.next()
    oh.write("readid gene1 gene2\n")
    for line in ih:
        tokens = [x.strip() for x in line.split()]
        qid = tokens[0]
        outstr = " ".join([qid, tokens[4], tokens[26]]) + "\n"
        oh.write(outstr)
logger.info("Processing complete for %s" % readpairkey)
