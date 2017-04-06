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
from operator import itemgetter
import logging

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

#Retrieve fastq files from the input
if args.r1:
	rp1 = args.r1
if args.r2:
	rp2= args.r2
#retreive blast db
if args.blastdb:
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
if not os.path.exists(tcr_outdir):
	os.makedirs(tcr_outdir)
if not os.path.exists(tcr_errordir):
	os.makedirs(tcr_errordir)
if not os.path.exists(tcr_errordir):
	os.makedirs(tcr_errordir)

# setup command logger
commandlogfile = os.path.join(outdir, readpairkey, "tcrOutput", "tcrCommands.txt")
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

biterdir_exists_msg = ("The results folder already exists, this may "
                     "happen if you've already run the script, however it may "
                     "be due to the fastq files being named without an _R1 in "
                     "them. Consider renaming your fastq files if this becomes "
                     "an issue.")
biterdir = os.path.join(outdir, readpairkey)
if not os.path.exists(biterdir):
	os.mkdir(biterdir)
else:
    logger.warning(biterdir_exists_msg)

tcr_errordir = os.path.join(biterdir, "tcrError")
tcr_outputdir = os.path.join(biterdir, "tcrOutput")
if not os.path.exists(tcr_errordir):
	os.mkdir(tcr_errordir)
if not os.path.exists(tcr_outputdir):
	os.mkdir(tcr_outputdir)


logger.info("Starting analysis.")

# mixcr alignments to T cell receptor beta sequences
logger.info("Starting mixcr align for %s." % readpairkey)
cline = ("mixcr align -l TRB -OvParameters.geneFeatureToAlign=Vgene -s hsa "
         "--save-description --report %s/%s/tcrOutput/%s.alignmentReport.log "
         "%s %s %s/%s/tcrOutput/%s.vdjca"
         %(outdir, readpairkey, readpairkey, rp1, rp2, outdir,
           readpairkey, readpairkey))
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
cline = "mixcr exportAlignments -s -pf %s/myFields.alignmentExport.txt %s/%s/tcrOutput/%s.vdjca %s/%s/tcrOutput/%s.results.txt" %(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
csvout = os.path.join(tcr_outputdir, readpairkey + ".results.txt")
convertToCsv(csvout, tcr_outputdir)
cmdlogger.info("convertToCsv (internal to script)")
logger.info("Through conversion to CSV for %s." % readpairkey)

logger.info("Starting filtering of mixcr results for %s." % readpairkey)
cline = "Rscript --vanilla %s/mixcrFiltering.R %s/%s/tcrOutput/%s.results.csv %s/%s/tcrOutput/%s.filteredResults.csv"%(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
fhw=open(('%s/%s/tcrOutput/%s.filteredread1.readkey.txt' %(outdir, readpairkey, readpairkey)), "w")
fhfasta = open("%s/%s/tcrOutput/%s.read1.filtered.fasta"%(outdir, readpairkey, readpairkey), "w")
cline = "gzip -cd %s |grep -A1 -Ff %s/%s/tcrOutput/%s.filteredread1.txt" % (rp1, outdir, readpairkey, readpairkey )
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
for line in pout.split('\n'):
	if "@" in line:
		readid = line.lstrip("@").strip()
	elif not line.strip() == '' and not "--" in line:
		fhw.write("%s\t%s\n"%(readid, len(line.strip())))
		fhfasta.write(">%s\n%s\n"%(readid, line))
fhw.close()
fhfasta.close()
logger.info("Through fastq for Read 1 of the pair for %s." % readpairkey)
# read2
fhw=open(('%s/%s/tcrOutput/%s.filteredread2.readkey.txt'%(outdir, readpairkey, readpairkey)), "w")
fhfasta = open("%s/%s/tcrOutput/%s.read2.filtered.fasta"%(outdir, readpairkey, readpairkey), "w")
cline = "gzip -cd %s |grep -A1 -Ff %s/%s/tcrOutput/%s.filteredread2.txt" % (rp2, outdir, readpairkey, readpairkey)
cmdlogger.info(cline)
child = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
pout, perr = child.communicate()
rc = child.returncode
if child.returncode:
    logger.error("Error running %s." % cline)
    sys.exit(1)
for line in pout.split('\n'):
	if "@" in line:
		readid = line.lstrip("@").strip()
	elif not line.strip() == '' and not "--" in line:
		fhw.write("%s\t%s\n"%(readid, len(line.strip())))
		fhfasta.write(">%s\n%s\n"%(readid, line))
fhw.close()
fhfasta.close()
logger.info("Through fastq for Read 2 of the pair. Fastq keys have been output "
            "and fasta files created for each of the filtered sets.")

logger.info("Starting BLAST confirmation for %s." % readpairkey)
cline = 'blastn -db %s -query %s/%s/tcrOutput/%s.read1.filtered.fasta -outfmt 6 -num_threads 6 > %s/%s/tcrOutput/%s.read1.blast.txt'%(blastdb, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
cline = 'blastn -db %s -query %s/%s/tcrOutput/%s.read2.filtered.fasta -outfmt 6 -num_threads 6 > %s/%s/tcrOutput/%s.read2.blast.txt'%(blastdb, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
cline = 'grep -vwE "HLA" %s/%s/tcrOutput/%s.read1.blast.txt > %s/%s/tcrOutput/%s.read1.blastClean.txt'%(outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
cline = 'grep -vwE "HLA" %s/%s/tcrOutput/%s.read2.blast.txt > %s/%s/tcrOutput/%s.read2.blastClean.txt'%(outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
blast2bed('%s/%s/tcrOutput/%s.read1.blastClean.txt'%(outdir, readpairkey, readpairkey), '%s/%s/tcrOutput/%s.read1.bed'%(outdir, readpairkey, readpairkey))
blast2bed('%s/%s/tcrOutput/%s.read2.blastClean.txt'%(outdir, readpairkey, readpairkey), '%s/%s/tcrOutput/%s.read2.bed'%(outdir, readpairkey, readpairkey))
logger.info("Completed conversion of BLAST hits to BED for %s." % readpairkey)

logger.info("Starting intersecting BED files for %s." % readpairkey)
cline = 'bedtools intersect -a %s -b %s/%s/tcrOutput/%s.read1.bed -wo > %s/%s/tcrOutput/%s.read1.intersectBed.txt'%(bedfile, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
cline = 'bedtools intersect -a %s -b %s/%s/tcrOutput/%s.read2.bed -wo > %s/%s/tcrOutput/%s.read2.intersectBed.txt'%(bedfile, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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

logger.info("Starting finalizing results for %s." % readpairkey)
cline = 'Rscript --vanilla %s/intersectBlastmerging.R %s/%s/tcrOutput/%s.read1.intersectBed.txt %s/%s/tcrOutput/%s.read2.intersectBed.txt %s/%s/tcrOutput/%s.read1.blastClean.txt %s/%s/tcrOutput/%s.read2.blastClean.txt' %(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
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
logger.info("Finished finalizing results for %s." % readpairkey)
logger.info("Process complete.")
