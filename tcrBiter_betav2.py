#! python

#TCRbiter: Detecting T cell rearrangements in silico from non-targeted DNA Sequencing (WGS/WES)
#Authors: Tristan Lubinski, Lara McGrath
#Date: Sept 12, 2016
#Purpose: Detect and count reads supporting T cell rearrangements 

#Built for MiXCR v1.7
#For download, use, and license information for MiXCR visit https://github.com/milaboratory/mixcr

import sys, copy, csv, os, glob, subprocess, argparse
from operator import itemgetter

#print "This script can be run from anywhere but must be in the same folder as the files: myFields.alignmentExport.txt, mixcrFiltering.R, TRBsequences.bed, and intersectBlastmerging.R"
print "Additional Requirements: Python, mixcr, BLAST, R, and bedtools must be installed and in the path correctly."
print "Additional Requirements: must have BLAST database of hg38 or create one using makeblastdb."
parser = argparse.ArgumentParser(prog='tcrBiter', usage='python /path/to/%(prog)s.py --r1 path/to/some_R1.fastq.gz --r2 path/to/some_r2.fastq.gz --blastdb /path/to/blast/db/for/hg38')
parser.add_argument('--build', help='human build to use', default='hg38')
parser.add_argument('--r1', help='Read 1 of gzipped Fastq pair', required=True)
parser.add_argument('--r2', help='Read 2 of gzipped Fastq pair', required=True)
parser.add_argument('--blastdb', help="This is the BLAST database made for hg38 (if you don't have one please run makeblastdb and create it first)", required=True)

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
	fhw = open(outfolder+Filename.split('/')[-1].rstrip('txt')+"csv",'w')
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
args=parser.parse_args()
# figure out if the requirements are actually installed
if not which("mixcr"):
    print "mixcr not found, please install it or put it in your path."
    sys.exit(1)
if not which("blastn"):
    print "blastn not found, please install it or put it in your path."
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
outdir = os.getcwd()
scriptfolder = "/".join(sys.argv[0].split('/')[:-1])
scriptfolder = sys.path[0]
bedfile = os.path.join(scriptfolder, "bed", "TRB-" + args.build + ".bed")
if not glob.glob('%s/intersectBlastmerging.R'%(scriptfolder)):
	print "Missing intersectBlastmerging.R (this should have been included in the git package, please try cloning the package again if you continue to get this error)"
if not os.path.exists(bedfile):
	print "Missing the TRB sequences bedfile (this should have been included in the git package, please try cloning the package again if you continue to get this error)"
if not glob.glob('%s/mixcrFiltering.R'%(scriptfolder)):
	print "Missing mixcrFiltering.R(this should have been included in the git package, please try cloning the package again if you continue to get this error)"
if not glob.glob('%s/myFields.alignmentExport.txt'%(scriptfolder)):
	print "Missing myFields.alignmentExport.txt(this should have been included in the git package, please try cloning the package again if you continue to get this error)"
#next check if folder structure for output is there, if not make it
if not os.path.exists(outdir+"/"+readpairkey):
	os.mkdir(outdir+"/"+readpairkey)
else:
	print "Warning: The results folder already exists, this may happen if you've already run the script, however it may be due to the fastq files being named without an _R1 in them. Consider renaming your fastq files if this becomes an issue."
if not os.path.exists(outdir+"/"+readpairkey+"/tcrError"):
	os.mkdir(outdir+"/"+readpairkey+"/tcrError")
if not os.path.exists(outdir+"/"+readpairkey+"/tcrOutput"):
	os.mkdir(outdir+"/"+readpairkey+"/tcrOutput")
#record the commands
fhwc = open(outdir+"/"+readpairkey+"/tcrOutput/tcrCommands.txt","w")
print "==================\n"
print "Starting analysis..."
#begin commands
#mixcr alignments to T cell receptor beta sequences
cline = "mixcr align -l TRB -OvParameters.geneFeatureToAlign=Vgene -s hsa --save-description --report %s/%s/tcrOutput/%s.alignmentReport.log %s %s %s/%s/tcrOutput/%s.vdjca"%(outdir, readpairkey, readpairkey, rp1, rp2, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrAlign.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrAlign.stdout.txt",'w')
fhw.write(pout)
fhw.close()

print "Through mixcr align for "+readpairkey
#export mixcr alignments
cline = "mixcr exportAlignments -s -pf %s/myFields.alignmentExport.txt %s/%s/tcrOutput/%s.vdjca %s/%s/tcrOutput/%s.results.txt" %(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey) 
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrExport.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrExport.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print "Through mixcr ExportAlign for "+readpairkey
#process the file
convertToCsv(outdir+"/"+readpairkey+"/tcrOutput/"+readpairkey+".results.txt", outdir+"/"+readpairkey+"/tcrOutput/")
fhwc.write("convertToCsv (internal to script)\n\n")
print "Through conversion to CSV for "+readpairkey
cline = "Rscript --vanilla %s/mixcrFiltering.R %s/%s/tcrOutput/%s.results.csv %s/%s/tcrOutput/%s.filteredResults.csv"%(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrFiltering.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/MixcrFiltering.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print "Completed filtering Rscript for "+readpairkey
## now we're going to parse the output into sets and filter the fastq files
fhw=open(('%s/%s/tcrOutput/%s.filteredread1.readkey.txt' %(outdir, readpairkey, readpairkey)), "w")
fhfasta = open("%s/%s/tcrOutput/%s.read1.filtered.fasta"%(outdir, readpairkey, readpairkey), "w")
cline = "zcat %s |grep -A1 -Ff %s/%s/tcrOutput/%s.filteredread1.txt" % (rp1, outdir, readpairkey, readpairkey )
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
for line in pout.split('\n'):
	if "@" in line:	
		readid = line.lstrip("@").strip()
	elif not line.strip() == '' and not "--" in line:
		fhw.write("%s\t%s\n"%(readid, len(line.strip())))
		fhfasta.write(">%s\n%s\n"%(readid, line))
fhw.close()
fhfasta.close()
print "Through fastq for Read 1 of the pair"
#repeat for read 2
fhw=open(('%s/%s/tcrOutput/%s.filteredread2.readkey.txt'%(outdir, readpairkey, readpairkey)), "w")
fhfasta = open("%s/%s/tcrOutput/%s.read2.filtered.fasta"%(outdir, readpairkey, readpairkey), "w")
cline = "zcat %s |grep -A1 -Ff %s/%s/tcrOutput/%s.filteredread2.txt" % (rp2, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
for line in pout.split('\n'):
	if "@" in line:	
		readid = line.lstrip("@").strip()
	elif not line.strip() == '' and not "--" in line:
		fhw.write("%s\t%s\n"%(readid, len(line.strip())))
		fhfasta.write(">%s\n%s\n"%(readid, line))
fhw.close()
fhfasta.close()	
print "Through fastq for Read 2 of the pair\nFastq keys have been output and fasta files created for each of the filtered sets\n\n==================\n"
print "Starting BLAST confirmation..."
#now perform BLAST to verify alignments
cline = 'blastn -db %s -query %s/%s/tcrOutput/%s.read1.filtered.fasta -outfmt 6 -num_threads 6 > %s/%s/tcrOutput/%s.read1.blast.txt'%(blastdb, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastread1.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastread1.stdout.txt",'w')
fhw.write(pout)
fhw.close()
cline = 'blastn -db %s -query %s/%s/tcrOutput/%s.read2.filtered.fasta -outfmt 6 -num_threads 6 > %s/%s/tcrOutput/%s.read2.blast.txt'%(blastdb, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastread2.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastread2.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print 'Blast completed.'
#Remove all hits with HLA as it prohibits uploading in R later
cline = 'grep -vwE "HLA" %s/%s/tcrOutput/%s.read1.blast.txt > %s/%s/tcrOutput/%s.read1.blastClean.txt'%(outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastcleanread1.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastcleanread1.stdout.txt",'w')
fhw.write(pout)
fhw.close()
cline = 'grep -vwE "HLA" %s/%s/tcrOutput/%s.read2.blast.txt > %s/%s/tcrOutput/%s.read2.blastClean.txt'%(outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastcleanread2.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/blastcleanread2.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print "Cleaned blast results."
#now modify this to a bed file
blast2bed('%s/%s/tcrOutput/%s.read1.blastClean.txt'%(outdir, readpairkey, readpairkey), '%s/%s/tcrOutput/%s.read1.bed'%(outdir, readpairkey, readpairkey))
blast2bed('%s/%s/tcrOutput/%s.read2.blastClean.txt'%(outdir, readpairkey, readpairkey), '%s/%s/tcrOutput/%s.read2.bed'%(outdir, readpairkey, readpairkey))
print "Bed conversion complete."
#intersect the bedfiles
cline = 'bedtools intersect -a %s -b %s/%s/tcrOutput/%s.read1.bed -wo > %s/%s/tcrOutput/%s.read1.intersectBed.txt'%(bedfile, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/bedtoolsintersect.read1.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/bedtoolsintersect.read1.stdout.txt",'w')
fhw.write(pout)
fhw.close()
cline = 'bedtools intersect -a %s -b %s/%s/tcrOutput/%s.read2.bed -wo > %s/%s/tcrOutput/%s.read2.intersectBed.txt'%(bedfile, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/bedtoolsintersect.read2.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/bedtoolsintersect.read2.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print "Bedfile intersection completed."
#finalize the results
cline = 'Rscript --vanilla %s/intersectBlastmerging.R %s/%s/tcrOutput/%s.read1.intersectBed.txt %s/%s/tcrOutput/%s.read2.intersectBed.txt %s/%s/tcrOutput/%s.read1.blastClean.txt %s/%s/tcrOutput/%s.read2.blastClean.txt' %(scriptfolder, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey, outdir, readpairkey, readpairkey)
fhwc.write(cline+"\n\n")
pout, perr = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True ).communicate()
fhw = open(outdir+"/"+readpairkey+"/tcrError/intersectBlastMerge.error.txt",'w')
fhw.write(perr)
fhw.close()
fhw = open(outdir+"/"+readpairkey+"/tcrError/intersectBlastMerge.stdout.txt",'w')
fhw.write(pout)
fhw.close()
print "Process complete." 
fhwc.close()
