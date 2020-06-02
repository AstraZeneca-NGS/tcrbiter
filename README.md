# tcrBiter usage and output
Author: Lara McGrath
Date: 8/29/2017 (last update: 6/2/2020)

## Requirements
tcrBiter is a python script that can be run from anywhere but must be in the same folder as the included bed files (`bed/TRB-hg38.bed` and `myFields.alignmentExport.txt`). 

Additional Requirements: 
* `Python` (if using Python v2, use `tcrBiter_betav2.py`; if using Python v3, use `tcrBiter_betav3.py`)
* `mixcr` (v1.7)
* `BLAST` and `bedtools` must be installed and in the path correctly
* must have BLAST database of hg38 or create one using `makeblastdb`.

## Usage
```
python /path/to/tcrBiter_betav3.py --r1 /path/to/some_R1.fastq.gz --r2 /path/to/some_R2.fastq.gz --blastdb /path/to/BLAST/db/for/hg38
```

## Output
tcrBiter creates two directories: `tcrError/` and `tcrOutput/`

The following files are made in `tcrError/`:
```
bedtoolsintersect.read1.error.txt   blastcleanread1.error.txt   blastread1.error.txt   intersectBlastMerge.error.txt   MixcrExport.error.txt
bedtoolsintersect.read1.stdout.txt  blastcleanread1.stdout.txt  blastread1.stdout.txt  intersectBlastMerge.stdout.txt  MixcrExport.stdout.txt
bedtoolsintersect.read2.error.txt   blastcleanread2.error.txt   blastread2.error.txt   MixcrAlign.error.txt            MixcrFiltering.error.txt
bedtoolsintersect.read2.stdout.txt  blastcleanread2.stdout.txt  blastread2.stdout.txt  MixcrAlign.stdout.txt           MixcrFiltering.stdout.txt
```

The following files are made in `tcrOutput/`:
```
filename.alignmentReport.log        filename.filteredStats.pdf     filename.read1.filtered.fasta    filename.results.csv
filename.description.txt            filename.intersectMerged.txt   filename.read1.intersectBed.txt  filename.results.curated.csv
filename.filteredread1.readkey.txt  filename.intersectStats.txt    filename.read2.bed               filename.results.txt
filename.filteredread1.txt          filename.intersectVDJ.txt      filename.read2.blastClean.txt    filename.statsFiltered.csv
filename.filteredread2.readkey.txt  filename.read1.bed             filename.read2.blast.txt         filename.vdjca
filename.filteredread2.txt          filename.read1.blastClean.txt  filename.read2.filtered.fasta    tcrCommands.txt
filename.filteredResults.csv        filename.read1.blast.txt       filename.read2.intersectBed.txt
```

Important results files are:
* `filename.intersectStats.txt`: Defines number of read pairs in: (1) total, (2) with VDJ rearrangment, (3) with only V genes, (4) with only J genes, (5) with one V gene, (6) with one J gene, (7) rearrangement found in single read, (8) rearrangement split across pair
* `filename.intersectVDJ.txt`: Breaks down features of each rearranged read

Other results files include:
* `filename.alignmentReport.log`: MiXCR alignment output report
* `filename.description.txt`: Information on how many reads were aligned, filtered, and had defined CDR3 regions
* `filename.filteredStats.pdf`: Displays histograms of features of aligned reads



	

