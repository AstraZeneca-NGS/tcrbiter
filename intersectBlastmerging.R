#!/usr/bin/env Rscript
#Merging and analyzing results of posthoc blast and intersectBed from filtered reads
#Author: Lara Lewis McGrath
#Date: January 29, 2016

#Usage: Rscript --vanilla script.R intersect.r1.txt intersect.r2.txt blast.r1.txt blast.r2.txt

args= commandArgs(trailingOnly=TRUE)

if(length(args)!=4){
  stop("Incorrect arguments. Please input 4 files. \n", call.=FALSE)
} else if (length(args)==4){
  print("Starting analysis.")
}

filename <- strsplit(args[1],".read1.intersectBed.txt",fixed=TRUE)
filename <- sapply(filename,"[",1)

#Install necessary packages
library("data.table")

##Make formatting functions

formatBLAST <- function(y){
  colnames(y)<- c("qseqid","sseqid","pident","length","mismatch","gapopen",
                  "qstart","qend","sstart","send","evalue","bitscore") 
  return(y)}

formatINTERSECT <- function(x){
  colnames(x)<- c("targetchr","targetstart","targetend","gene","querychr",
                  "querystart","queryend","qseqid","alignmentlength") 
  return(x)}

whichgene <- function(w) {
  ifelse((grepl("TRBV",w$gene)=="TRUE"),paste("V"),paste("J"))
}


#Load data into dataframe

intersect.r1 <- formatINTERSECT(read.table(args[1],header=F, stringsAsFactors=FALSE))
intersect.r2 <- formatINTERSECT(read.table(args[2],header=F, stringsAsFactors=FALSE))

blast.r1 <- formatBLAST(read.table(args[3],header=F, stringsAsFactors=FALSE))
blast.r2 <- formatBLAST(read.table(args[4],header=F, stringsAsFactors=FALSE))

print("Files loaded, beginning formatting and filtering...")

#How many chromosomes did each read map to?
blast.r1 <- data.table(blast.r1)
blast.r1 <- blast.r1[, blastchrcount := length(unique(sseqid)), by=qseqid]
blast.r1 <- as.data.frame(blast.r1)
blast.r2 <- data.table(blast.r2)
blast.r2 <- blast.r2[, blastchrcount := length(unique(sseqid)), by=qseqid]
blast.r2 <- as.data.frame(blast.r2)

r1 <- merge(intersect.r1,blast.r1,by="qseqid",all.x=T,all.y=F)
r1$read <- "1"
r2 <- merge(intersect.r2,blast.r2,by="qseqid",all.x=T,all.y=F)
r2$read <- "2"

#Combine reads 1 and 2
results <- rbind(r1,r2)
results$type  <- whichgene(results)


#Take all V hits and reduce to longest one per read
withv <- subset(results,grepl("TRBV",results$gene))
withv <- withv[order(-withv$alignmentlength),]
withv <- withv[!duplicated(withv[c("qseqid")]),]

#Take all J hits and reduce to longest one per read
withj <- subset(results,grepl("TRBJ",results$gene))
withj <- withj[order(-withj$alignmentlength),]
withj <- withj[!duplicated(withj[c("qseqid")]),]

#Combine reads with v and j
merged <- merge(withv,withj,by="qseqid",all=TRUE)
merged$sameread <- merged$read.x==merged$read.y
#Remove any reads with NA results 
final <- merged[complete.cases(merged),]

#Write files
write.table(merged,file=paste(filename,".intersectMerged.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
write.table(final,file=paste(filename,".intersectVDJ.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

#Run descriptive statistics on merged reads
stats <- data.frame(
  sample = c(paste(filename)),
  total = c(nrow(merged)),
  VDJ = c(nrow(subset(merged,merged$type.x=="V" & merged$type.y=="J" | merged$type.x=="J" & merged$type.y=="V"))),
  VV = c(nrow(subset(merged,merged$type.x=="V" & merged$type.y=="V"))),
  JJ = c(nrow(subset(merged,merged$type.x=="J" & merged$type.y=="J"))),
  VNA = c(nrow(subset(merged,is.na(merged$type.x) & merged$type.y=="V"|merged$type.x=="V" & is.na(merged$type.y)))),
  JNA = c(nrow(subset(merged,is.na(merged$type.x) & merged$type.y=="J"|merged$type.x=="J" & is.na(merged$type.y)))),
  singleread = c(nrow(subset(merged,merged$read.x==merged$read.y))),
  acrossreads = c(nrow(subset(merged,merged$read.x!=merged$read.y)))
)
if (rowSums(stats[,c(3:7)])==sum(stats[2])){
  print ("Descriptive stats agree and are complete.")
}
#Write files
write.table(stats,file=paste(filename,".intersectStats.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


print("Finished.")