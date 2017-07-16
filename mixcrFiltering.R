#!/usr/bin/env Rscript

##For command line interpretation of MiXCR output
#by Lara Lewis McGrath, January 22, 2016
#Edited 3/3/2016: added writing of curated output file
#Edited 4/1/2016: edited to include additional column

#Input is results.csv from MiXCR alignment

args= commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("Atleast one argument must be supplied. \n", call.=FALSE)
} else if (length(args)==1){
  args[2]=paste(sapply((strsplit(args[1],".results.csv",fixed=TRUE)),"[",1),".filtered.csv",sep="")
}

filename <- strsplit(args[1],".results.csv",fixed=TRUE)
filename <- sapply(filename,"[",1)

print(paste("Loading ", args[1], sep=""))
#Load data into dataframe
df <- read.csv(args[1], header=T,stringsAsFactors=FALSE)
print("File loaded, beginning formatting and filtering...")



##Make function for formatting input
###Create column header function
formatPE<- function(y){
  colnames(y)<- c("readId","read1","read2","sequence","targets","CDR3","CDR3sequence","bestVhit",
                  "allVhits","Vstart","Vend","Vgenelength","Vquerystart","Vqueryend",
                  "Vmut","Vscore","bestJhit","allJhits","Jstart","Jend","Jgenelength",
                  "Jquerystart","Jqueryend","Jmut", "Jscore")
  y$alignmentlength <- nchar(y$sequence)
  y$Jcliplength <- y$Jend-y$Jstart
  y$Vcliplength <- y$Vend-y$Vstart
  y$Jcoverage <- (y$Jcliplength/y$Jgenelength)
  y$Vcoverage <- y$Vcliplength/y$Vgenelength
  y$Vdistancefromend <- (y$Vgenelength-y$Vend)
  y$Jdistancealignmentend <- abs(y$Jqueryend-y$alignmentlength)
  y <- y[c("readId","read1","read2","sequence","alignmentlength","targets","CDR3","bestVhit",
           "allVhits","Vstart","Vend","Vcliplength","Vgenelength","Vquerystart","Vqueryend",
           "Vcoverage","Vdistancefromend","Vmut","Vscore","bestJhit","allJhits","Jstart","Jend","Jcliplength","Jgenelength",
           "Jquerystart","Jqueryend","Jcoverage","Jdistancealignmentend","Jmut", "Jscore")]

  return(y)
}

#Create function with filter parameters
filter <- function(z){
  z <- subset(z,(z$CDR3!="NA" & z$CDR3 < 60 & z$CDR3 > 20) | z$Vdistancefromend < 20 & z$Jstart < 15 & z$Vquerystart <20)
  z <- z[order(-z$Jcoverage),]
  return(z)
}


#push data through functions for formatting and filtering
df <- formatPE(df)
df_out <- filter(df)

#Run descriptive statistics on filtered reads
category <- c("mean","min","max","median","var","sd")
mean <-lapply(df_out,mean,na.rm=T)
min <-lapply(df_out,min,na.rm=T)
max <- lapply(df_out,max,na.rm=T)
median <- lapply(df_out,median,na.rm=T)
var <- lapply(df_out,var,na.rm=T)
sd <- lapply(df_out,sd,na.rm=T)
stats <- as.data.frame(rbind(mean, min,max,median,var,sd))
stats <- cbind(category,stats)
stats <- subset(stats, select=c("category","alignmentlength","CDR3","Vstart","Vend",
                                "Vcliplength","Vgenelength","Vquerystart","Vqueryend",
                                "Vcoverage","Vdistancefromend","Jstart","Jcliplength",
                                "Jgenelength","Jquerystart","Jqueryend","Jcoverage"))
stats_out<- data.frame(lapply(stats,as.character),stringsAsFactors=F)


#Create description of aligned file and new filtered file
library(plyr)
aligned <- length(df$CDR3)
filtered <- length(df_out$CDR3)
CDR3 <- sum(!is.na(df_out$CDR3))
CDR3isNA <- sum(is.na(df_out$CDR3))
meanCDR3 <- round(mean(df_out$CDR3,na.rm=T),0)
minCDR3 <- min(df_out$CDR3,na.rm=T)
maxCDR3 <- max(df_out$CDR3,na.rm=T)
desc <- data.frame(rbind(aligned,filtered,CDR3,CDR3isNA,meanCDR3,minCDR3,maxCDR3))
colnames(desc)<- filename

#Create graphs for visualizing data
if(nrow(df_out) > 0)  {
  pdf(paste(filename,".filteredStats.pdf",sep=""),width=7,height=11)
  par(mfrow=c(3,2))
  par(mar=c(5,4,7,2))
  hist(df_out$Jstart, main="Location of J rearrangement",xlab="Distance (bp)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  hist(df_out$Vdistancefromend, main="Location of V rearrangement",xlab="Distance (bp)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  if(is.finite(minCDR3)) {
      hist(df_out$CDR3, main="CDR3",xlab="Length (bp)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  }
  hist(df_out$Jcoverage, main="J alignment",xlab="Coverage of germline sequence (%)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  hist(df_out$alignmentlength, main="Alignment Length",xlab="Length (bp)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  hist(df_out$Vcliplength, main="V alignment",xlab="Length (bp)",cex.main=1,cex.lab=0.9, col="darkmagenta")
  mtext(paste(filename," - Statistics on Filtered Reads",sep=""),side=3,padj=2,outer=TRUE)
  dev.off()
  }



#write output
write.table(df_out,file=args[2],row.names=FALSE,sep=",",quote=FALSE)
write.table(stats_out,file=paste(filename,".statsFiltered.csv",sep=""),row.names=FALSE,sep=",",quote=FALSE)
write.table(desc,file=paste(filename,".description.txt",sep=""),row.names=TRUE,sep="\t",quote=FALSE)
write.table(df_out$read1,file=paste(filename,".filteredread1.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(df_out$read2,file=paste(filename,".filteredread2.txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
write.table(df,file=paste(filename,".results.curated.csv",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")

print("Finished.")
