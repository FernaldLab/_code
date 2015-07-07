setwd('~/Documents/_Fernald_lab/_dyadStats');
library(WGCNA);
source('../_code/_bootLib/bootstrapFunctions_6-16-13.R');

w2 = read.csv('3setsExplorationData_forR_setw2.csv');##WRONG FILES???!!!
w2[,1]=as.character(w2[,1]);

w2entries = w2[grepl('entered pot', w2[,1]),2:5];
w2entries=cbind(subj=paste('subj',1:8,sep=''),env=c(rep('P',4),rep('NP',4)),w2entries);
names(w2entries)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {w2entries[,bin]=as.numeric(as.character(w2entries[,bin]))}; rm(bin);
w2entriesP = w2entries[1:4,3:6]; rownames(w2entriesP) = w2entries$subj[1:4];
w2entriesNP = w2entries[5:8,3:6]; rownames(w2entriesNP) = w2entries$subj[5:8];

w2lines = w2[grepl('zone lines crossed', w2[,1]),2:5];
w2lines=cbind(subj=paste('subj',1:8,sep=''),env=c(rep('P',4),rep('NP',4)),w2lines);
names(w2lines)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {w2lines[,bin]=as.numeric(as.character(w2lines[,bin]))}; rm(bin);
w2linesP = w2lines[1:4,3:6]; rownames(w2linesP) = w2lines$subj[1:4];
w2linesNP = w2lines[5:8,3:6]; rownames(w2linesNP) = w2lines$subj[5:8];


w1 = read.csv('3setsExplorationData_forR_setw1.csv',header=F);
w1[,1]=as.character(w1[,1]);

w1entries = w1[grepl('entered pot', w1[,1]),2:5];
w1entries=cbind(subj=paste('subj',1:10,sep=''),env=c(rep('P',5),rep('NP',5)),w1entries);
names(w1entries)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {w1entries[,bin]=as.numeric(as.character(w1entries[,bin]))}; rm(bin);
w1entriesP = w1entries[1:5,3:6]; rownames(w1entriesP) = w1entries$subj[1:5];
w1entriesNP = w1entries[6:10,3:6]; rownames(w1entriesNP) = w1entries$subj[6:10];

w1lines = w1[grepl('zone lines crossed', w1[,1]),2:5];
w1lines=cbind(subj=paste('subj',1:10,sep=''),env=c(rep('P',5),rep('NP',5)),w1lines);
names(w1lines)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {w1lines[,bin]=as.numeric(as.character(w1lines[,bin]))}; rm(bin);
w1linesP = w1lines[1:5,3:6]; rownames(w1linesP) = w1lines$subj[1:5];
w1linesNP = w1lines[6:10,3:6]; rownames(w1linesNP) = w1lines$subj[6:10];


s1 = read.csv('3setsExplorationData_forR_setS1.csv',header=F);
s1[,1]=as.character(s1[,1]);

s1entries = s1[grepl('entered pot', s1[,1]),2:5];
s1entries=cbind(subj=paste('subj',1:10,sep=''),env=c(rep('NP',5),rep('P',5)),s1entries);
names(s1entries)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {s1entries[,bin]=as.numeric(as.character(s1entries[,bin]))}; rm(bin);
s1entriesNP = s1entries[1:5,3:6]; rownames(s1entriesNP) = s1entries$subj[1:5];
s1entriesP = s1entries[6:10,3:6]; rownames(s1entriesP) = s1entries$subj[6:10];

s1lines = s1[grepl('zone lines crossed', s1[,1]),2:5];
s1lines=cbind(subj=paste('subj',1:10,sep=''),env=c(rep('NP',5),rep('P',5)),s1lines);
names(s1lines)[3:6]=c('0-5min','5-10min','10-15min','total');
for (bin in 3:6) {s1lines[,bin]=as.numeric(as.character(s1lines[,bin]))}; rm(bin);
s1linesNP = s1lines[1:5,3:6]; rownames(s1linesNP) = s1lines$subj[1:5];
s1linesP = s1lines[6:10,3:6]; rownames(s1linesP) = s1lines$subj[6:10];
#######################################
par(mfrow=c(1,3))
verboseScatterplot(w2entries$total,w2lines$total,abline=T,abline.color='red',type='n',xlab='total pot entries',ylab='total line crosses',main='w2 dyads\n',frame.plot=F,xlim=c(0,15),ylim=c(0,80));
text(w2entries$total,w2lines$total,paste(gsub('subj','',w2entries$subj),'-',w2entries$env,sep=''),cex=1.3);

verboseScatterplot(w1entries$total,w1lines$total,abline=T,abline.color='red',type='n',xlab='total pot entries',ylab='total line crosses',main='w1 dyads\n',frame.plot=F,xlim=c(0,15),ylim=c(0,80));
text(w1entries$total,w1lines$total,paste(gsub('subj','',w1entries$subj),'-',w1entries$env,sep=''),cex=1.3);

verboseScatterplot(s1entries$total,s1lines$total,abline=T,abline.color='red',type='n',xlab='total pot entries',ylab='total line crosses',main='s1 dyads\n',frame.plot=F,xlim=c(0,15),ylim=c(0,80));
text(s1entries$total,s1lines$total,paste(gsub('subj','',s1entries$subj),'-',s1entries$env,sep=''),cex=1.3);

BIN = 1;
out=bootstrap.2independent(w2entriesP[,BIN],w2entriesNP[,BIN], Func='median',
						   dataDescriptor=paste('pot entries,',names(w2entriesP)[BIN]), 
						   groupNames=c('P','NP'), cex=1.3
						   ,trials=100000);
						   
BIN = 1;
out=bootstrap.2independent(w2linesP[,BIN],w2linesNP[,BIN], Func='median',
						   dataDescriptor=paste('line crosses,',names(w2linesP)[BIN]), 
						   groupNames=c(' P','NP'), cex=1.3
						   ,trials=100000);
						   
for (BIN in 1:4) {print(wilcox.test(w2entriesP[,BIN],w2entriesNP[,BIN]))}; rm(BIN);

#######################################

w2lines.df = data.frame(subj=as.factor(rep(w2lines$subj,3)),
						env=as.factor(rep(w2lines$env,3)),
						time=as.factor(c(rep('t0-5',8),rep('t5-10',8),rep('t10-15',8))),
						lines=c(w2lines$'0-5min',w2lines$'5-10min',w2lines$'10-15min')
						);
						
w2entries.df = data.frame(subj=as.factor(rep(w2entries$subj,3)),
						env=as.factor(rep(w2entries$env,3)),
						time=as.factor(c(rep('t0-5',8),rep('t5-10',8),rep('t10-15',8))),
						entries=c(w2entries$'0-5min',w2entries$'5-10min',w2entries$'10-15min')
						);