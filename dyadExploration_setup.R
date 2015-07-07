

df.z=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t0-5',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,1]
		   );
		   
tmp=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t5-10',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,2]
		   );
df.z=rbind(df.z,tmp);

tmp=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t10-15',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,3]
		   );
df.z=rbind(df.z,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t0-5',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,4]
		   );
df.z=rbind(df.z,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t5-10',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,5]
		   );
df.z=rbind(df.z,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t10-15',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   lines=z[,6]
		   );
df.z=rbind(df.z,tmp);
##########
#########
########
######
####
##
#

df.entries=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t0-5',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=NPentries[,1]
		   );
		   
tmp=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t5-10',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=NPentries[,2]
		   );
df.entries=rbind(df.entries,tmp);

tmp=data.frame(env=as.factor(rep('NP',14)), 
		   time=as.factor(rep('t10-15',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=NPentries[,3]
		   );
df.entries=rbind(df.entries,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t0-5',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=Pentries[,1]
		   );
df.entries=rbind(df.entries,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t5-10',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=Pentries[,2]
		   );
df.entries=rbind(df.entries,tmp);

tmp=data.frame(env=as.factor(rep('P',14)), 
		   time=as.factor(rep('t10-15',14)), 
	       subj=as.factor(paste('subj',1:14,sep='')), 
		   entries=Pentries[,3]
		   );
df.entries=rbind(df.entries,tmp);

######

par(mfrow=c(2,3));
boxplot(NPentries[,1:3],ylim=c(0,9),main='NP',ylab='entries');
boxplot(Pentries[,1:3],ylim=c(0,9),main='P',ylab='entries');
boxplot(rbind(NPentries[,1:3],Pentries[,1:3]),ylim=c(0,9),main='all',ylab='entries');

boxplot(NPlines[,1:3],ylim=c(0,40),main='NP',ylab='lines');
boxplot(Plines[,1:3],ylim=c(0,40),main='P',ylab='lines');
boxplot(rbind(NPlines[,1:3],Plines[,1:3]),ylim=c(0,40),main='all',ylab='lines');
