rm(list=ls());
setwd('~/Documents/_Fernald_lab/_dyadStats');


db.agg = read.csv('DominanceAggression-AggressionDmales_sheet.csv');
names(db.agg)[c(8, 11, 13, 14, 16, 18, 19)] = c('Fish.ID', '1st.kill.date', 'Age.at.1st.kill', '2nd.kill.date', '3rd.kill.date', 'num.females.killed', 'num.males.killed');

db.testes = read.csv('DominanceAggression-extractDBtestes_sheet.csv');
names(db.testes)[c(4, 5, 15)] = c('Weight_fish.gr', 'Length_Fish.mm', 'tissue_weight.gr');
db.testes = db.testes[match(db.agg$Fish.ID, db.testes$FISH_ID), ];

# checks
dim(db.agg)[[1]] == dim(db.testes)[[1]];
sum(db.agg$Fish.ID == db.testes$FISH_ID) == dim(db.agg)[[1]];
sum(as.character(db.agg$Dissected) == as.character(db.testes$DISSECTION_DATE)) == dim(db.agg)[[1]];
diff = which(as.character(db.agg$Dissected) != as.character(db.testes$DISSECTION_DATE));
cbind(as.character(db.agg$Dissected)[diff], as.character(db.testes$DISSECTION_DATE)[diff]);
rm(diff);

# some discrepancies in dissection dates, use ones from db.testes
db = cbind(db.agg[, 1:5], db.testes[, c(1, 4, 5, 15, 18)], db.agg[, 7:20]);
names(db) = toupper(names(db));

###
db = cbind(db, DENSITY.HIGH = db$RUN >= 12);

densityMat = matrix(nrow = 2, 
					ncol = 5,		# ncol = 4 if no FKT
					dimnames = list(c('high', 'low'), 
									names(table(db$STATUS))
									)
					);	
# DIFFERENT RESULTS? BUILDING densityMat WRONG?	
for (s in names(table(db$STATUS)))
{
	densityMat[1, colnames(densityMat)==s] = sum(db$DENSITY.HIGH[db$STATUS==s]);
	densityMat[2, colnames(densityMat)==s] = sum(!db$DENSITY.HIGH[db$STATUS==s]);
}
rm(s);

# statuses = names(table(db$STATUS));
# status.lookup = data.frame(statuses, c(2, 5, 3, 4, 1)); #rm(statuses);
# statusvec = c();
# for (fish in 1:nrow(db))
# {
	# statusvec[fish] = status.lookup[status.lookup[, 1] == as.character(db$STATUS[fish]), 2];
# }
# rm(fish);

# db = cbind(db, statusvec);
# names(db)[25] = 'STATUS.CODED';

####################
library(WGCNA); allowWGCNAThreads(); options(stringsAsFactors=F);

toTest = 'AGE.AT.DYAD.CREATION';
toTestCol = match(toTest, names(db));

verboseBoxplot(db[, toTestCol], db$STATUS, xlab = 'Status', ylab = toTest);


par(mfrow=c(2, 5));
for (s1 in 1:length(statuses))
{
	temp1 = db[as.character(db$STATUS) == statuses[s1], toTestCol];
	for (s2 in (s1+1):length(statuses))
	{
		temp2 = db[as.character(db$STATUS) == statuses[s2], toTestCol];
		verboseBoxplot(c(temp1, temp2), 
					   c(rep(statuses[s1], length(temp1)), rep(statuses[s2], length(temp2))), 
					   xlab = 'Status', 
					   ylab = toTest,
					   ylim = c(min(db[, toTestCol]), max(db[, toTestCol]))
					   );
	}
}
rm(toTest, toTestCol, s1, s2, temp1, temp2);


######################
par(mfrow=c(2, 4));
cex.lab=1.2;
ylim=c(0, 60);
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[!db$STATUS=='T'], 
			   db$STATUS[!db$STATUS=='T'], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('AT', 'FKT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'FKT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('AT', 'KT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'KT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('AT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'PKT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('FKT', 'KT')], 
			   db$STATUS[db$STATUS %in% c('FKT', 'KT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('FKT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('FKT', 'PKT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(db$DAYS.IN.DYAD.UNTIL.1ST.KILL[db$STATUS %in% c('KT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('KT', 'PKT')], 
			   ylab='DAYS.IN.DYAD.UNTIL.1ST.KILL', 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
			   
##################
par(mfrow=c(2, 4));
cex.lab=1.2;
ylim=c(100, 250);
ylab='AGE.AT.1ST.KILL';
dat=db$AGE.AT.1ST.KILL;
verboseBoxplot(dat[!db$STATUS=='T'], 
			   db$STATUS[!db$STATUS=='T'], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('AT', 'FKT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'FKT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('AT', 'KT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'KT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('AT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('AT', 'PKT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('FKT', 'KT')], 
			   db$STATUS[db$STATUS %in% c('FKT', 'KT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('FKT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('FKT', 'PKT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );
verboseBoxplot(dat[db$STATUS %in% c('KT', 'PKT')], 
			   db$STATUS[db$STATUS %in% c('KT', 'PKT')], 
			   ylab=ylab, 
			   xlab='', 
			   notch=F,
			   cex.lab=cex.lab,
			   ylim=ylim
			   );