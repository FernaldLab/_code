par(mfrow=c(1,2))

verboseScatterplot(datNP2$Age.NP.fish, datNP2$Cortisol, abline=T, xlab='Age', ylab='Cortisol',main='NP winners (n=15)\n',type='n', frame.plot=T, abline.col='red');
text(datNP2$Age.NP.fish, datNP2$Cortisol, gsub('Dyad ','',rownames(datNP2)));



verboseScatterplot((datNP2$Age.P.fish - datNP2$Age.NP.fish), datNP2$Length.diff, ylab='Length diff (P - NP)', xlab='Age diff (P - NP)', abline=T, main='NP winners (n=15)\n', type='n', frame.plot=T, abline.col='red', xlim=c(-20, 0));
text((datNP2$Age.P.fish - datNP2$Age.NP.fish), datNP2$Length.diff, gsub('Dyad ','',rownames(datNP2)))

######################

win.ratio.length = c();
win.ratio.age = c();

for(row in 1:nrow(dat))
{
	if (dat$Winner[row] == 1)
	{
		win.ratio.length = c(win.ratio.length, dat$P.fish.length.cm[row]/dat$NP.fish.length.cm[row]);
		win.ratio.age = c(win.ratio.age, dat$Age.P.fish[row]/dat$Age.NP.fish[row]);
	}
	else if (dat$Winner[row] == 0)
	{
		win.ratio.length = c(win.ratio.length, dat$NP.fish.length.cm[row]/dat$P.fish.length.cm[row]);
		win.ratio.age = c(win.ratio.age, dat$Age.NP.fish[row]/dat$Age.P.fish[row]);
	}
}
rm(row);

par(mfrow=c(1, 2));
verboseScatterplot(win.ratio.length,win.ratio.age,abline=T, xlab = 'Length D/ND', ylab = 'Age D/ND', abline.col='red', type = 'n', main='All dyads\n');
text(win.ratio.length,win.ratio.age, gsub('Dyad ','',rownames(dat)));

verboseScatterplot(win.ratio.length[dat$Winner==0],win.ratio.age[dat$Winner==0],abline=T, xlab = 'Length D/ND', ylab = 'Age D/ND', abline.col='red', type = 'n', main='NP winners\n');
text(win.ratio.length[dat$Winner==0],win.ratio.age[dat$Winner==0], gsub('Dyad ','',rownames(dat)));
#or
par(mfrow=c(1, 2));
verboseScatterplot(win.ratio.length,win.ratio.age,abline=T, xlab = 'Length D/ND', ylab = 'Age D/ND', abline.col='red', main='All dyads (n=21)\n',frame.plot=F,corOptions="method='p'",pch=21,col='darkgrey',bg='grey',cex=1.5);
verboseScatterplot(win.ratio.length[dat$Winner==0],win.ratio.age[dat$Winner==0],abline=T, xlab = 'Length D/ND', ylab = 'Age D/ND', abline.col='red', main='NS winners (n=15)\n',frame.plot=F,corOptions="method='p'",pch=21,col='darkgrey',bg='grey',cex=1.5);

# verboseScatterplot(win.ratio.length[dat$Winner==1],win.ratio.age[dat$Winner==1],abline=T, xlab = 'Length D/ND', ylab = 'Age D/ND', abline.col='red', type = 'n', main='P winners\n');
# text(win.ratio.length[dat$Winner==1],win.ratio.age[dat$Winner==1], gsub('Dyad ','',rownames(dat)));

################

winner.age=c()
for (r in 1:nrow(dat))
{
	if (dat$Winner[r]==1)
	{
		winner.age=c(winner.age, dat$Age.P.fish[r]);
	}
	else if (dat$Winner[r]==0)
	{
		winner.age=c(winner.age, dat$Age.NP.fish[r]);
	}
}
rm(r);

################

verboseBoxplot(win.ratio.length,as.factor(win.ratio.age),notch=F,ylab='Length (D/ND)',xlab='Age difference (days)',names=c('0 (n=11)','6 (n=5)','19 (n=5)'),border = 'black',main='All winners\nkruskal-wallis, ');
stripchart(as.data.frame(cbind(dat[,3],dat[,2],dat[,1])), vertical=T, add=T,pch=21,col='black',bg='grey',cex=1.5);

difflengthmat=matrix(nrow=11,ncol=3);
difflengthmat[1:5,1]=win.diff.length[1:5];
difflengthmat[1:5,2]=win.diff.length[6:10];
difflengthmat[1:11,3]=win.diff.length[11:21];

par(oma=c(0,2,0,0));
verboseBoxplot(win.diff.length,as.factor(win.ratio.age),notch=F,ylab='D - ND Length diff (cm)',xlab='Age difference (days)',names=c('0 (n=11)','6 (n=5)','19 (n=5)'),border = 'black',main='All winners, ');
stripchart(as.data.frame(cbind(difflengthmat[,3], difflengthmat[,2], difflengthmat[,1])), vertical=T, method='jitter',add=T,pch=21,col='black',bg='grey',cex=1.5,jitter=.2);
  

################
win.diff.length=c();
for (row in 1:nrow(alldat)) {
	if (alldat$Winner[row]==0) {
		win.diff.length=c(win.diff.length, alldat$NP.fish.length.cm[row]-alldat$P.fish.length.cm[row]);
	} else {
		win.diff.length=c(win.diff.length, alldat$P.fish.length.cm[row]-alldat$NP.fish.length.cm[row]);
	}
	
}; rm(row)

par(oma=c(0,2,0,0))
verboseScatterplot(win.diff.age, win.diff.length, 
				   abline=T, abline.col='red',
				   pch=21, bg='grey', col='black',
				   frame.plot=F, 
				   xlab='Winner - loser age (days)', ylab='Winner - loser length (cm)',
				   cex=2, 
				   xlim=c(-5, 20), ylim=c(-.3,.3), 
				   main='All dyads (n=21)\n'
				   );
abline(h=0,lty='dashed',col='darkgrey');
abline(v=0,lty='dashed',col='darkgrey');
text(1.2, .2, 'n=6');
text(1.2, .095, 'n=3');
text(6, .03, 'n=2');
text(19, -.02, 'n=2');

###################
# new figure 3

#par(mfrow=c(1,2),oma=c(0,1,2,0));
par(mfrow=c(1,3));
cex.lab=1.3;
cex.axis=1.3;
data1=datNPnomatch[,11:12];
data2=datNPmatch[,11:12]

boxplot(data1, ylab = 'Aggressive behaviors / min', xlab='',
				names=c('1st day','2nd day'),
				main = 'NS older, p < 1e-5 (n=10)', frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,ylim=c(0,7)
				);
				
# add data points
stripchart(data1, vertical = T, add = T, pch =21, bg ='grey', cex = 1.5);
for (row in 1:nrow(data1)) {segments(1, data1[row, 1], 2, data1[row, 2], col = 'darkgrey')}


boxplot(data2, ylab = 'Aggressive behaviors / min',xlab='',
				names=c('1st day','2nd day'),
				main = 'Age matched, p = 0.66 (n=5)', frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,ylim=c(0,7)
				);
				
# add data points
stripchart(data2, vertical = T, add = T, pch =21, bg ='grey', cex = 1.5);
for (row in 1:nrow(data2)) {segments(1, data2[row, 1], 2, data2[row, 2], col = 'darkgrey')}

#title('NS winners', outer=T)

data1=data1[,1]-data1[,2];
data2=data2[,1]-data2[,2];
boxplot(-cbind(data1,c(data2,rep(NA,5))), ylab = 'Aggressive behaviors / min, 2nd - 1st day', xlab='',
				names=c('NS older (n=10)','Age matched (n=5)'),
				main = 'p = 0.023', frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis#,ylim=c(0,7)
				);
				
# add data points
stripchart(list(-data1,-data2), vertical = T, add = T, pch =21, bg ='grey', cex = 1.5);

#############


par(mfrow=c(1,2));
cex.lab=1.3;
cex.axis=1.3;
#data1=datNPnomatch[,11:12];
#data2=datNPmatch[,11:12];
data1=datNotMatched[,10:11];
data2=datMatched[,10:11];

boxplot(data1, ylab = 'Aggressive behaviors / min', xlab='',
				names=c('1st day','2nd day'),
				main = 'Older, p < 1e-5 (n=10)', frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,ylim=c(0,7)
				);
				
# add data points
stripchart(data1, vertical = T, add = T, pch =21, bg ='grey', cex = 1.5);
for (row in 1:nrow(data1)) {segments(1, data1[row, 1], 2, data1[row, 2], col = 'darkgrey')}


boxplot(data2, ylab = 'Aggressive behaviors / min',xlab='',
				names=c('1st day','2nd day'),
				main = 'Age matched, p = 0.23 (n=11)', frame.plot = F,
				cex.lab = cex.lab,
				cex.axis = cex.axis,ylim=c(0,7)
				);
				
# add data points
stripchart(data2, vertical = T, add = T, pch =21, bg ='grey', cex = 1.5);
for (row in 1:nrow(data2)) {segments(1, data2[row, 1], 2, data2[row, 2], col = 'darkgrey')}