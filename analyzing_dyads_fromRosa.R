##########################
### load and prep data ###
##########################
# load data
dat = read.csv('~/Documents/_Fernald_lab/fish_lengths_dyads_clean.csv', row.names=1);

# view data
dat;

# remove redundant column
dat = dat[, -7];

# view column names
names(dat);

# clean up some names
names(dat)[c(1,2,5,6)] = c('P.fish.length.cm', 'NP.fish.length.cm', 'Winner', 'Density.Grown');

# add column of length differences
dat = cbind(dat, Length.diff = dat$P.fish.length.cm - dat$NP.fish.length.cm);

# load additional info on GSI, cortisol, aggressive behaviors
dat2 = read.csv('~/Documents/_Fernald_lab/Totaltable dyads.csv', sep='\t');
names(dat2)[7] = '11.KT';
dat2winners = dat2[1:21, ];
rownames(dat2winners) = paste('Dyad ', dat2winners[, 1], sep = '');
dat2winners = dat2winners[, -1];

# combine with existing data
dat = cbind(dat, dat2winners[match(rownames(dat), rownames(dat2winners)), c(4:6, 8:10)]);


#####################################
### define groups and subset data ###
#####################################
# get rows where NP fish won
NPwinner = dat$Winner == 0;

# get rows where fish were age matched
ageMatched = dat$Age.P.fish == dat$Age.NP.fish;

# subset dat to rows where NP fish won
datNP = dat[NPwinner, ];

# subset datNP to rows where fish were age matched
datNPmatch = datNP[ageMatched[NPwinner], ];

# subset datNP to rows where fish were not age matched
datNPnomatch = datNP[!ageMatched[NPwinner], ];

###########################################################################
### investigate length differences using bootstrap2independent function ###
###########################################################################
# read in source code for function, may need to specify directory
source('_code/bootstrap_testMay2013.R');

# run function bootstrap2independent() - only works/is valid with 2 independent datasets
# Arguments are:
# 'group1' and 'group2': data, must be numeric vectors
# 'trials': number of resampling runs
# 'Func': function used to compute representative number for each dataset
#     Func(group1) - Func(group2) = test statistic
# 'replace': boolean indicating whether to resample with/without replacement
# 'plots': boolean indicating whether histograms and boxplot should be produced
# 'col': color of histogram bars and data points in boxplot
# 'border': color of histogram bar borders
# 'col.line': color of mean/median/statistic lines in histograms
# 'dataDescriptor': optional text to label axes in plots
# 'groupNames': vector containing group names
# 'jitter': determines spacing of individual data points in boxplot
# 'pch': determines shape of individual data points in boxplot
# 'boxLineMedian': boolean indicating whether boxplot should show median/mean
# 'printResults': boolean indicating whether to print summary of output to console
# 'verbose': boolean indicating whether to print resampling run number every 1000x

# don't need to explicitly set argument if ok with default
# only 'group1', 'group2', 'trials', 'Func', and 'replace' will actually affect output
# rest of arguments only deal with plots or printing info to console

boot.out = bootstrap2independent(group1 = datNPmatch$Length.diff, 
								 group2 = datNPnomatch$Length.diff,
								 trials = 100000,
								 Func = 'mean',
								 replace = T,
								 plots = T,
								 col = 'grey',
								 border = 'darkgrey',
								 col.line = 'red',
								 dataDescriptor = 'P - NP length difference (cm)',
								 groupNames = c('Age matched', 'NP older'),
								 jitter = .15,
								 pch = 21,
								 boxLineMedian = F,
								 printResults = T,
								 verbose = T
								 );
								 
# store results of multiple tests in list
# first check that dataframes match
check1 = sum(names(datNP) == names(datNPmatch));
check2 = sum(names(datNP) == names(datNPnomatch));
if (check1 == ncol(datNP) & check1 == check2)
{
	cat('\nOK!\n');
} else {
	stop('DATAFRAMES DON\'T MATCH');
}
rm(check1, check2);

# skip 'Winner' column
toTest = c(1:4, 6:12);

boot.outLIST = list();
for (column in toTest)
{
	boot.outLIST[[column]] = bootstrap2independent(group1 = datNPmatch[, column], 
								 group2 = datNPnomatch[, column],
								 trials = 100000,
								 Func = 'mean',
								 replace = T,
								 plots = F,
								 col = 'grey',
								 border = 'darkgrey',
								 col.line = 'red',
								 dataDescriptor = names(datNP)[column],
								 groupNames = c('Age matched', 'NP older'),
								 jitter = .15,
								 pch = 21,
								 boxLineMedian = F,
								 printResults = T,
								 verbose = T
								 );
	names(boot.outLIST)[column] = names(datNP)[column];
}							 
rm(toTest, column);
								 
### OR ###
#######################################################################
### investigate length differences using Wilcoxon and plot manually ###
#######################################################################
# open plotting window with 3 cells in one row
par(mfrow=c(1,3)); 

# plot histograms
hist(datNPmatch$Length.diff, 
	 main = paste('age matched (n=', length(datNPmatch$Length.diff), ')', sep = ''),
	 xlab = 'P length - NP length (cm)', 
	 xlim = c(-.3, .3),
	 ylim = c(0, 3),
	 breaks = length(datNPmatch$Length.diff),
	 col = 'grey',
	 border = 'darkgrey'
	 ); 
abline(v = mean(datNPmatch$Length.diff));
text(0.1, 2, 'solid line = mean');
abline(v = median(datNPmatch$Length.diff), lty = 'dashed');
text(0.1, 1.9, 'dashed line = median');
hist(datNPnomatch$Length.diff, 
	 main = paste('NP fish older (n=', length(datNPnomatch$Length.diff), ')', sep = ''),
	 xlab = 'P length - NP length (cm)', 
	 xlim = c(-.3, .3),
	 ylim = c(0, 3),
	 breaks = length(datNPnomatch$Length.diff),
	 col = 'grey',
	 border = 'darkgrey'
	 );
abline(v = mean(datNPnomatch$Length.diff));
text(0.15, 2, 'solid line = mean');
abline(v = median(datNPnomatch$Length.diff), lty = 'dashed');
text(0.15, 1.9, 'dashed line = median');

# compute wilcoxon p-value
wilcoxP = wilcox.test(datNPmatch$Length.diff, datNPnomatch$Length.diff)$p.value;

# boxplot
toPlot = c(datNPmatch$Length.diff, datNPnomatch$Length.diff);
grp = c(rep('age matched', length(datNPmatch$Length.diff)), 
		rep('NP fish older', length(datNPnomatch$Length.diff))
		);
boxplot(toPlot ~ grp, 
	    ylab = 'P length - NP length (cm)',
	    main = paste('p = ', signif(wilcoxP, 2), ' (wilcoxon rank sum test)', sep = ''),
	    cex.axis = 1.2,
	    cex.lab = 1.2
	    );
stripchart(toPlot ~ grp, 
		   vertical = T,
		   add = T, 
		   method = 'jitter',
		   jitter = .15,
		   pch = 21,
		   bg = 'grey',
		   cex = 1.5
		   );
		   
############
### baseline aggression data
############

bagg = read.csv('baselineAggression.csv',row.names=1);
baggP = bagg[,1:6];