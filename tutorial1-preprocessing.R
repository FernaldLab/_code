##### Microarray pre-processing tutorial
##### Austin Hilliard, White lab UCLA, August 2011

####################
##### sections #####
####################
#####	1) start up
#####	2) load raw expression data
#####	3) remove outlier probe measurements 
#####		- try 3 different outlier thresholds
#####		- every subsequent step will be done on results from each of the 3 thresholds
#####	4) remove probes with too many outlier measurements removed
#####	5) remove outlier samples
#####	6) quantile normalization
#####	7) remove unannotated probes
#####	8) restrict to one representative probe for each gene

##### at the end of most sections is a call to the save() function
##### this is so you can load the data at that particular stage of pre-processing later using load()
##### each filename I used is in a comment but you should use whatever name is most sensible to you

####################
##### start up #####
####################

# clear workspace
rm( list = ls() );

# enter the name of your working direstory inside the quotes
setwd('~/Documents/_analysis_compare');

# load WGCNA functions and set options
library('WGCNA');
options(stringsAsFactors = F);

# load preprocessing functions
# enter the name of directory and file within quotes
# if you have this file in your working directory you can just enter the filename
# e.g. I use > source('_code/preProcATH-July2011.R')
source('');

####################################
##### load raw expression data #####
####################################

# 'S_no21.csv' is the raw expression data from VSP
# 'X_no21.csv' is the raw expression data from area X
# rows are probes (n=42921), columns are samples (n=26)
# we'll work on VSP so I've commented out lines below pertaining to area X data to avoid mistakes

datVSP0 = read.csv('_expressionData/S_no21.csv', row.names = 1);
#datX0 = read.csv('_expressionData/X_no21.csv', row.names = 1);

# strip leading 'X' from sample names
names(datVSP0) = substring( names(datVSP0), 2);
#names(datX0) = substring( names(datX0), 2);

# check dimensions and view part of the data
# should by 42921 rows by 26 columns
dim(datVSP0); head(datVSP0);
#dim(datX0); head(datX0);

# store whichever dataset you plan to work on in the variable DATA
# instead of editing everything below each time you want to use different data you just have to change this one line

DATA = datVSP0;
## OR ##
#DATA = datX0;

#################################
##### remove outlier probes #####
#################################

# try different values for stdevs from mean a probe has to be to be an outlier
# i.e. removeOutlierProbesIterate() parameter 'deviate' = 2, 2.5, and 3
deviate_vals = c(2, 2.5, 3);

# list to hold outputs using different values of deviate parameter
outProbes = list();

# loop through deviate_vals, find probe outlier measurements for each value
# store results in outProbes, which will be 3 component list after loop
# each component is output list from removeOutlierProbesIterate() using certain value for deviate
# thus outProbes is a list where each component is itself a list with multiple components

for( j in 1:3 ) {
	outProbes[[j]] = removeOutlierProbesIterate(DATA, deviate = deviate_vals[j], rowORcol = 1);
	}
	
	# always clean up variables defined during loop
	rm(j);	
	
# name components of outProbes based on which value of deviate they were computed using
# paste() is a very useful function
names(outProbes) = paste('deviate', deviate_vals, sep = '');

# clean up memory, usually good idea after working on multiple large datasets
collectGarbage();

# save for later use, I named the file 'VSPprobes_removed_deviate-2-2.5-3'
save(outProbes, file = '');

###################################################
##### remove probes/samples with too many NAs #####
###################################################

# list to hold outputs, analagous to outProbes
# will be 3 component list, each holding results based on corresponding outProbes component
outNAs = list();

# loop through outProbes 
# check the cleaned data from each run of removeOutlierProbesIterate for too many NAs

for( j in 1:3 ) {
	
	# extract cleaned expression data from outProbes
	dat = outProbes[[j]]$dataClean;
	
	# print deviate value used during outlier probe removal
	# cat() allows printing of strings and variable values together
	# '\n' is the newline character, i.e. same as hitting return
	cat('\ndeviate = ', outProbes[[j]]$sd_cutoff, sep = '');
	
	# set thresholds for how many NAs are too many 
	# floor() rounds input down to nearest integer
	# could try different thresholds to see how network construction would be affected
	probe_thresh = floor( ncol(dat) / 2 );
	sample_thresh = floor( nrow(dat) / 2 );
	outNAs[[j]] = removeTooManyNAs(dat, probe_thresh = probe_thresh, sample_thresh = sample_thresh);
	
	}
	
	# always clean up variables defined during loop
	rm(j, dat, probe_thresh, sample_thresh);
	
# name components of outNAs to match outProbes
names(outNAs) = names(outProbes);

# clean up memory, usually good idea after working on multiple large datasets
collectGarbage();

# save for later use, I named the file 'VSPprobes_removed_deviate-2-2.5-3_tooManyNAs'
save(outNAs, file = '');

##################################
##### remove outlier samples #####
##################################

# list to hold outputs, analagous to outProbes and outNAs
# will be 3 component list, each holding results based on corresponding outNAs component
outSamples = list();

# don't loop because of issue with outlierSamplesIterate2 needing user input on sample removal
# remove outlier samples for each cleaned dataset (outlier probes removed with deviate = 2, 2.5, or 3, and probes with too many NAs removed)
# store in corresponding component of outSamples 

# run these next 3 lines ONE AT A TIME
outSamples[[1]] = outlierSamplesIterate2( outNAs[[1]]$dataClean ); # I only removed 4S
outSamples[[2]] = outlierSamplesIterate2( outNAs[[2]]$dataClean ); # I removed 3S, 1S 14S
outSamples[[3]] = outlierSamplesIterate2( outNAs[[3]]$dataClean ); # I removed 1S, 3S

# name components of outSamples to match outProbes and outNAs
names(outSamples) = names(outProbes);

# clean up memory, usually good idea after working on multiple large datasets
collectGarbage();

# save for later use, I named the file 'VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed'
save(outSamples, file = '');

##################################
##### quantile normalization #####
##################################

# list to hold results of quantile normalization on each cleaned dataset (outlier probes removed with deviate = 2, 2.5, or 3, probes with too many NAs removed, outlier samples removed)
# analogous to previous output lists
outQnorm = list();

# loop through outSamples

for( j in 1:3 ) {
	
	# extract cleaned expression data from outSamples
	dat = outSamples[[j]]$dataClean;
	
	# normalize.quantiles() needs matrix as input, and returns a matrix as output
	# convert dat to matrix to normalize, then convert back to data frame
	# store in corresponding component of outQnorm
	outQnorm[[j]] = as.data.frame( normalize.quantiles( as.matrix(dat) ) );
	
	# restore probe and sample names since they were lost in conversions
	# dimnames() returns a list with 2 components: rownames and colnames
	dimnames( outQnorm[[j]] ) = dimnames(dat);
	
	}
	
	# always clean up variables defined during loop
	rm(j, dat);
	
# name components of outQnorm to match previous output lists	
names(outQnorm) = names(outSamples);

# check dimensions of output
lapply(outQnorm, dim);

# clean up memory, usually good idea after working on multiple large datasets
collectGarbage();

# save for later use, I named the file 'VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm'
save(outQnorm, file = '');

########################################################
##### remove probes not annotated with gene symbol #####
########################################################

# reminder: we tried 3 different thresholds for removing outlier probes
#			then, each of 3 resulting datasets had probes with too many NAs removed
#			then, each dataset had outlier samples removed, and was quantile normalized
#
#			now, we'll filter the normalized data 
#			first, we'll throw out probes that aren't annotated with a gene symbol
#			then, for each gene on the array we'll choose 1 representative probe

# load the file 'annotationInfo' that I sent you
# loads data frame called 'annos' that holds annotation info
# if 'annotationInfo' is in your working directory just use > load('annotationInfo')
# otherwise you need to enter the full path
load('');

# restrict to probes annotated with gene symbol, store in new data frame
annos_noNA = annos[ !( is.na(annos$geneSymbol) ) , ];

# list to hold output from each dataset
DATAlist = list();

# loop through outQnorm
# restrict each normalized dataset to annotated probes and store in corresponding component of DATAlist

for( j in 1:3 ) {
	
	# extract cleaned and normalized dataset from outQnorm
	dat = outQnorm[[j]];
	
	# restrict dat to annotated probes
	# match() is very useful function, read the help file for it!
	dat = dat[ match( annos_noNA$probeID, rownames(dat) ), ];
	
	# check if NAs were introduced by matching
	# this will probably happen if probes were removed by removeTooManyNAs after outlier probe removal
	# there are a couple concepts and some syntax that we haven't discussed yet in the following code block
	# I won't try to explain here in comments
	if ( sum( grepl( 'NA', rownames(dat) ) ) > 0 ) {
		toRemove = ( 1:nrow(dat) )[ grepl( 'NA', rownames(dat) ) ];
		cat('NA rows:', toRemove, '\n', sep=' ');
		DATAlist[[j]] = dat[-toRemove, ];
		} else {
			print('no NAs introduced');
			DATAlist[[j]] = dat;
			}
			
	}
	
	# always clean up variables defined during loop
	rm(j, dat, toRemove);
	
# name components of DATAlist to match previous output lists
names(DATAlist) = names(outQnorm);

# check dimensions of output
lapply(DATAlist, dim);

#########################################################
##### restrict to one probe for each annotated gene #####
#########################################################

# remove duplicate probes using collapseRows()
# collapseRows() chooses one representative probe for each gene using various metrics
# we'll stick with default settings

# list to hold outputs
DATAcollapsedList = list();

# loop through DATAlist

for( j in 1:3 ) {
	
	# report which dataset is currently being processed
	# each dataset will take a few mintues to run
	cat('working on ', names(DATAlist)[j], '\n', sep = '');
	
	# extract filtered data from DATAlist
	dat = DATAlist[[j]];
	
	# store all genes represented by probes in filtered dataset
	# %in% is similar to match() but returns T/F instead of numbers
	genes = annos_noNA$geneSymbol[ annos_noNA$probeID %in% rownames(dat) ];
	
	# store all probe ID #'s
	probeIDs = annos_noNA$probeID[ annos_noNA$probeID %in% rownames(dat) ];
	
	# store output of collapseRows() in corresponding component of DATAcollapsedList
	DATAcollapsedList[[j]] = collapseRows(dat, rowGroup = genes, rowID = probeIDs);
	
	}
	
	# always clean up variables defined during loop
	rm(j, dat, genes, probeIDs);
	
# name components of DATAcollapsedList to match previous output lists
names(DATAcollapsedList) = names(DATAlist);

# check dimensions and head of output
lapply(DATAcollapsedList, lapply, dim);
lapply(DATAcollapsedList, lapply, head);

# save for later use, I named the file 'VSPprobes_removed_deviate-2-2.5-3_tooManyNAs_samples_removed_Qnorm_collapsed'
save(DATAcollapsedList, file = '');

# clean up
collectGarbage();