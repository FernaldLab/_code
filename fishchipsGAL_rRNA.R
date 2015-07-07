# library(limma);
# gal=readGAL('../Fishchip4.03_annotatedGAL_20090730HEM.gal.gal');
# rgal=gal[grepl('ribosomal', ignore.case=T, gal$EST_DFCI_HIT_2009),];
# rgal28=rgal[grepl('28S ribosomal', ignore.case=T, rgal$EST_DFCI_HIT_2009),];
# rgal18=rgal[grepl('18S ribosomal', ignore.case=T, rgal$EST_DFCI_HIT_2009),];
# rgal60=rgal[grepl('60S ribosomal', ignore.case=T, rgal$EST_DFCI_HIT_2009),];
# rgal40=rgal[grepl('40S ribosomal', ignore.case=T, rgal$EST_DFCI_HIT_2009),];

# DAT = rgal28;


# # viewXXS = function(DAT)
# {
	# XXS = c();
	# for (row in 1:nrow(DAT)) {
		# thisrow = strsplit(DAT[row, 9], ' ')[[1]];
		# thisrow = thisrow[grep('[0-9]{1,2}S', thisrow)];
		# if (length(thisrow) > 0 && nchar(thisrow) == 3) {XXS = c(XXS, thisrow)}
	# }
	# return(XXS);
# }



checkHomology = function(DAT)
{
	df = data.frame(matrix(nrow=nrow(DAT), ncol=ncol(DAT)+1));
	
	for (row in 1:nrow(DAT)) {
		#print(row)
		row_split = strsplit(DAT[row, 9], ' ')[[1]];
		check_complete = grep('complete', row_split);#print(check_complete)
		check_partial = grep('partial', row_split);#print(check_partial)
		
		if (length(check_complete) > 0) {
			df[row, ] = c(DAT[row, ], 'complete');
		}
		else if (length(check_partial) > 0) {
			follow_partial = row_split[check_partial + 1];#print(follow_partial)
			check_percent = grep('\\([0-9]{1,2}%\\)', follow_partial);#print(check_percent)
			if (length(check_percent) == 1) {
				toPaste = paste('partial ', follow_partial[check_percent], sep = '');
				df[row, ] = c(DAT[row, ], toPaste);
			}
			else {
				cat('SOMETHING WRONG AT ROW ', row, '\n', 
				    '  row_split=', row_split, '\n',
				    '  check_partial=', check_partial, '\n',
				    '  follow_partial=', follow_partial, '\n',
				    '  check_percent=', check_percent, '\n', 
				    sep = '');
			}
		}
		else {
			cat('SOMETHING WRONG AT ROW ', row, '\n  NO PARTIAL OR COMPLETE?', sep = '');
		}
	}
	colnames(df) = c(colnames(DAT), 'homology');
	completes = df[df$homology=='complete', ];
	partials = df[df$homology!='complete', ];
	
	forRank = c();
	for (row in 1:nrow(partials)) {
		percent = partials$homology[row];
		percent = gsub('partial (', '', percent, fixed=T);
		percent = gsub('%)', '', percent, fixed=T);#print(percent)
		forRank = c(forRank, as.numeric(percent))
	}
	partials = partials[order(forRank, decreasing=T),];
	
	out = list(completes=completes, partials=partials, all=df);
	return(out);	
}
