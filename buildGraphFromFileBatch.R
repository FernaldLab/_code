buildGraphFromFileBatch = function(dir = 'folder of log files', 
								   writeOutProbMats = T,
								   sigFigs = 6,
								   groupChar = 4
								   )
{
	
	# check working directory and get filenames
	wd = strsplit(getwd(), split = '/', fixed = T)[[1]];
	if (wd[length(wd)] == dir)
	{
		log_files = list.files();
	}
	else
	{
		stop('CHANGE WORKING DIRECTORY TO FOLDER CONTAINING LOGS...USE setwd()');
	}
	
	# get group for each file and group names
	grouping = substr(log_files, groupChar, groupChar);
	groupNames = names(table(grouping));
	
	# lists to hold probability matrices for each group
	gp1mats = list();
	gp2mats = list();
	
	# loop through logs
	for (log in 1:length(log_files))	
	{
		# get filename
		this_log = log_files[log];
		
		# check that log is .txt file
		this_ext = substring(this_log, first = (nchar(this_log)-3));
		if (this_ext != '.txt')
		{
			stop('WRONG FILE TYPE DETECTED: ', 
				  gsub(paste(dir, '/', sep = ''), '', this_log), 
				  '\n...LOGS CAN ONLY BE .txt FILES!'
				  );
		}
		
		# build .dot file and compute probability matrix
		temp = buildGraphFromFile(this_log, dotFileBase=NULL);
		
		# get group
		this_group = grouping[log];
		
		# store probability matrix in group list
		if (this_group == groupNames[1])
		{
			gp1mats[[log]] = temp$matrix;
			names(gp1mats)[log] = log_files[log];
		}
		else
		{
			gp2mats[[log]] = temp$matrix;
			names(gp2mats)[log] = log_files[log];
		}
		
		# write .csv file of probability matrix
		if (writeOutProbMats)
		{
			write.csv(signif(temp$matrix, sigFigs), 
				      file = gsub('txt', 'csv', this_log)
				      );
		}
		
	}
	
	probMats = list(gp1mats, gp2mats);
	names(probMats) = groupNames;
	
	OUT = list(files = log_files,
			   groups = grouping,
			   probMats = probMats
			   );
	
	return(OUT);
}