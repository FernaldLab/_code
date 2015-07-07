checkGeneListEnrichment = function(dataset1, dataset2, refset, alt='two.sided') {
	
	data1name = as.character(match.call()[2]);
	data2name = as.character(match.call()[3]);

	incommon = sum(dataset1 %in% dataset2);
	outof = length(dataset2);
	innet = sum(dataset1 %in% refset);
	
	countdimnames = list(c("Y", "N", ""), c("Y", "N", ""));
	names(countdimnames)[[1]] = data2name;
	names(countdimnames)[[2]] = data1name;
	counts0 = matrix(ncol = 3, nrow = 3, dimnames = countdimnames);
	
	counts0[1,] = c(incommon, outof-incommon, outof);
	counts0[3,] = c(innet, length(refset)-innet, length(refset));
	counts0[2,] = counts0[3,] - counts0[1,];
	counts = counts0[-3,-3];
	
	#if (any(counts==0)) {
	#	return('no overlap')
	#}
	
	test = fisher.test(counts, alternative = alt);
	out = list(counts0, test);
	
	return(out);
}
	
#####
#####

checkGeneListEnrichmentList = function(dataset1, 
									   dataset2list, 
									   refset, 
									   only_sig = 0, 
									   thresh = .05, 
									   alt = 'two.sided', 
									   order = T
									   ) {
									   	
	data1name = as.character(match.call()[2]);
	out=list();
	sigs=c(); 
	pvals = vector(mode = "numeric");  
	modpvals = matrix(nrow = length(dataset2list), ncol = 2);
	
	#compute enrichments, store count tables and results of fisher test in list
	for (set in 1:length(dataset2list)) {
		temp = checkGeneListEnrichment(dataset1, dataset2list[[set]], refset, alt = alt);
		
		#print(temp)
		#print(names(dataset2list[set])); print(data1name)
		
		names(dimnames(temp[[1]])) = c(names(dataset2list[set]), data1name);
		out[[set]] = temp;
		names(out)[set] = names(dataset2list[set]);
	}
		
	#store pvals only in data frame for printing to screen and output
	#user chooses to include all pvals, or only pvals below significance threshold 
	module = names(out);
	for (mod in module) {
		if (only_sig == 0) {
			pvals=c(pvals, format(out[[mod]][[2]]$p.value, digits = 4));
		}
		if (only_sig == 1) {
			sigs = c(sigs, out[[mod]][[2]]$p.value < thresh);
			pvals = c(pvals,
				format(out[[mod]][[2]]$p.value[out[[mod]][[2]]$p.value<thresh], digits = 4));
		}
	}
	
	#build pvals data frame	
	if (only_sig == 0){modpvals = cbind(module, pvals)};
	if (only_sig == 1){modpvals = cbind(module[sigs], pvals)};
	modpvals = as.data.frame(modpvals);
	names(modpvals) = c("module", "pval");
	#print(modpvals);
	#modpvals=cbind(modpvals,qvalue(as.numeric(modpvals$pval))$qvalue);
	if(order){modpvals = modpvals[order(as.numeric(modpvals$pval)), ]};
		#names(modpvals)[3]='qval'};
	
	#sorted_pvals=sort(as.numeric(modpvals[,2]));upto=length(sorted_pvals);
#	corrected_pvals=c();
#	corrected_pvals=for(p in 1:upto){corrected_pvals=c(corrected_pvals,p/upto*thresh)};
#	corrected_alpha=corrected_pvals[sorted_pvals<corrected_pvals];
	
	return(list(details = out, pvals = modpvals));
	
}