addGeneSymToDAVIDChart = function(chart, IDs) {
	if (ncol(chart) != 13) {
		stop('chart must be raw output from DAVID Functional Annotation');
	}
	chart_new = cbind(chart, Genes=rep(NA, nrow(chart)));
	names(chart_new)[14] = 'GeneSym';
	termsIDs = strsplit(chart$Genes, ', ');
	for (term in 1:length(termsIDs)) {
		tmp = unlist(IDs[IDs %in% termsIDs[[term]]]);
		syms = names(tmp);
		names(syms) = tmp;	
		chart_new$Gene.Syms[term] = paste(tmp[match(termsIDs[[term]], names(tmp))], collapse=', ');		
	}
	return(chart_new);
}