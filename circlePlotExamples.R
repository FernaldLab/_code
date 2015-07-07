exn.circlePlotGenes = function (genes, data, colors, type='signed', power=18, returnAdj=F, linecol.base=.9, linecol.gamma=1, cex.labels=c(.5,1.2), min.cex.points=1, max.cex.points=2, orderBy=NULL) {
	if (!is.null(orderBy) & length(orderBy)!=length(genes)) {
		stop('Genes and orderBy vectors must be of equal length');
	}
	if (length(colors)==length(genes)) {
		colors = colors;
	} else if (length(colors)==ncol(data)) {
		colors = colors[match(genes, names(data))];
	} else {
		stop('Colors vector must be same length as genes vector or number of columns in data');
	}
	gdata = data[, match(genes, names(data))];
	adj = adjacency(gdata, type=type, power=power);
	if (is.null(orderBy)) {
		orderBy = order(-apply(adj,1,sum));
	} else {
		orderBy = as.numeric(orderBy);
	}
	circlePlot(adj, colnames(adj), orderBy,
			   lineColors=grey2red(100, linecol.base, linecol.gamma),
			   pointBg=colors,
			   max.cex.labels=cex.labels[2], min.cex.labels=cex.labels[1],
			   min.cex.points=min.cex.points, max.cex.points=max.cex.points
			   );
	if (returnAdj) {
		return(adj);
	}
}

exn.circlePlotCompareGenesAcrossNetworks =  function (genes, data1, data2, colors1, colors2, type='signed', power=18, linecol.base=.9, linecol.gamma=1, cex.labels=c(.5,1.2), min.cex.points=1, max.cex.points=2) {
	gdata1 = data1[, match(genes, names(data1))];
	gdata2 = data2[, match(genes, names(data2))];
	adj1 = adjacency(gdata1, type=type, power=power);
	adj2 = adjacency(gdata2, type=type, power=power);
	#colors1 = colors1[match(genes, names(data1))];
	#colors2 = colors2[match(genes, names(data2))];
	par(mfrow=c(1, 2));
	circlePlot(adj1, colnames(adj1), order(-apply(adj1,1,sum)),
		   	   lineColors=grey2red(100, linecol.base, linecol.gamma),
		  	   pointBg=colors1,
		   	   max.cex.labels=cex.labels[2], min.cex.labels=cex.labels[1],
		   	   min.cex.points=min.cex.points, max.cex.points=max.cex.points
		   	   );
	circlePlot(adj2, colnames(adj2), order(-apply(adj1,1,sum)),
	   	   	   lineColors=grey2red(100, linecol.base, linecol.gamma),
	  	  	   pointBg=colors2,
	   	       max.cex.labels=cex.labels[2], min.cex.labels=cex.labels[1],
	   	       min.cex.points=min.cex.points, max.cex.points=max.cex.points
	   	       );
}
