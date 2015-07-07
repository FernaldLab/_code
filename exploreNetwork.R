########################################################################
#####
##### Parsers for cichlid gene/transcript id formats
#####
########################################################################
#############   	extract scaffold and gene #'s from different cichlid id formats
# Depends on:		NA 
#
## Arguments
### transcripts:	character vector of transcript ids formatted as species.'mrna'.scaffold#.gene#.transcript#
### genes:			character vector of gene ids formatted as species.'gene'.scaffold#.gene#
### subfeature_id:	character vector of subfeature ids (e.g. exons) formatted as species.subfeature.scaffold#.gene#.transcript#.subfeature#
### prefix:			character string representing species and subfeature type, e.g. mz.exon
#
## Value
### returns character string formatted as scaffold#.gene#

exn.mRNAParser = function(transcripts) {
	return(gsub('.[0-9]+$', '', gsub('mz.mrna.', '', transcripts)));
}

exn.geneParser = function(genes) {
	return(gsub('mz.gene.', '', genes));
}

exn.subfeatureParser = function(subfeature_id, prefix='mz.exon.') {
	return(gsub('.[0-9]+$', '', gsub('.[0-9]+$', '', gsub(prefix, '', subfeature_id))));
}

########################################################################
#############   parse cichlid ids of different formats
# Depends on:	NA 
#
## Arguments
### ids:		character vector of ids, each id containing up to five fields separated by common delimiter
### splitOn:	delimiter separating fields in each id, '.' by default
#
## Value
### mat:		character matrix with 5 columns (species, feature, geneID, transcriptID, subfeatureID), one row per id

exn.generalParser = function(ids, splitOn='\\.') {
	mat = matrix(nrow=length(ids), ncol=5);
	colnames(mat) = c('species', 'feature', 'geneID', 'transcriptID', 'subfeatureID');
	tmp = strsplit(ids, splitOn);
	if (max(as.numeric(unlist(lapply(tmp, length)))) > 6) {
		warning('BE CAREFUL, MORE SUBFEATURE FILEDS THAN USUAL');
	}
	#mat[, 1] = unlist(lapply(tmp, function(vec){vec[2]} ));
	for (i in 1:length(tmp)) {
		this = tmp[[i]];
		mat[i, 1] = this[1];
		mat[i, 2] = this[2];
		mat[i, 3] = paste(this[3:4], collapse='.');
		if (mat[i, 2]=='gene') {
			mat[i, 4:5] = NA;
		} else if (mat[i, 2]=='mrna') {
			mat[i, 4] = paste(this[3:5], collapse='.');
			mat[i, 5] = NA;
		} else {
			mat[i, 4] = paste(this[3:5], collapse='.');
			mat[i, 5] = paste(this[3:length(this)], collapse='.');
		}		
	}
	return(mat);
}

########################################################################
#####
##### Organize gene/transcript ids, annotations, and Entrez IDs by module or for whole network
#####
########################################################################
#############   get names of genes in each module
# Depends on:	NA 
#
## Arguments
### DATA:		data frame of numeric expression values, rows=samples, cols=genes/transcripts
### colors:		character vector of module colors, length(colors) must equal ncol(DATA)		
#
## Value
### modGenes:	list of character vectors containing module gene ids	

exn.getModuleGenes = function(DATA, colors) {
	modNames = names(table(colors));
	modGenes = list();
	for (mod in 1:length(modNames)) {
		modGenes[[mod]] = names(DATA)[colors==modNames[mod]];
		names(modGenes)[mod] = modNames[mod];
	}
	return(modGenes);
}

########################################################################
#############   	get annotations for genes in each module
# Depends on:		exn.getModuleGenes()
#					exn.mRNAParser()
#					exn.geneParser()
#
## Arguments
### DATA:			data frame of numeric expression values, rows=samples, cols=genes/transcripts
### colors:			character vector of module colors, length(colors) must equal ncol(DATA)	
### annos:			cichlid annotation table, first column must be gene ids formatted as species.'gene'.scaffold#.gene#
#
## Value
### modGeneAnnos:	list where each element is a data frame holding annotations for a given module

exn.getModuleGeneAnnotations = function(DATA, colors, annos, addCleanSymCol=T, symCol=3) {
	modGenes = exn.getModuleGenes(DATA, colors);
	check_mRNA = sum(grepl('mrna',unlist(modGenes))) == length(unlist(modGenes));
	if (check_mRNA) {
		parser = 'exn.mRNAParser';
	} else {
		parser = 'exn.geneParser';
	}
	modGeneAnnos = list();
	for (mod in 1:length(modGenes)) {
		modGenesParsed = eval(call(parser, modGenes[[mod]]));
		#modGeneAnnos[[mod]] = annos[exn.geneParser(annos[,1]) %in% modGenesParsed, ];print(modGenesParsed)
		modGeneAnnos[[mod]] = annos[match(modGenesParsed, exn.geneParser(annos[,1])), ];
		names(modGeneAnnos)[mod] = names(modGenes)[mod];
		names(modGeneAnnos[[mod]]) = c('geneID', 'ensemblID', 'geneSymbol', 'geneName');
		#modGeneAnnos[[mod]] = modGeneAnnos[[mod]][!(is.na(modGeneAnnos[[mod]][, 1])), ];
		if (addCleanSymCol) {
			tmp = toupper(as.character(modGeneAnnos[[mod]][,symCol]));
			tmp = strsplit(gsub('\\s', '', tmp), '\\(');
			tmpclean = c();
			for (gene in 1:length(tmp)) {
				if (length(tmp[[gene]]) == 0) {
					tmpclean = c(tmpclean, NA);
					msg = paste('MISSING GENE SYMBOL IN ', names(modGenes)[mod], ': ', modGeneAnnos[[mod]][gene, 1], sep='');
					warning(msg);
					next;
				} else {
					tmpclean = c(tmpclean, tmp[[gene]][1]);#print(tmpclean)
				}
			}
			modGeneAnnos[[mod]] = cbind(modGeneAnnos[[mod]], geneSymbolClean=tmpclean);
		}
	}
	return(modGeneAnnos);
}

########################################################################
#############   	get gene symbols for genes in each module
# Depends on:		exn.getModuleGeneAnnotations()
#
## Arguments
### DATA:			data frame of numeric expression values, rows=samples, cols=genes/transcripts
### colors:			character vector of module colors, length(colors) must equal ncol(DATA)	
### annos:			cichlid annotation table, first column must be gene ids formatted as species.'gene'.scaffold#.gene# 
### symCol:			number indicating column that holds gene symbols in annotation table
#
## Value
### modGeneSyms:	list where each element is a character vector of gene symbols for a given module

exn.getModuleGeneSymbols = function(DATA, colors, annos, modGeneAnnos=NULL, addCleanSymCol=NULL, symCol=3) {
	if (!(is.null(modGeneAnnos)) & any(grepl('geneSymbolClean', names(modGeneAnnos[[1]])))) {
		cat('...using provided modGeneAnnos to get gene symbols\n');
		modGeneSyms = list();
		for (mod in 1:length(modGeneAnnos)) {
			modGeneSyms[[mod]] = modGeneAnnos[[mod]]$geneSymbolClean;
			names(modGeneSyms)[mod] = names(modGeneAnnos)[mod];
		}
	} else {
		cat('...running exn.getModuleGeneAnnotations to get gene symbols\n');
		modGeneAnnos = exn.getModuleGeneAnnotations(DATA, colors, annos);
		modGeneSyms = list();
		for (mod in 1:length(modGeneAnnos)) {
			tmp = toupper(as.character(modGeneAnnos[[mod]][,symCol]));
			tmp = strsplit(gsub('\\s', '', tmp), '\\(');
			tmpclean = c();
			for (gene in 1:length(tmp)) {
				if (length(tmp[[gene]]) == 0) {
					next;
				} else {
					tmpclean = c(tmpclean, tmp[[gene]][1]);
				}
			}
			modGeneSyms[[mod]] = tmpclean;
			names(modGeneSyms)[mod] = names(modGeneAnnos)[mod];
		}
	}
	return(modGeneSyms);
}

########################################################################
#############   	translate list of module gene symbols to human Entrez IDs 
# Depends on:		library(org.Hs.eg.db); IDs=as.list(org.Hs.egSYMBOL2EG);
#
## Arguments
### modGeneSyms:	list where each element is a character vector of gene symbols for a given module
### IDs:			list where gene symbols are names of elements and Entrez IDs are values
#
## Value
### modIDs:			list where each element is a character vector of Entrez IDs for a given module

exn.moduleGeneSymbolsToHumanEntrezIDs = function(modGeneSyms, IDs) {
	modIDs = list();
	for (mod in 1:length(modGeneSyms)) {
		modIDs[[mod]] = unlist(IDs[names(IDs) %in% modGeneSyms[[mod]]]);
		names(modIDs)[mod] = names(modGeneSyms)[mod];
	}
	return(modIDs);
}

########################################################################
#############   get annotations, gene symbols, and Entrez IDs for all genes
# Depends on:	exn.mRNAParser()
#				exn.geneParser()
#				library(org.Hs.eg.db); IDs=as.list(org.Hs.egSYMBOL2EG);
#
## Arguments
### DATA:		data frame of numeric expression values, rows=samples, cols=genes/transcripts
### annos: 		cichlid annotation table, first column must be gene ids formatted as species.'gene'.scaffold#.gene# 
### IDs:		list where gene symbols are names of elements and Entrez IDs are values
### write:		boolean indicating whether to write .txt files of gene symbols and IDs
### path:		path of directory to write to (need trailing /), which must already exist, or by default writes to working directory
### suffix:		character string to be appended to filename for each module, '' by default
### symCol:		number indicating column that holds gene symbols in annotation table
#
## Value
### list with 3 elements: data frame of annotations for all genes, character vector of all gene symbols, character vector of all gene Entrez IDs

exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork = function(DATAorNames, annos, IDs, write=T, path='', suffix='', symCol=3) {
	if (is.vector(DATAorNames)) {
		allGenes = DATAorNames;
	} else {
		allGenes = names(DATAorNames);
	}
	#print(head(allGenes))
	check_mRNA = sum(grepl('mrna',allGenes)) == length(allGenes);
	if (check_mRNA) {
		parser = 'exn.mRNAParser';
	} else {
		parser = 'exn.geneParser';
	}
	allGenesAnnos = annos[exn.geneParser(annos[,1]) %in% eval(call(parser, allGenes)),];#print(head(allGenesAnnos))
	names(allGenesAnnos) = c('geneID', 'ensemblID', 'geneSymbol', 'geneName');
	allGenesSym = toupper(as.character(allGenesAnnos[, symCol]));
	allGenesSym = strsplit(gsub('\\s', '', allGenesSym), '\\(');
	tmp = c();
	for (gene in 1:length(allGenesSym)) {
		if (length(allGenesSym[[gene]])==0) {
			next;
		} else {
			tmp = c(tmp, allGenesSym[[gene]][1]);
		}
	}; 
	allGenesSym = tmp;
	allGenesAnnos = cbind(allGenesAnnos, geneSymbolClean=allGenesSym);
	allGenesIDs = unlist(IDs[names(IDs) %in% allGenesSym]);
	# IDvec = c();
	# for (gene in 1:length(allGenesIDs)) {
		# IDvec = c(IDvec, allGenesIDs[[gene]]);
	# }
	if (write) {
		write.table(allGenesSym, file=paste(path, 'allGeneSymbols', suffix, '.txt', sep=''), 
					quote=F, row.names=F, col.names=F);
		write.table(allGenesIDs, file=paste(path, 'allGeneEntrezIDs', suffix, '.txt', sep=''), 
					quote=F, row.names=F, col.names=F);
	}
	return(list(geneAnnos=allGenesAnnos, geneSym=allGenesSym, geneIDs=allGenesIDs));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.reOrderEntrez = function(vecSyms, vecIDs) {
	return (vecIDs[match(vecSyms, names(vecIDs))]);
}

########################################################################
#############   
# Depends on:	NA
#
## Arguments
### IDvec:
### IDs:
#
## Value
### out:

exn.EntrezIDsToGeneSymbols = function(IDvec, IDlist=IDs) {
	tmp = unlist(IDs[IDs %in% IDvec]);
	out = names(tmp);
	names(out) = tmp;
	return(out);
}

########################################################################
#############   get genes and their annotations, gene symbols, and Entrez IDs for each module
# Depends on:	exn.getModuleGenes()
#				exn.getModuleGeneAnnotations()
#				exn.getModuleGeneSymbols()
#				exn.moduleGeneSymbolsToHumanEntrezIDs()
#
## Arguments
### DATA:		data frame of numeric expression values, rows=samples, cols=genes/transcripts
### colors:		character vector of module colors, length(colors) must equal ncol(DATA)	
### annos:		cichlid annotation table, first column must be gene ids formatted as species.'gene'.scaffold#.gene#
### IDs:		list where gene symbols are names of elements and Entrez IDs are values
#
## Value
### list with 4 elements: $genes = character vector of gene names in each module, $geneAnnos = list of data frames containing annotation information for each module, $geneSyms = character vector of gene symbols in each module, $geneIDs = character vector of Entrez IDs in each module

exn.buildModuleGeneInfo = function(DATA, colors, annos, IDs, ...) {
	modGenes = exn.getModuleGenes(DATA, colors);
	modGeneAnnos = exn.getModuleGeneAnnotations(DATA, colors, annos, ...);
	modGeneSyms = exn.getModuleGeneSymbols(DATA, colors, annos, ...);
	#modGeneSyms = exn.getModuleGeneSymbols(modGeneAnnos=modGeneAnnos, ...);
	modGeneIDs = exn.moduleGeneSymbolsToHumanEntrezIDs(modGeneSyms, IDs);
	for (mod in 1:length(modGeneIDs)) {
		modGeneIDs[[mod]] = exn.reOrderEntrez(modGeneSyms[[mod]], modGeneIDs[[mod]]);
	}
	return(list(genes=modGenes, geneAnnos=modGeneAnnos, geneSyms=modGeneSyms, geneIDs=modGeneIDs));
}

########################################################################
#############   		examine numbers of transcripts, genes, annotations, gene symbols, and Entrez IDs in each module
# Depends on:			exn.mRNAParser()
#						output from exn.buildModuleGeneInfo()
#
## Arguments
### moduleGeneInfoList:	list structured as output from exn.buildModuleGeneInfo()
#
## Value
### ratioTable:	matrix with numbers of transcripts, genes, annotated genes, gene symbols, and Entrez IDs (columns) in each module (rows)

exn.buildModuleGeneInfoRatioTable = function(moduleGeneInfoList) {
	if (any(names(moduleGeneInfoList) != c('genes', 'geneAnnos', 'geneSyms', 'geneIDs'))) {
		stop('moduleGeneInfoList must be output from exn.buildModuleGeneInfo() with elements named $genes, $geneAnnos, $geneSyms, and $geneIDs');
	}
	modGenes = moduleGeneInfoList$genes;
	modGeneAnnos = moduleGeneInfoList$geneAnnos;
	modGeneSyms = moduleGeneInfoList$geneSyms;
	modGeneIDs = moduleGeneInfoList$geneIDs;
	modLengths = unlist(lapply(modGenes,length));
	modGenesParsed = lapply(modGenes, exn.mRNAParser);
	modUniqueLengths = unlist(lapply(lapply(modGenesParsed, unique), length));
	# modAnnoDims = lapply(modGeneAnnos, dim);
	# modAnnoRows = c();
	# for (mod in 1:length(modAnnoDims)) {
		# modAnnoRows = c(modAnnoRows, modAnnoDims[[mod]][1]);
	# }
	modAnnoRows = c();
	for (mod in 1:length(modGeneAnnos)) {
		modAnnoRows = c(modAnnoRows, sum(!(is.na(modGeneAnnos[[mod]][, 1]))));
	}
	#modSymLengths = unlist(lapply(modGeneSyms, length));
	modSymLengths = modAnnoRows;		# maybe remove altogether
	# modIDLengths = unlist(lapply(modGeneIDs, length));
	modIDLengths = c();
	for (mod in 1:length(modGeneIDs)) {
		modIDLengths = c(modIDLengths, sum(!(is.na(modGeneIDs[[mod]]))));	
	}
	ratioTable = matrix(nrow=length(modLengths), ncol=5);
	colnames(ratioTable) = c('num_transcripts', 'num_genes', 'annotated', 'symbols', 'HsEntrezIDs');
	rownames(ratioTable) = names(modLengths);
	ratioTable[, 1] = modLengths;
	ratioTable[, 2] = modUniqueLengths;
	ratioTable[, 3] = modAnnoRows;
	ratioTable[, 4] = modSymLengths;
	ratioTable[, 5] = modIDLengths;
	return(ratioTable);
}

########################################################################
#############   		view genes, annotations, gene symbols, and Entrez IDs for one module
# Depends on:			output from exn.buildModuleGeneInfo()
#
## Arguments
### moduleGeneInfoList:	list structured as output from exn.buildModuleGeneInfo()
### module:				character string naming module to view
#
## Value
### out:				list containing output from exn.buildModuleGeneInfo() for a single module

exn.viewMod = function(moduleGeneInfoList, module) {
	check = names(moduleGeneInfoList) == c('genes', 'geneAnnos', 'geneSyms', 'geneIDs');
	if (!all(check)) {
		stop('WRONG INPUT DATA \n  MUST USE OUTPUT FROM buildModuleGeneInfo()');
	}
	ind = match(module, names(moduleGeneInfoList$genes));
	out = list();
	for (m in 1:4) {
		#cat('ind=',ind,'\n');cat('m=',m,'\n');cat('module=',module,'\n'); cat('names(mod[[m]])[ind]',names(mod[[m]])[ind],'\n')
		if (!(names(moduleGeneInfoList[[m]])[ind] == module)) {
			stop('MODULE NAMES DO NOT MATCH ACROSS SUBLISTS');
		}
		out[[m]] = moduleGeneInfoList[[m]][[ind]];
	}
	#out = out[2:5];
	names(out) = c('genes', 'geneAnnos', 'geneSyms', 'geneIDs');
	return(out);
} 

########################################################################
#############   	get distributions of transcripts across modules, organized by gene
# Depends on:		exn.mRNAParser()
#
## Arguments
### transcripts:	character vector of transcript ids formatted as species.'mrna'.scaffold#.gene#.transcript#
### colors:			character vector of module colors
### sorted:			boolean indicating whether to sort genes in output by number of transcripts 
#
## Value
### list with 3 elements: $counts = numeric vector of transcript counts for each gene in the network, $modules = list of character vectors containing modules of transcripts for a given gene, $isoforms = list of character vectors containing names of transcripts for each gene

exn.transcriptsAcrossModules = function(transcripts, colors, sorted=T) {
	parsedTranscripts = exn.mRNAParser(transcripts);
	if (any(grepl('A|B$', parsedTranscripts))) {
		parsedTranscripts = gsub('.[0-9]+[AB]$', '', parsedTranscripts);
	}
	transcriptTable = table(parsedTranscripts);
	if (sorted) {
		transcriptTable = sort(transcriptTable, decreasing=T);
	}
	transcriptMods = list();
	transcriptsByGene = list();
	for (tr in 1:length(transcriptTable)) {
		transcriptMods[[tr]] = colors[parsedTranscripts==names(transcriptTable)[tr]];
		transcriptsByGene[[tr]] = transcripts[parsedTranscripts==names(transcriptTable)[tr]];
	}
	names(transcriptMods) = names(transcriptTable);
	names(transcriptsByGene) = names(transcriptTable);
	return(list(counts=transcriptTable, modules=transcriptMods, isoforms=transcriptsByGene));
}

########################################################################
#####
##### File input/output - write modules to files, read gtf/gff files, read DAVID output
#####
########################################################################
#############   write names/symbols/IDs of module genes to txt files
# Depends on:	NA
#
## Arguments
### modList:	list where each element is a character vector for a given module
### path:		path of directory to write to (need trailing /), which must already exist, or by default writes to working directory
### suffix:		character string to be appended to filename for each module, '' by default
### combine:	boolean indicating whether to write additional tab-delimited file containing one module per column
#
## Value
### filenames:	character vector of filenames for each module

exn.writeModules = function(modList, path='', suffix='', combine=F, ...) {
	filenames = c();
	for (mod in 1:length(modList)) {
		thisfile = paste(path, names(modList)[mod], suffix, '.txt', sep='');
		filenames = c(filenames, thisfile);
		write.table(modList[[mod]], file=thisfile, quote=F, row.names=F, col.names=F, ...);
	}
	if (combine) {
		tmplengths = unlist(lapply(modList, length));
		tmpMat = matrix(nrow=max(tmplengths), ncol=length(tmplengths));
		colnames(tmpMat) = names(tmplengths);
		for (mod in 1:ncol(tmpMat)) {
			thismod = modList[[mod]];
			diff = nrow(tmpMat) - length(thismod);
			if (diff>0) {
				thismod = c(thismod, rep('', diff));
			} 
			tmpMat[, mod] = thismod;
		}
		combofilename = paste(path, 'allMods', suffix, '.txt', sep='');
		write.table(tmpMat, file=combofilename, quote=F, row.names=F, sep='\t');
	}
	return(filenames);
}

########################################################################
#############   get transcripts from network that are in a given gtf/gff file
# Depends on:	NA
#
## Arguments
### gtf_file:	character string name of input feature file
### IDfield:	character string name of id type of interest, usually 'transcript_id' or 'gene_id'
### fieldNum:	number of column containing ids in feature file
### uniqueOnly:	boolean
#
## Value
### list with 2 elements: $gtf = input file stored as data frame, $hits = names of genes/transcripts that hit in feature file

exn.getIDsFromGTF = function(gtf_file, IDfield='transcript_id', fieldNum=9, uniqueOnly=T) {
	gtf = read.table(gtf_file, header=F, sep='\t');
	hits = c();
	for (row in 1:nrow(gtf)) {
		thisrow = strsplit(as.character(gtf[row, fieldNum]), ';')[[1]];
		ind = grep(IDfield, thisrow);
		hit = gsub('\\s', '', gsub(IDfield, '', thisrow[ind]));
		hits = c(hits, hit);
	}
	if (uniqueOnly) {
		hits = unique(hits);
	}
	return(list(gtf=gtf, hits=hits));
}

########################################################################
#############   get subset of gtf/gff file containing rows for gene/transcript
# Depends on:	NA
#
## Arguments
### gtf:		data frame containing gtf/gff file, usually get from output of exn.getIDsFromGTF() 
### ID:			character string transcript/gene id 
### fieldNum:	number of column containing ids in gtf
#
## Value
### data frame containing rows of input gtf that contain gene/transcript id

exn.getGTFRowsFromID = function(gtf, ID, fieldNum=9) {
	IDfields = strsplit(gtf[, fieldNum], ';');
	keepMe = c();
	for (row in 1:length(IDfields)) {
		if (any(grepl(ID, IDfields[[row]]))) {
			keepMe = c(keepMe, row);
		}
	}
	return(gtf[keepMe, ]);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# expects files to be formatted as output from exn.getGOfromEntrezByGene()
exn.readGeneGOfromFilesIntoList = function(dir) {
	geneGO = list();
	files = list.files(dir);
		for (f in 1:length(files)) {
			geneGO[[f]] = read.table(paste(dir, '/', files[f], sep=''), sep='\t', header=T);
			names(geneGO)[f] = gsub('entrez', '', gsub('_GOterms.txt', '', files[f]));
		}	
	return(geneGO);
}

########################################################################
#####
##### Eigengene and kME related functions
#####
########################################################################
#############   compute eigengene based connectivity (kME) for all genes/transcripts in all modules
# Depends on:	WGCNA::orderMEs()
#				WGCNA::corPvalueStudent()
#
## Arguments
### DATA:		data frame of numeric expression values, rows=samples, cols=genes/transcripts
### MEs:		data frame of module eigengenes, as computed by WGCNA::blockwiseModules() or WGCNA::moduleEigengenes()
### combine:
#
## Value
### list with 2 elements: $k = data frame of kME for each network member (rows) in each module (columns), $pval = data frame of pvalues corresponding to values in $k

exn.computekME = function(DATA, MEs, modules=F, colors=NULL) {
	MEs = orderMEs(MEs);
	kME = as.data.frame(cor(DATA, MEs, use='p'));
	names(kME) = paste('k', names(kME), sep='');
	kMEpval = as.data.frame(corPvalueStudent(as.matrix(kME), nrow(DATA)));
	names(kMEpval) = paste('p.', names(kMEpval), sep='');
	
	out = as.data.frame(matrix(nrow=nrow(kME), ncol=(2*ncol(kME))));
	kcols = seq(1, ncol(out), 2); 
	pcols = kcols+1;
	out[, kcols] = kME;
	names(out)[kcols] = names(kME);
	out[, pcols] = kMEpval;
	names(out)[pcols] = names(kMEpval);
	rownames(out) = rownames(kME);
	
	modList = NULL;
	if (modules) {
		modList = list();
		modNames = names(table(colors));
		for (mod in 1:length(modNames)) {
			modList[[mod]] = out[colors==modNames[mod], ];
			names(modList)[mod] = modNames[mod];
		}
	}

	return(list(all=out, mods=modList));
}

########################################################################
#############   	plot eigengene dendrogram and heatmap with predefined margins and labels
# Depends on:		WGCNA::plotEigengeneNetworks()
#					WGCNA::blueWhiteRed()
#
## Arguments
### MEs:			data frame of module eigengenes, as computed by WGCNA::blockwiseModules() or WGCNA::moduleEigengenes()
### heatmapColors:	character vector of colors, e.g. those returned by WGCNA::blueWhiteRed()
#
## Value
### besides plotting in Quartz window, returns matrix of ME-ME correlations

exn.plotEigengeneNetworks2 = function(MEs, heatmapColors=blueWhiteRed(50), returnCors=F) {
	plotEigengeneNetworks(orderMEs(MEs), 
						  setLabels='', marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2), 
						  heatmapColors=heatmapColors);
	if (returnCors) {
		return(cor(MEs));
	}
}

########################################################################
#############   	get kMEs and pvals across all modules for a given transcript
# Depends on:		NA
#
## Arguments
### transcript_id:	character string naming gene/transcript, must match rowname in kMEtable and kMEtablePval
### kMEtable:		data frame containing kME values for all network members (rows) in all modules (columns), e.g. exn.computekME$k
### kMEtablePval:	data frame containing pvals that correspond to kMEtable, e.g. exn.computekME$pval
#
## Value
### mat:			numeric matrix where rows are modules and columns are kME values and pvals 

exn.getTranscriptkME = function(transcript_id, kMEtable, kMEtablePval=NULL) {
	if (is.null(kMEtablePval)) {
		kcols = seq(1, ncol(kMEtable), 2);
		mat = matrix(nrow=(ncol(kMEtable)/2), ncol=2, dimnames=list(gsub('kME', '', colnames(kMEtable))[kcols], c('kME','pval')));
		mat[, 1] = t(kMEtable[rownames(kMEtable)==transcript_id, kcols]);
		mat[, 2] = t(kMEtable[rownames(kMEtable)==transcript_id, kcols+1]);
	} else {
		mat = matrix(nrow=ncol(kMEtable), ncol=2, dimnames=list(gsub('kME', '', colnames(kMEtable)), c('kME','pval')));
		mat[, 1] = t(kMEtable[rownames(kMEtable)==transcript_id, ]);
		mat[, 2] = t(kMEtablePval[rownames(kMEtablePval)==transcript_id, ]);
	}
	return(mat);
}

########################################################################
#############   	get kME and pval for a given transcript in a module
# Depends on:		NA
#
## Arguments
### transcript_id:	character string naming gene/transcript, must match rowname in kMEtable and kMEtablePval
### DATA:			data frame of numeric expression values, rows=samples, cols=genes/transcripts
### colors:			character vector of module colors, length(colors) must equal ncol(DATA)
### kMEtable:		data frame containing kME values for all network members (rows) in all modules (columns), e.g. exn.computekME$k
### kMEtablePval:	data frame containing pvals that correspond to kMEtable, e.g. exn.computekME$pval
### module:			optional character string specifying module
#
## Value
### out:			numeric vector containing kME and pval for given transcript, either in its assigned module (if module=NULL) or a specified module

exn.getTranscriptModulekME = function(transcript_id, DATA, colors, kMEtable, kMEtablePval=NULL, module=NULL) {
	if (is.null(module)) {
		tmod = colors[names(DATA)==transcript_id];
	} else {
		tmod = module;
	}
	row = rownames(kMEtable)==transcript_id;
	if (is.null(kMEtablePval)) {
		#columns = grep(tmod, names(kMEtable));
		columns = which(gsub('kME|p.kME', '', names(kMEtable)) == tmod);
		out = c(kMEtable[row, columns[1]], kMEtable[row, columns[2]]);
		names(out) = names(kMEtable)[columns];
	} else {
		column = gsub('kME', '', colnames(kMEtable))==tmod;
		kMEval = kMEtable[row, column];
		kMEpval = kMEtablePval[row, column];
		out = c(kMEval, kMEpval);
		names(out) = c(colnames(kMEtable)[column], colnames(kMEtablePval)[column]);
	}
	return(out);
}

exn.getTranscriptModulekMEAll = function(transcript_ids, DATA, colors, kMEtable) {
	tmp = c();
	for (tr in 1:length(transcript_ids)) {
		tmp = c(tmp, exn.getTranscriptModulekME(transcript_ids[tr], DATA, colors, kMEtable)[1]);
	}
	return(tmp);
}
########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.getTopkMEFromEachModule = function(kMEtable, colors) {
	if (nrow(kMEtable)!=length(colors)) {
		stop();
	}
	modNames = names(table(colors));
	hubs = as.data.frame(matrix(nrow=length(modNames), ncol=2));
	rownames(hubs) = modNames;
	names(hubs) = c('id', 'kME');
	for (mod in 1:length(modNames)) {
		ind = match(paste('kME', modNames[mod], sep=''), names(kMEtable));
		tmpkME = kMEtable[colors==modNames[mod], ];
		tmpkME = tmpkME[order(tmpkME[, ind], decreasing=T), ];
		hubs[mod, ] = c(rownames(tmpkME)[1], tmpkME[1, ind]); 
	} 
	return(hubs);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getModulekME = function(module, colors, kMEtable, kMEpvals=NULL, orderByRank=T) {
	modGenes = colors==module;
	if (is.null(kMEpvals)) {
		#modCols = grep(module, names(kMEtable));#######
		modCols = which(gsub('ME|kME|p.kME', '', names(kMEtable)) == module);
		out = kMEtable[modGenes, modCols];
	} else {
		modCol = match(module, gsub('kME', '', names(kMEtable)));
		out = cbind(kMEtable[modGenes, modCol], kMEpvals[modGenes, modCol]);
		dimnames(out) = list(rownames(kMEtable)[modGenes], c(names(kMEtable)[modCol], names(kMEpvals)[modCol]));
	}
	if (orderByRank) {
		out = out[order(out[, 1], decreasing=T), ];
	}
	return(out);
}

########################################################################
#####
##### Simulate and work with allele specific networks
#####
########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.countsToFreqs = function(countMat) {
	if (ncol(countMat)!=2) {
		stop('ONLY 2 ALLELES ALLOWED');
	}
	freqMat = countMat;
	for (row in 1:nrow(countMat)) {
		f1 = countMat[row, 1]/sum(countMat[row, ]);
		f2 = countMat[row, 2]/sum(countMat[row, ]);
		freqMat[row, ] = c(f1, f2);
	}
	return(list(counts=countMat, freqs=freqMat));
}

########################################################################

# exn.simulateFreqs = function(DATA, nTranscripts=floor(.05*ncol(DATA)), subjRows) {
	# tsb = list();
	# tr = sample(names(DATA), nTranscripts);
	# for (this_tr in 1:length(tr)) {
		# mat = matrix(nrow=length(tr), ncol=2, dimnames=list(c(), c()));
		# rownames(mat)[this_tr] = tr[this_tr];
		# colnames(mat) = c('A', 'B');
	# }
	# return(mat);
# }

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.applyAlleleFreqsToTranscript = function(transcript_id, DATA, freqs=matrix(rep(.5, 4), ncol=2), rSamples=c(5, 6)) {
	ind = match(transcript_id, names(DATA));
	newExpr = data.frame(A=DATA[, ind], B=DATA[, ind]);
	dimnames(newExpr) = list(c(rownames(DATA)), c(paste(transcript_id, names(newExpr), sep='')));
	if (any(is.na(newExpr[rSamples, ]))) {
		warning('NAs detected for ', transcript_id);
	}
	newExpr[rSamples, ] = newExpr[rSamples, ] * freqs;
	if (ind==1) {
		DATA.new = as.data.frame(cbind(newExpr, DATA[, 2:ncol(DATA)]));
	} else if (ind==ncol(DATA)) {
		DATA.new = as.data.frame(cbind(DATA[, 1:(ncol(DATA)-1)], newExpr));
	} else {
		DATA.new = as.data.frame(cbind(DATA[, 1:(ind-1)], newExpr, DATA[, (ind+1):ncol(DATA)]));
	}
	return(DATA.new);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# transcriptsBySubjectList is list where each component is m x 2 df/matrix of allele frequencies with rownames corresponding to colnames in DATA (transcripts) and columns correspond to A and B allele frequencies, list component names must match rownames in DATA (subject ids)
exn.applyAlleleFreqsToMutipleTranscripts = function(transcriptsBySubjectList, DATA, createPaddedDATA=T) {
	tbs = transcriptsBySubjectList;
	tr_ids = unique(rapply(tbs, rownames));
	tr_ids = tr_ids[order(match(tr_ids, names(DATA)))];
	trList = list(); 
	for (tr in 1:length(tr_ids)) { 
		nrow = sum(rapply(tbs, rownames)==tr_ids[tr]);
		tmpFreqs = matrix(nrow=nrow, ncol=2, dimnames=list(rep('x', nrow), c('A', 'B')));
		tmpFreqsRow = 1;
		for (sub in 1:length(tbs)) {
			for (row in 1:nrow(tbs[[sub]])) { 
				if (tr_ids[tr]==rownames(tbs[[sub]])[row]) {
					tmpFreqs[tmpFreqsRow, ] = tbs[[sub]][match(tr_ids[tr], rownames(tbs[[sub]])), ];
					rownames(tmpFreqs)[tmpFreqsRow] = names(tbs)[sub];
					tmpFreqsRow = tmpFreqsRow + 1; 
				}
			}
		}
		trList[[tr]] = tmpFreqs;
		names(trList)[tr] = tr_ids[tr];
	}
	tmpDATA = DATA;
	for (tr in 1:length(trList)) {
		rSamples = match(rownames(trList[[tr]]), rownames(tmpDATA));
		tmpDATA = exn.applyAlleleFreqsToTranscript(transcript_id=names(trList)[tr], DATA=tmpDATA, freqs=trList[[tr]], rSamples=rSamples);
	}
	padDATA = NULL;
	if (createPaddedDATA) {
		padtrList = trList;
		padDATA = DATA;
		for (tr in 1:length(padtrList)) {
			for (row in 1:nrow(padtrList[[tr]])) {
				padtrList[[tr]][row, ] = c(0.5, 0.5);
			}
			rSamples = match(rownames(padtrList[[tr]]), rownames(padDATA));
			padDATA = exn.applyAlleleFreqsToTranscript(transcript_id=names(padtrList)[tr], DATA=padDATA, freqs=padtrList[[tr]], rSamples=rSamples);
		}
	}
	return(list(DATA.new=tmpDATA, DATA.padded=padDATA, trList=trList));
}
 
########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.applyFreqsToRawExpr = function(DATA, freqMat, subjNum) {
	dat = DATA;
	freqMat = freqMat[order(match(rownames(freqMat), names(dat))), ];
	for (f in 1:nrow(freqMat)) {
		this_col = match(rownames(freqMat)[f], names(dat));#print(this_col)
		mat = cbind(dat[, this_col], dat[, this_col]);
		rownames(mat) = rownames(dat);
		colnames(mat) = c(paste(rownames(freqMat)[f], 'A', sep=''), paste(rownames(freqMat)[f], 'B', sep=''));
		mat[subjNum, ] = mat[subjNum, ] * freqMat[f, ];
		if (this_col==1) {
			dat = as.data.frame(cbind(mat, dat[, 2:ncol(dat)]));	
		} else if (this_col==ncol(dat)) {
			dat = as.data.frame(cbind(dat[, 1:(ncol(dat)-1)], mat));
		} else {
			dat = as.data.frame(cbind(dat[, 1:(this_col-1)], mat, dat[, (this_col+1):ncol(dat)]));
		}
	}
	return(list(DATA=DATA, newDATA=dat, freqMat=freqMat));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.applyFreqsToPreparedExpr = function(newDATA, freqMat, subjNum) {
	dat = newDATA;
	parsedNames = gsub('[AB]', '', names(dat));
	freqMat = freqMat[order(match(rownames(freqMat), parsedNames)), ];
	for (f in 1:nrow(freqMat)) {
		these_cols = which(parsedNames == rownames(freqMat)[f]);
		dat[subjNum, these_cols] = dat[subjNum, these_cols] * freqMat[f, ];
	}
	return(list(DATA=newDATA, newDATA=dat, freqMat=freqMat));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.assignAlleleFreqs = function(DATA, transcriptsBySubjectList) {
	dat = DATA;
	tsb = transcriptsBySubjectList;
	tsbFreqs = tsb;
	subjRows = match(names(tsb), rownames(dat));
	freqs = exn.countsToFreqs(tsb[[1]])$freqs;
	tsbFreqs[[1]] = freqs;
	dat = exn.applyFreqsToRawExpr(dat, freqs, subjRows[1])$newDATA;
	if (length(tsb)>1) {
		for (s in 2:length(tsb)) {
			freqs = exn.countsToFreqs(tsb[[s]])$freqs;
			tsbFreqs[[s]] = freqs;
			dat = exn.applyFreqsToPreparedExpr(dat, freqs, subjRows[s])$newDATA;
		}
	}
	return(list(DATA=DATA, newDATA=dat, counts=tsb, freqs=tsbFreqs));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# exn.padTranscriptsToMatchAlleleFreqOutput = function(DATA, DATA.alleles, transcripts, rSamples, divide=T) {
	# DATA.padded = DATA;
	# for (tr in 1:length(transcripts)) {
		# DATAcol = match(transcripts[tr], names(DATA.padded));
		# if (divide) {
			# new_tr = DATA.padded[, DATAcol];
			# new_tr[rSamples] = new_tr[rSamples] / 2;
		# } else {
			# new_tr = DATA.padded[, DATAcol];
		# }
		# DATA.padded = cbind(DATA.padded[, -DATAcol], new_tr, new_tr); 		
		# names(DATA.padded)[c(ncol(DATA.padded)-1, ncol(DATA.padded))] = paste(transcripts[tr], c('A', 'B'), sep='');
	# }
	# DATA.padded = DATA.padded[, match(names(DATA.alleles), names(DATA.padded))];
	# return(DATA.padded);
# }

exn.padTranscriptsToMatchAlleleFreqOutput = function(DATA.alleles, transcripts, rSamples, divide=T) {
	DATA.padded = DATA.alleles;
	if (all(grepl('A|B$', transcripts))) {
		transcripts = transcripts
	} else if (all(!grepl('A|B$', transcripts))) {
		tr1 = paste(tr, 'A', sep='');
		tr2 = paste(tr, 'B', sep='');
		tr.alleles = vector(mode='character', length=(2*length(transcripts)));
		tr.alleles[seq(1, length(tr.alleles), 2)] = tr1;
		tr.alleles[seq(2, length(tr.alleles), 2)] = tr2;
		transcripts = tr.alleles;
	} else {
		stop('something wrong with transcript ids');
	}
	inds = match(transcripts, names(DATA.padded));
	for (tr in seq(1, length(transcripts), by=2)) {
		for (row in rSamples) {
			total = sum(DATA.padded[row, inds[c(tr, tr+1)]]);
			if (divide) {
				DATA.padded[row, inds[c(tr, tr+1)]] = c(total/2, total/2);
			} else {
				DATA.padded[row, inds[c(tr, tr+1)]] = c(total, total);
			}
		}
	}
	return(DATA.padded);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.simulateAlleleFreqs = function(DATA, nTranscripts=floor(.01*ncol(DATA)), rSamples=c(5,6), jitter=0.2) {
	# randomly sample transcripts from data names and order them
	tr = sample(names(DATA), nTranscripts);
	checkNA = any(is.na(DATA[rSamples, names(DATA) %in% tr]));
	if (checkNA) {
		while (checkNA) {
			posNA = which(is.na(DATA[rSamples, names(DATA) %in% tr]), arr.ind=T)[, 2];
			tr[posNA] = sample(names(DATA), length(posNA)); 
			checkNA = any(is.na(DATA[rSamples, names(DATA) %in% tr]));
		}
	}
	tr = tr[order(match(tr, names(DATA)))];
	
	# randomly generate frequency pairs for each transcript in each specified sample
	
	# freqs = list();
	# for (f in 1:length(rSamples)) {
		# tmp = sample(1:100, nTranscripts);
		# freqs[[f]] = cbind(A=tmp, B=100-tmp)/100;
		# rownames(freqs[[f]]) = tr;
	# }
	# names(freqs) = rownames(DATA)[rSamples];
	
	tmp = sample(1:100, nTranscripts);
	freqs = list();
	freqs[[1]] = cbind(A=tmp, B=100-tmp)/100;
	rownames(freqs[[1]]) = tr;
	for (f in 2:length(rSamples)) {
		freqs[[f]] = freqs[[1]];
		for (ff in 1:nrow(freqs[[f]])) {
			dom = max(freqs[[f]][ff, ]);#print(dom)
			domCol = match(dom, freqs[[f]][ff, ]);
			range = dom * jitter;
			add = sample(seq(.01, range, length.out=100), 1) * sample(c(-1, 1), 1);
			domNew = dom + signif(add, 1);#print(domNew)
			if (domNew > 1) {
				domNew = 1 - abs(signif(add, 1));
			}
			if (domNew < (1-domNew)) {
				domNew = 1-domNew;
			}
			freqs[[f]][ff, domCol] = domNew;
			if (domCol==1) {
				freqs[[f]][ff, 2] = 1 - domNew;
			} else if (domCol==2) {
				freqs[[f]][ff, 1] = 1 - domNew;
			} else {
				stop('problem with column selection of dominant allele');
			}
		}
	}
	names(freqs) = rownames(DATA)[rSamples];
	
	# generate matrix (nrow(DATA) x 2) for each transcript with weighted expression values
	matList = list();
	for (t in 1:length(tr)) {
		# get transcript data and build matrix
		DATAcol = match(tr[t], names(DATA));
		this_tr = DATA[, DATAcol];
		dimnames = list(rownames(DATA), c(paste(tr[t], 'A', sep=''), paste(tr[t], 'B', sep='')));
		mat = matrix(nrow=nrow(DATA), ncol=2, dimnames=dimnames);
		mat[, 1] = this_tr;
		mat[, 2] = this_tr;
		# weight columns in matrix by simulated allele frequencies
		for (s in 1:length(rSamples)) {
			mat[rSamples[s], 1] = this_tr[rSamples[s]] * freqs[[s]][t, 1];
			mat[rSamples[s], 2] = this_tr[rSamples[s]] * freqs[[s]][t, 2];
		}
		matList[[t]] = mat;
		names(matList)[t] = tr[t];
	}
	matList0 = matList;
	DATAnew = DATA;
	if (names(matList)[1]==names(DATAnew)[1]) {
		DATAnew = as.data.frame(cbind(matList[[1]], DATAnew[, 2:ncol(DATAnew)]));
		matList = matList[-1];
	}
	if (names(matList)[length(matList)]==names(DATAnew)[ncol(DATAnew)]) {
		DATAnew = as.data.frame(cbind(DATAnew[, 1:(length(DATAnew)-1)], matList[[length(matList)]]));
		matList = matList[-(length(matList))];
	}
	if (length(matList)>0) {
		for (m in 1:length(matList)) {
			ind = match(names(matList)[m], names(DATAnew));
			DATAnew = as.data.frame(cbind(DATAnew[, 1:(ind-1)], matList[[m]], DATAnew[, (ind+1):ncol(DATAnew)]));
			if (ncol(DATAnew[, 1:(ind-1)])==1) {
				names(DATAnew)[1] = names(DATA)[1];
				warning('check first name of DATA.alleles')
			}
			if (ncol(DATAnew[, (ind+1):ncol(DATAnew)])==1) {
				names(DATAnew)[ncol(DATAnew)] = names(DATA)[ncol(DATA)];
				warning('check last name of DATA.alleles');
			}
		}
	}
	
	tr1 = paste(tr, 'A', sep='');
	tr2 = paste(tr, 'B', sep='');
	tr.alleles = vector(mode='character', length=(2*length(tr)));
	tr.alleles[seq(1, length(tr.alleles), 2)] = tr1;
	tr.alleles[seq(2, length(tr.alleles), 2)] = tr2;
	#DATApad = exn.padTranscriptsToMatchAlleleFreqOutput(DATA=DATA, DATA.alleles=DATAnew, transcripts=tr.alleles, rSamples=rSamples);
	DATApad = exn.padTranscriptsToMatchAlleleFreqOutput(DATA.alleles=DATAnew, transcripts=tr.alleles, rSamples=rSamples);
	return(list(DATA.padded=DATApad, DATA.alleles=DATAnew, tr=tr, tr.alleles=tr.alleles, freqs=freqs, matList=matList0));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getAlleleModsAndkME = function(DATA.padded, colors.padded, kME.padded, DATA.alleles, colors.alleles, kME.alleles) {
	if (sum(grepl('A|B$', names(DATA.padded)))<2 | sum(grepl('A|B$', names(DATA.alleles)))<2) {
		stop();
	}
	if (!all(grep('A|B$', names(DATA.padded)) == grep('A|B$', names(DATA.alleles)))) {
		stop();
	} else {
		inds = grep('A|B$', names(DATA.padded));
	}
	if (!(all(grep('A|B$', rownames(kME.padded))==inds)) | !(all(grep('A|B$', rownames(kME.alleles))==inds))) {
		stop();
	}
	kME = list(padded=kME.padded[inds, ], alleles=kME.alleles[inds, ]);
	colors = data.frame(padded=as.character(colors.padded[inds]), alleles=as.character(colors.alleles[inds]));
	rownames(colors) = names(DATA.padded)[inds];
	k.p = c(); k.a = c(); r.p = c(); c.p = c(); r.a = c(); c.a = c();
	for (tr in 1:nrow(colors)) {
		k.p = c(k.p, kME$padded[tr, match(colors$padded[tr], gsub('kME', '', names(kME$padded)))]);
		k.a = c(k.a, kME$alleles[tr, match(colors$alleles[tr], gsub('kME', '', names(kME$alleles)))]);
		
		tmpkME.p = kME.padded[colors.padded==colors$padded[tr], ];
		ind.p = match(colors$padded[tr], gsub('kME', '', names(tmpkME.p)));
		tmpkME.p = tmpkME.p[order(tmpkME.p[, ind.p], decreasing=T), ];
		r.p = c(r.p, (1:nrow(tmpkME.p))[rownames(tmpkME.p)==rownames(colors)[tr]]);
		c.p = c(c.p, r.p[tr]/nrow(tmpkME.p));
		
		tmpkME.a = kME.alleles[colors.alleles==colors$alleles[tr], ];
		ind.a = match(colors$alleles[tr], gsub('kME', '', names(tmpkME.a)));
		tmpkME.a = tmpkME.a[order(tmpkME.a[, ind.a], decreasing=T), ];
		r.a = c(r.a, (1:nrow(tmpkME.a))[rownames(tmpkME.a)==rownames(colors)[tr]]);
		c.a = c(c.a, r.a[tr]/nrow(tmpkME.a));
	}
	modkME = data.frame(mod.padded=as.character(colors$padded), mod.kME=as.numeric(k.p), rank=as.numeric(r.p), centile=as.numeric(c.p), mod.alleles=as.character(colors$alleles), mod.kME=as.numeric(k.a), rank=as.numeric(r.a), centile=as.numeric(c.a));
	names(modkME)[6:8] = c('mod.kME', 'rank', 'centile');
	rownames(modkME) = rownames(colors);
	return(list(kME=kME, colors=colors, modkME=modkME));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# exn.getAllelekMEs = function(transcripts, kMEs) {
	# out = kMEs[rownames(kMEs) %in% paste(transcripts, 'A', sep=''), ];
	# out = rbind(out, kMEs[rownames(kMEs) %in% paste(transcripts, 'B', sep=''), ]);
	# return(out[order(rownames(out)), ]);
# }

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# exn.getAlleleModColors = function(transcripts, colors, DATAorNames) {
	# if (is.vector(DATAorNames)) {
		# all_tr = DATAorNames;
	# } else {
		# all_tr = names(DATAorNames);
	# }
	# tr_colors = colors[match(transcripts, gsub('A|B$', '', all_tr))];
	# return(tr_colors);
# }

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getAlleleNetworks = function(DATA.padded, DATA.alleles, annos, IDs, matchLabels=T, networkType='signed', power=18, ...) {
	if (!(all(names(DATA.padded)==names(DATA.alleles)))) {
		stop('Both datasets must contain same transcripts in same order');
	} else {
		allele_names = names(DATA.padded)[grep('A|B$', names(DATA.padded))];
		tr_names = gsub('A|B', '', names(DATA.padded));
		#all = exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork(gsub('A|B', '', tr_names), annos, IDs, F);
		all = exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork(tr_names, annos, IDs, F);
	}
	cat('Building network 1\n');
	net.padded = blockwiseModules(DATA.padded, networkType=networkType, power=power, ...);
	mod.padded = exn.buildModuleGeneInfo(DATA.padded, net.padded$colors, annos, IDs);
	tr.padded = exn.transcriptsAcrossModules(names(DATA.padded), net.padded$colors);
	kME.padded = exn.computekME(DATA.padded, net.padded$MEs);
	padded = list(net=net.padded, mod=mod.padded, tr=tr.padded, kME=kME.padded);
	
	cat('Building network 2\n');
	net.alleles = blockwiseModules(DATA.alleles, networkType=networkType, power=power, ...);
	if (matchLabels) {
		net.alleles$colors = matchLabels(net.alleles$colors, net.padded$colors)[, 1];
		net.alleles$MEs = moduleEigengenes(DATA.alleles, net.alleles$colors)$eigengenes;
	}
	mod.alleles = exn.buildModuleGeneInfo(DATA.alleles, net.alleles$colors, annos, IDs);
	tr.alleles = exn.transcriptsAcrossModules(names(DATA.alleles), net.alleles$colors);
	kME.alleles = exn.computekME(DATA.alleles, net.alleles$MEs);
	alleles = list(net=net.alleles, mod=mod.alleles, tr=tr.alleles, kME=kME.alleles);
	
	return(list(padded=padded, alleles=alleles, all=all, ase=allele_names));
}

########################################################################
#####
##### DAVID 
#####
########################################################################
#############   		initialize and connect a DAVIDWebService object, then add network background list and gene lists from each module
# Depends on:			RDAVIDWebService::DAVIDWebService$new
#						RDAVIDWebService::addList()
#						RDAVIDWebService::setCurrentBackgroundPosition()
#						RDAVIDWebService::getBackgroundListNames()
#
## Arguments
### login:				character string DAVID Webservice login
### moduleGeneInfoList:	list structured as output from exn.buildModuleGeneInfo()
### idType:				character string naming type of gene identifier
### BG:					character vector holding unique gene identifiers, e.g. from exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork$allGeneIDs
### BGname:				character string naming background list
### setBG:				boolean indicating whether to set DAVIDWebService object's current background to BG 
#
## Value
### david:				connected DAVIDWebService object holding network background list and module gene lists

exn.startAndUploadBgAndModListsToDAVID = function(login='ahilliar@stanford.edu', moduleGeneInfoList=mod, idType='ENTREZ_GENE_ID', BG=all$geneIDs, BGname='allBG', setBG=T, verbose=T) {
	if (any(names(moduleGeneInfoList) != c('genes','geneAnnos','geneSyms','geneIDs'))) {
		stop('moduleGeneInfoList should be output from exn.buildModuleGeneInfo()');
	}
	# initialize DAVID object
	david = DAVIDWebService$new(email=login);
	# upload background
	cat('Uploading background... ', BGname, '\n', sep='');
	addList(david, inputIds=BG, idType=idType, listName=BGname, listType='Background');
	# add list for each module
	cat('Uploading modules... ', sep='');
	modNames = names(moduleGeneInfoList$geneIDs);
	for (m in 1:length(modNames)) {
		if (verbose) {cat(modNames[m], ' ', sep='')};
		cleanIDs = moduleGeneInfoList$geneIDs[[m]][!is.na(moduleGeneInfoList$geneIDs[[m]])];
		addList(david, inputIds=cleanIDs, idType=idType, listName=modNames[m], listType='Gene');
	}
	if (setBG) {
		cat('\nSetting background to ', BGname, '\n', sep='');
		setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == BGname));
	}
	return(david);
}

########################################################################

exn.startAndUploadBgAndModListsToDAVID2 = function (colors, login, idType, ids) {
	mods = names(table(colors));
	david = DAVIDWebService$new(email=login);
	addList(david, inputIds=ids, idType=idType, listName='BG', listType='Background');
	for (m in 1:length(mods)) {
		print(mods[m]);
		modGenes = ids[colors==mods[m]];
		addList(david, inputIds=modGenes, idType=idType, listName=mods[m], listType='Gene');
	}
	setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
	setCurrentGeneListPosition(david, 1);
	return(david);
}

.startAndUploadBgAndModGeneListToDAVID = function (modGenes, login, idType) {
	david = DAVIDWebService$new(email=login);
	addList(david, inputIds=unlist(modGenes), idType=idType, listName='BG', listType='Background');
	for (m in 1:length(modGenes)) {
		cat(m, ':', names(modGenes)[m], ' ', sep='');
		addList(david, inputIds=modGenes[[m]], idType=idType, listName=names(modGenes)[m], listType='Gene');
	}
	setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
	setCurrentGeneListPosition(david, 1);
	return(david);
}

########################################################################
#############   	get functional annotation charts from DAVID for each module in the network
# Depends on:		RDAVIDWebService::is.connected()
#					RDAVIDWebService::setCurrentBackgroundPosition()
#					RDAVIDWebService::getBackgroundListNames()
#					RDAVIDWebService::getCurrentBackgroundListPosition()
#					RDAVIDWebService::getGeneListNames()
#					RDAVIDWebService::setCurrentGeneListPosition()
#					RDAVIDWebService::getCurrentGeneListPosition()
#					RDAVIDWebService::getFunctionalAnnotationChart()
#
## Arguments
### david:			DAVIDWebService object, e.g. from output of exn.startAndUploadBgAndModListsToDAVID()
### setBgPosition:	optional number indicating desired background list position in david
#
## Value
### modDAVID:		list where each element is a data frame holding a given module's DAVID functional annotation chart

exn.getModChartsFromDAVID = function(david, setBgPosition=NULL, factorsToChars=T, verbose=T) {
	if (!is.connected(david)) {
		stop('DAVID object not connected to web server');
	}
	if (is.numeric(setBgPosition)) {
		setCurrentBackgroundPosition(david, setBgPosition);
	}
	cat('Using ', getBackgroundListNames(david)[getCurrentBackgroundListPosition(david)], ' as background\n', sep='');
	modDAVID = list();
	geneListNames = getGeneListNames(david);
	cat('Getting charts for module lists... ', sep='');
	for (m in 1:length(geneListNames)) {
		setCurrentGeneListPosition(david, m);
		if (verbose) {cat(getCurrentGeneListPosition(david), ':', geneListNames[getCurrentGeneListPosition(david)], ' ', sep='')};
		modDAVID[[m]] = getFunctionalAnnotationChart(david);
		if (factorsToChars) {
			checkFactor = which(lapply(modDAVID[[m]], class)=='factor');
			if (length(checkFactor)>0) {
				modDAVID[[m]][, checkFactor] = as.character(modDAVID[[m]][, checkFactor]);
			}
		}
		names(modDAVID)[m] = geneListNames[m];
	}
	# check for modules with 0 hits and fix colnames
	nrows = unlist(lapply(modDAVID, nrow));
	goodColnames = names(modDAVID[nrows>0][[1]]);
	if (goodColnames[4]=='X.') {
		goodColnames[4] = 'Percent';
	}
	for (m in 1:length(modDAVID)) {
		if (nrows[m]==0) {
			modDAVID[[m]] = as.data.frame(matrix(nrow=1, ncol=length(goodColnames)));
			names(modDAVID[[m]]) = goodColnames;
		} else {
			names(modDAVID[[m]]) = goodColnames;
		}
	}
	return(modDAVID);
}

########################################################################

.getGenesAvgkMEInMod = function(genes, modkME) {
	return(mean(modkME[rownames(modkME) %in% genes, 1]));
}
.addkMEColToDAVID = function(modDAVIDmod, modkME) {
	modDAVIDmod = cbind(modDAVIDmod, avgkME=1:nrow(modDAVIDmod));
	for (r in 1:nrow(modDAVIDmod)) {
		modDAVIDmod[r, ]$avgkME = .getGenesAvgkMEInMod(.getTermGenes(modDAVIDmod,r), modkME);
	}
	return(modDAVIDmod);
}
.addkMEColsToDAVIDList = function(modDAVID, modkME, order=T) {
	for (m in 1:length(modDAVID)) {
		if (nrow(modDAVID[[m]])==0 | names(modDAVID)[m]=='grey') {
			modDAVID[[m]] = cbind(modDAVID[[m]], avgkME=rep(NA,nrow(modDAVID[[m]])));
		} else {
			modDAVID[[m]] = .addkMEColToDAVID(modDAVID[[m]], modkME[[m]]);
			if (order) {
				modDAVID[[m]] = modDAVID[[m]][order(modDAVID[[m]]$avgkME, decreasing=T), ];
			}
		}
	}
	return(modDAVID);
}

.addkMEColsToDAVIDList2 = function(modDAVID, modkME, order=T) {
	for (m in 1:length(modDAVID)) {
		if (nrow(modDAVID[[m]])==0) {
			modDAVID[[m]] = modDAVID[[m]];
		} else {
			modDAVID[[m]] = .addkMEColToDAVID(modDAVID[[m]], modkME);
			if (order) {
				modDAVID[[m]] = modDAVID[[m]][order(modDAVID[[m]]$avgkME, decreasing=T), ];
			}
		}
	}
	return(modDAVID);
}

###############

.getGenesAvgGS = function (geneVec, GS, GScols=NULL) {
	if (is.null(GScols)) {
		GScols = grep('^GS', names(GS));
	} else if (is.numeric(GScols)) {
		GScols = GScols;
	} else {
		stop('Invalid GScols');
	}
	geneGS = GS[rownames(GS) %in% geneVec, GScols];
	return(apply(geneGS, 2, mean));
}

.addGSColsToDAVID = function(modDAVIDmod, GS, GScols=NULL) {
	if (is.null(GScols)) {
		GScols = grep('^GS', names(GS));
	} else if (is.numeric(GScols)) {
		GScols = GScols;
	} else {
		stop('Invalid GScols');
	};
	modDAVIDmod = as.data.frame(cbind(modDAVIDmod, matrix(nrow=nrow(modDAVIDmod),ncol=length(GScols))));
	names(modDAVIDmod)[(ncol(modDAVIDmod)-length(GScols)+1):ncol(modDAVIDmod)] = names(GS)[GScols];
	for (r in 1:nrow(modDAVIDmod)) {
		modDAVIDmod[r, (ncol(modDAVIDmod)-length(GScols)+1):ncol(modDAVIDmod)] = .getGenesAvgGS(.getTermGenes(modDAVIDmod,r), GS=GS, GScols=GScols);
	}
	return(modDAVIDmod);
}

.addGSColsToDAVIDList = function (modDAVID, GS, GScols=NULL) {
	newlist = modDAVID;
	for (m in 1:length(newlist)) {
		if (nrow(newlist[[m]]) == 0) {
			warning(paste('Skipping ', names(newlist)[m], ' because it is empty', sep=''));
			next;
		}
		newlist[[m]] = .addGSColsToDAVID(modDAVIDmod=newlist[[m]], GS=GS, GScols=GScols);
	}
	return(newlist);
}


###############
.convertTermGenes = function (modDAVIDmod,term,lookup,fromCol,toCol,lengthWarning=F) {
	genes = .getTermGenes(modDAVIDmod,term);
	newgenes = lookup[which(lookup[,fromCol] %in% genes), toCol];
	if (lengthWarning & length(genes) != length(newgenes)) {
		warning('Input and output gene lists are different lengths');
	}
	return(unique(newgenes));
}

.convertTermGenesForDAVIDChart = function (modDAVIDmod,lookup,fromCol,toCol,lengthWarning=F) {
	newchart = modDAVIDmod;
	for (r in 1:nrow(newchart)) {
		newgenes = .convertTermGenes(modDAVIDmod=newchart,term=r,lookup=lookup,fromCol=fromCol,toCol=toCol,lengthWarning=lengthWarning);
		newchart$Genes[r] = paste(newgenes, collapse=', ', sep='');
	}
	return(newchart);
}

.convertTermGenesForDAVIDChartList = function (modDAVID,lookup,fromCol,toCol,lengthWarning=F) {
	newlist = modDAVID;
	for (m in 1:length(newlist)) {
		if (nrow(newlist[[m]]) == 0) {
			warning(paste('Skipping ', names(newlist)[m], ' because it is empty', sep=''));
			next;
		}
		newlist[[m]] = .convertTermGenesForDAVIDChart(modDAVIDmod=newlist[[m]],lookup=lookup,fromCol=fromCol,toCol=toCol,lengthWarning=lengthWarning);
	}
	return(newlist);
}



.getkMEsFromAssignedModules = function (modkMEList, matchToData=T, data) {
	mat = modkMEList[[1]];
	names(mat) = c('kME','p.kME');
	for (m in 2:length(modkMEList)) {
		names(modkMEList[[m]]) = c('kME','p.kME');
		mat = rbind(mat, modkMEList[[m]])
	}
	if (matchToData) {
		mat = mat[match(names(data), rownames(mat)), ];
	}
	return(mat);
}

.getTermGenes = function(modDAVIDmod,term) {
	if (is.numeric(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[term], ', ')));
	} else if (is.character(term)) {
		return(unlist(strsplit(modDAVIDmod$Genes[modDAVIDmod$Term==term], ', ')));
	} else {
		stop()
	}
}
# assumes nums has gene names
.getAvgNumForGenes = function (genes, nums) {
	return(mean(nums[names(nums) %in% genes], na.rm=T));
}

.getSomeNumForGenes = function (genes, nums, f='mean') {
	n = nums[names(nums) %in% genes];
	if (all(is.na(n))) {
		return(NA);
	} else {
		return(eval(call(f, n, na.rm=T)));
	}
}

.getAvgNumForModules = function (modGenes, nums) {
	modNum = c();
	for ( m in modGenes ) {
		modNum = c(modNum, .getAvgNumForGenes(m, nums))
	}
	names(modNum) = names(modGenes);
	return(modNum)
}

.getSomeNumForGenesList = function (genesList, nums, f='mean') {
	vec = c();
	for (g in genesList) {
		vec = c(vec, .getSomeNumForGenes(genes=g, nums=nums, f=f));
	}
	names(vec) = names(genesList);
	return(vec);
}

.addGeneAvgColToDAVID = function(modDAVIDmod, nums, colname='colname') {
	modDAVIDmod = cbind(modDAVIDmod, colname=1:nrow(modDAVIDmod));
	for (row in 1:nrow(modDAVIDmod)) {
		g = .getTermGenes(modDAVIDmod, row);
		modDAVIDmod[row, ncol(modDAVIDmod)] = .getAvgNumForGenes(g, nums);
	}
	names(modDAVIDmod)[ncol(modDAVIDmod)] = colname;
	return(modDAVIDmod)
}

.addGeneAvgColToDAVIDList = function(modDAVID, nums, modGenesList, colname='colname') {
	newDAVID = modDAVID;
	for (m in 1:length(newDAVID)) {
		#mnums = nums[names(nums) %in% modGenesList[names(modGenesList)==names(newDAVID)[m]]]
		newDAVID[[m]] = .addGeneAvgColToDAVID(newDAVID[[m]], nums, deparse(substitute(nums)));
	}
	return(newDAVID);
}

.getModGenesRankedBykME = function(module_names,colors,kME) {
	outList = list();
	for (m in 1:length(module_names)) {
		outList[[m]] = exn.getModulekME(module_names[m],colors,kME)
	}
	names(outList) = module_names;
	return(outList);
}


.matColToNamedVec = function (mat, col) {
	if (is.character(col)) {
		col = match(col, colnames(mat));
	} else if (is.numeric(col)) {
		col = col;
	} else {
		stop('need to specify column with num or name');
	}
	vec = mat[, col];
	names(vec) = rownames(mat);
	return(vec);
}



.checkTermGeneEnrichment = function(modDAVIDmod, row, genes, bg=NULL) {
	termGenes = unlist(strsplit(modDAVIDmod[row,]$Genes, ', '));
	if (is.null(bg)) {
		bg = unique(unlist(strsplit(modDAVIDmod$Genes, ', ')));
	} else if (is.character(bg) & all(termGenes %in% bg)) {
		bg = bg;
	} else {
		warning('Using all term genes in module as background because bg is invalid');
		bg = unique(unlist(strsplit(modDAVIDmod$Genes, ', ')));
	}
	return(checkGeneListEnrichment(genes, termGenes, bg));
}

.checkTermGeneEnrichmentModule = function(modDAVIDmod, genes, bg=NULL, thresh=NULL) {
	outlist = list();
	for (row in 1:nrow(modDAVIDmod)) {
		tmp = .checkTermGeneEnrichment(modDAVIDmod=modDAVIDmod, row=row, genes=genes, bg=bg);
		if (is.numeric(thresh)) {
			if (tmp[[2]]$p.value <= thresh) {
				outlist[[length(outlist)+1]] = tmp;
				names(outlist)[length(outlist)] = modDAVIDmod$Term[row];
			}
		} else {
			outlist[[row]] = tmp;
			names(outlist)[row] = modDAVIDmod$Term[row];
		}
	}
	return(outlist);
}

.checkTermGeneEnrichmentAllModules = function(modDAVID, genes, modGenes, bg=NULL, thresh=NULL) {
	if (any(names(modDAVID)!=names(modGenes))) {
		stop('modDAVID and modGenes must have same modules in same order');
	}
	outlist = list();
	for ( mod in 1:length(modDAVID) ) {
		#print(names(modDAVID)[mod])
		if (sum(genes %in% modGenes[[mod]]) == 0) { outlist[[mod]] = 'None in module'; next };
		if (nrow(modDAVID[[mod]]) == 0 | all(is.na(modDAVID[[mod]]))) { outlist[[mod]] = 'No terms'; next };
		outlist[[mod]] = .checkTermGeneEnrichmentModule(modDAVID[[mod]], genes=genes, bg=modGenes[[mod]], thresh=thresh);
	}
	names(outlist) = names(modDAVID);
	return(outlist);
}


.printTermGeneEnrichment = function(modDAVIDmod, genes, thresh=.05) {
	for(row in 1:nrow(modDAVIDmod)){  
		tmp = .checkTermGeneEnrichment(modDAVIDmod,row,genes); 
		if(tmp[[2]]$p.value<thresh){
			print(modDAVIDmod[row,c(1,2,3,4,5,11:ncol(modDAVIDmod))]);
			print(tmp)
		}   
	}
}

.countTermsWithPatternInModules = function (pattern, modDAVID) {
	terms = exn.modDAVIDByTerm(modDAVID);
	termsU = unlist(terms$uniqueToMod);
	all = sort(table(unlist(terms[[1]][grepl(pattern, names(terms[[1]]))])));
	uni = sort(table(termsU[grepl(pattern, names(termsU))]));
	return(list(all=all,uniqueToMod=uni));
}

# nums must be named vec
.comparekMEWithNumsForFilteredGenesInModTerm = function (modDAVIDmod, term, nums, modkME, filter, ...) {
	y=gsub('()','',match.call()[4]); 
	g = .getTermGenes(modDAVIDmod, term);
	g.kME = modkME[rownames(modkME) %in% g, ];
	g.num = nums[names(nums) %in% g];
	tmp = cbind(g.kME, nums=g.num[match(rownames(g.kME), names(g.num))]);
	f = rownames(tmp) %in% filter;
	fcol=as.numeric(f);
	fcol[fcol==0] = 'grey'; fcol[fcol==1] = 'black';
	par(mfrow=c(1,3),oma=c(0,0,2,0));
	verboseScatterplot(tmp[,1], tmp[,3], abline=T, xlab='kME', ylab=y, frame.plot=F, col=fcol, pch=19, cex=2);
	suppressWarnings(verboseBoxplot(tmp[,1], as.factor(f), xlab='', ylab='kME', frame.plot=F, col='grey', ...));
	suppressWarnings(verboseBoxplot(tmp[,3], as.factor(f), xlab='', ylab=y, frame.plot=F, col='grey', ...));
	toTitle = paste(gsub('()','',match.call()[2]), ', ', gsub('()','',match.call()[3]), ', ', gsub('()','',match.call()[6]), sep='');
	title(toTitle, outer=T);
}


.plotTermCors = function(modDAVIDmod,cols=c(14,15), main='',bg='black',col='black', ...) {
	okrows1 = !is.nan(modDAVIDmod[, cols[1]]);
	okrows2 = !is.nan(modDAVIDmod[, cols[2]]);
	verboseScatterplot(modDAVIDmod[okrows1 & okrows2, cols[1]], modDAVIDmod[okrows1 & okrows2, cols[2]],abline=T,xlab=colnames(modDAVIDmod)[cols[1]],ylab=colnames(modDAVIDmod)[cols[2]],frame.plot=F,abline.col='red',bg=bg,pch=21,main=main,col=col, ...)
}

.plotTermCorsAllModulesThresh = function (modDAVID, cols=c(14,15), thresh=10, mfrow, bg='', col='', verbose=T, horizontal=NULL, vertical=NULL, ...) {
	underThresh = sapply(modDAVID, nrow) < thresh;
	tmp = modDAVID[!underThresh];
	xmax = max(sapply(tmp, function(f) max(as.numeric(f[,cols[1]]))), na.rm=T);
	xmin = min(sapply(tmp, function(f) min(as.numeric(f[,cols[1]]))), na.rm=T);
	ymax = max(sapply(tmp, function(f) max(as.numeric(f[,cols[2]]))), na.rm=T);
	ymin = min(sapply(tmp, function(f) min(as.numeric(f[,cols[2]]))), na.rm=T);
	cors = rep(0, length(modDAVID));
	names(cors) = names(modDAVID);
	par(mfrow=mfrow, oma=c(0,0,2,0));
	for (m in 1:length(modDAVID)) {
		if (nrow(modDAVID[[m]]) < thresh) {
			if (verbose) {
				cat('skipping the ', names(modDAVID)[m], ' module because it has <', thresh, ' enriched terms\n', sep='');
			}
			cors[m] = NA;
			next;
		}
		.plotTermCors(modDAVID[[m]], cols=cols, main=paste(names(modDAVID)[m], '\n', sep=''), bg=names(modDAVID)[m], col='black', xlim=c(xmin,xmax), ylim=c(ymin,ymax), ...);
		cors[m] = cor.test(modDAVID[[m]][,cols[1]], modDAVID[[m]][,cols[2]])$estimate;
		if (is.numeric(horizontal)) {abline(h=horizontal, col='darkgrey', lty='dashed')}
		if (is.numeric(vertical)) {abline(v=vertical, col='darkgrey', lty='dashed')}
	}
	title(paste(deparse(substitute(modDAVID)), ', thresh=', thresh, ' terms', sep=''), outer=T);
	return(cors);
}

# num1 and num2 must be vectors ordered the same as data and colors
.verboseScatterplotAcrossModules = function(num1, num2, colors, mfrow, xlab=NULL, ylab=NULL, colors.plot=NULL, horiz=NULL, vert=NULL, ...) {
	if (is.null(xlab)) {
		xlab=deparse(substitute(num1))
	} else {
		xlab=xlab;
	}
	if (is.null(ylab)) {
		ylab=deparse(substitute(num2))
	} else {
		ylab=ylab
	}
	par(mfrow=mfrow, ...);
	mods = names(table(colors));
	for (m in 1:length(mods)) {
		modgenes = colors==mods[m];
		if (is.null(colors.plot)) {
			col = mods[m];
		} else {
			col = colors.plot[modgenes];
		}
		verboseScatterplot(num1[modgenes], num2[modgenes],
						   abline=T, frame.plot=F,
						   col=col, main=mods[m],
						   pch=19,
						   xlab=xlab,
						   ylab=ylab, ...
						   );
		if (is.numeric(horiz)) {
			abline(h=horiz);
		}
		if (is.numeric(vert)) {
			abline(v=vert);
		}
	}
	
}

.convertmodDAVIDlistToDataFrame = function (modDAVID) {
	df = as.data.frame(cbind(modDAVID[[1]], module=rep(names(modDAVID)[1], nrow(modDAVID[[1]]))));
	for (m in 2:length(modDAVID)) {
		df = as.data.frame(rbind(df, cbind(modDAVID[[m]], module=rep(names(modDAVID)[m], nrow(modDAVID[[m]])))))
	}
	return(df)
}

# nums must be ordered same as colors
.getNumsForModule = function (nums, module, colors, normalize=T, ...) {
	mnums = nums[colors==module];
	if (normalize) {
		mnums = mnums/max(mnums, na.rm=T);
	}
	return(mnums);
}

.getNumsForModulesList = function (nums, colors, normalize=T, ...) {
	modules = names(table(colors));
	modNums = list();
	for (m in 1:length(modules)) {
		modNums[[m]] = .getNumsForModule(nums=nums, module=modules[m], colors=colors, normalize=normalize, ...);
		names(modNums)[m] = modules[m];
	}
	return(modNums);
}


.vecListToVec = function (vecList, namesToMatch) {
	newVec = c();
	for (i in vecList) {
		newVec = c(newVec, i);
	}
	newVec = newVec[match(namesToMatch, names(newVec))];
	names(newVec) = namesToMatch;
	return(newVec);
}

# rows are samples, columns are genes
.permuteSamples = function (data1, data2) {
	data = rbind(data1, data2);
	numTotal = nrow(data);
	ind1 = sample(1:numTotal, nrow(data1), replace=F);
	ind2 = setdiff(1:numTotal, ind1);
	out = list(data[ind1, ], data[ind2, ]);
	names(out) = c(paste('p', deparse(substitute(data1)), sep=''), 
				   paste('p', deparse(substitute(data2)), sep=''));
	return(out);
}

.kINdiffPermutationTest = function (data1, data2, computekINdiffOutput, permutations=3, networkType='signed', power=18, ...) {
	if (any(names(data1) != names(data2))) {
		stop('datasets must contain same genes in same order');
	}
	check = sum(grepl('module', names(computekINdiffOutput)));
	if (check==1) {
		colors2 = NULL;
	} else if (check==2) {
		colors2 = computekINdiffOutput$module2;
	} else {
		stop('must use output from .computekINdiff');
	}

	diffNull = matrix(nrow=ncol(data1), ncol=permutations);
	for (p in 1:permutations) {
		# # # permutedData = .permuteSamples(data1, data2);
		# # # tmp = .computekINdiff(dataRef=permutedData[[1]], data2=permutedData[[2]], 
							  # # # colorsRef=computekINdiffOutput$moduleRef, colors2=colors2, 
							  # # # type=networkType, power=power, 
							  # # # ...);
		
	}
}

# tmp = table(unlist(terms.hsf$uniqueToMod));
# tmp = tmp[match(names(table(colors.hs)[-14]), names(tmp))];
# names(tmp) = names(table(colors.hs)[-14]);
# tmp[is.na(tmp)] = 0;


# par(mfrow=c(1,4));
# verboseScatterplot(table(colors.hs)[-14], sapply(modDAVID.hs,nrow),
				   # xlab='module size (genes)', ylab='number of enriched terms',
				   # abline=T, abline.lty='dashed', frame.plot=F,
				   # pch=19, col=names(table(colors.hs)[-14]), cex=1.8,
				   # ylim=c(0,200)
				   # );
# verboseScatterplot(table(colors.hs)[-14], sapply(modDAVID.hsf,nrow),
				   # xlab='module size (genes)', ylab='number of enriched terms with FDR<10',
				   # abline=T, abline.lty='dashed', frame.plot=F,
				   # pch=19, col=names(table(colors.hs)[-14]), cex=1.8,
				   # ylim=c(0,200)
				   # );
# verboseScatterplot(table(colors.hs)[-14], tmp,
				   # xlab='module size (genes)', ylab='number of unique enriched terms with FDR<10',
				   # abline=T, abline.lty='dashed', frame.plot=F,
				   # pch=19, col=names(table(colors.hs)[-14]), cex=1.8,
				   # ylim=c(0,200)
				   # );
# verboseScatterplot(table(colors.hs)[-14], tmp / sapply(modDAVID.hs,nrow),
				   # xlab='module size (genes)', ylab='% of terms both unique and have FDR<10',
				   # abline=T, abline.lty='dashed', frame.plot=F,
				   # pch=19, col=names(table(colors.hs)[-14]), cex=1.8,
				   # ylim=c(0,1)
				   # );
########################################################################

.getTermClustersFromDAVID = function (david, overlap=5L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=100L, EASEthresh=NULL) {
	cl = getClusterReport(david, overlap=overlap, initialSeed=initialSeed, finalSeed=finalSeed, linkage=linkage, kappa=kappa);
	cl = attributes(cl)$cluster;
	if (is.numeric(EASEthresh)) {
		cl = cl[sapply(cl, function(f) f$EnrichmentScore) >= EASEthresh];
	}
	return(cl);
}

.getModTermClustersFromDAVID = function (david, setBgPosition=NULL, overlap=5L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=100L, EASEthresh=NULL, verbose=T) {
	if (!is.connected(david)) {
		stop('DAVID object not connected to web server');
	}
	if (is.numeric(setBgPosition)) {
		setCurrentBackgroundPosition(david, setBgPosition);
	}
	cat('Using ', getBackgroundListNames(david)[getCurrentBackgroundListPosition(david)], ' as background\n', sep='');
	modDAVID = list();
	geneListNames = getGeneListNames(david);
	cat('Getting charts for module lists... ', sep='');
	
	for (m in 1:length(geneListNames)) {
		setCurrentGeneListPosition(david, m);
		if (verbose) {cat(getCurrentGeneListPosition(david), ':', geneListNames[getCurrentGeneListPosition(david)], ' ', sep='')};
		modDAVID[[m]] = .getTermClustersFromDAVID(david, overlap=overlap, initialSeed=initialSeed, finalSeed=finalSeed, linkage=linkage, kappa=kappa, EASEthresh=EASEthresh);
		names(modDAVID)[m] = geneListNames[m];
	}
	l = unlist(lapply(modDAVID, length)); 
	if (any(l == 0)) {
		warning('Check ', '  ', names(modDAVID)[l == 0], ' ', sep='');
	}
	return(modDAVID);
}

## assumes all terms within a cluster have the xact same genes, i.e. that .getTermClustersFromDAVID was ran with kappa=100
.collapseGenesFromDAVIDClusters = function (cl, EASEthresh=NULL) {
	if (is.numeric(EASEthresh)) {
		cl = cl[sapply(cl, function(f) f$EnrichmentScore) >= EASEthresh];
	}
	out = list();
	for (i in 1:length(cl)) {
		this = cl[[i]];
		out[[i]] = list(EASE=this$EnrichmentScore, Terms=this$Members$Term, Genes=.getTermGenes(this$Members, this$Members$Term[1]));
	}
	return(out);
}

.removeOverlapsInGeneListPermutations = function (listOfGeneVecs) {
	numLists = length(listOfGeneVecs);
	for (gl in 1:numLists) {
		incommon = c();
		for (ol in which(1:numLists != gl)) {
			incommon = c(incommon, intersect(listOfGeneVecs[[gl]], listOfGeneVecs[[ol]]));
		}
		listOfGeneVecs[[gl]] = setdiff(listOfGeneVecs[[gl]], incommon);
	}
	return(listOfGeneVecs);
}



########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getChartFromDAVID = function(geneList, listName='list', BG, BGname='BG', login='ahilliar@stanford.edu', idType='ENTREZ_GENE_ID', setBG=T, suppress=F) {
	d = DAVIDWebService$new(email=login);
	addList(d, inputIds=BG, idType=idType, listName=BGname, listType='Background');
	addList(d, inputIds=geneList, idType=idType, listName=listName, listType='Gene');
	setCurrentBackgroundPosition(d, which(getBackgroundListNames(d) == BGname));
	chart = getFunctionalAnnotationChart(d);
	if (!suppress) {
		print(d);
	}
	return(list(david=d, chart=chart));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getModuleClustCoef = function(module, colors, DATA, adjMat=NULL, networkType='signed', power=18, suppress=F) {
	modGenes = colors==module;
	if (is.null(adjMat)) {
		if (!suppress) {cat('...using adjacency(DATA[, colors==module]) to compute module adjacency matrix\n', sep='')};
		adjMat = adjacency(DATA[, modGenes], type=networkType, power=power);
	}
	ccoef = clusterCoef(adjMat);
	return(ccoef);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.getModuleClustCoefList = function(colors, DATA, adjMat=NULL, networkType='signed', power=18, suppress=T) {
	modClustCoefs = list();
	for (m in 1:length(table(colors))) {
		modClustCoefs[[m]] = exn.getModuleClustCoef(module=names(table(colors))[m], colors=colors, DATA=DATA, adjMat=adjMat, networkType=networkType, power=power, suppress=suppress);
		names(modClustCoefs)[m] = names(table(colors))[m];
	}
	return(modClustCoefs);
}
 
########################################################################
#############   
# Depends on:	exn.viewMod()
#				exn.getModulekME()
# 				exn.EntrezIDsToGeneSymbols()
#				exn.geneParser()
#				exn.mRNAParser()
#
## Arguments
### chart:
### IDs:
#
## Value
###

exn.addGeneInfoToDAVIDChart = function(chart, IDs, moduleGeneInfoList, module, colors, DATA, networkType='signed', power=18, kMEtable=kMEs, kMEpvals=NULL, reorderColNum=NULL) {
	if (ncol(chart) != 13) {
		stop('chart must be output from RDAVIDWebService::getFunctionalAnnotationChart()');
	}
	numColToAdd = 6;
	newColNames = c('Gene.Syms', 'avg.kME', 'avg.kME.p', 'avg.kME.FDR', 'avg.kME.p.FDR', 'avg.clustCoef');
	geneCol = match('Genes', names(chart));
	if (nrow(chart)==1 & all(is.na(chart[1, ]))) {
		cat('...skipping since chart is empty, but adding dummy columns\n', sep='');
		chart_new = cbind(chart[, 1:(geneCol-1)], chart[, (geneCol+1):ncol(chart)], Genes=chart[, geneCol], matrix(nrow=nrow(chart), ncol=numColToAdd));
		names(chart_new)[14:(ncol(chart)+numColToAdd)] = newColNames;
		return(chart_new);
	}
	modInfo = exn.viewMod(moduleGeneInfoList, module);
	modkME = exn.getModulekME(module, colors, kMEtable, kMEpvals=kMEpvals);#print(kMEpvals)
	modClustCoef = exn.getModuleClustCoef(module=module, colors=colors, DATA=DATA, networkType=networkType, power=power);
	chart_new = cbind(chart[, 1:(geneCol-1)], chart[, (geneCol+1):ncol(chart)], Genes=chart[, geneCol], matrix(nrow=nrow(chart), ncol=numColToAdd));
	names(chart_new)[14:(ncol(chart)+numColToAdd)] = newColNames;
	termsIDs = strsplit(chart$Genes, ', ');
	for (term in 1:length(termsIDs)) {
		tmp = exn.EntrezIDsToGeneSymbols(termsIDs[[term]], IDlist=IDs);
		chart_new$Gene.Syms[term] = paste(tmp[match(termsIDs[[term]], names(tmp))], collapse=', ');
		genes = exn.geneParser(modInfo$geneAnnos[match(tmp, modInfo$geneSyms), 1]); #print(genes)
		transcriptkME = modkME[exn.mRNAParser(rownames(modkME)) %in% genes, 1];
		transcriptClustCoef = modClustCoef[exn.mRNAParser(names(modClustCoef)) %in% genes];#print(mean(transcriptClustCoef))
		chart_new$avg.kME[term] = mean(transcriptkME);
		chart_new$avg.kME.p[term] = mean(transcriptkME) * (1-chart_new$PValue[term]);
		chart_new$avg.kME.FDR[term] = mean(transcriptkME) * ((100-chart_new$FDR[term])/100);
		chart_new$avg.kME.p.FDR[term] = mean(transcriptkME) * (1-chart_new$PValue[term]) * ((100-chart_new$FDR[term])/100);
		chart_new$avg.clustCoef[term] = mean(transcriptClustCoef);
	}
	if (is.numeric(reorderColNum)) {
		chart_new = chart_new[order(chart_new[, reorderColNum], decreasing=T), ];
	}
	collectGarbage();
	return(chart_new);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.addGeneInfoToDAVIDChart2 = function(chart, IDs, modGeneInfoList, module, colors, DATA, GStable, modFNCList, kMEtable, type='signed', power=18) {
	if (ncol(chart) != 13) {
		stop('chart must be output from RDAVIDWebService::getFunctionalAnnotationChart()');
	}
	if (length(modGeneInfoList[[1]])!=length(modFNCList)) {
		stop();
	}
	if (any(names(modFNCList[[1]])!=c('Connectivity', 'ScaledConnectivity', 'ClusterCoef', 'MAR', 'Density', 'Centralization', 'Heterogeneity'))) {
		stop();
	}
	GSnames = names(GStable)[seq(1, ncol(GStable), 3)];
	newNames = c('Gene.Syms', 'kME', 'kIN', 'clustCoef', 'MAR', 'density', 'centralization', 'heterogeneity', GSnames);
	numColToAdd = length(newNames);
	geneCol = match('Genes', names(chart));
	if (nrow(chart)==1 & all(is.na(chart[1, ]))) {
		cat('...skipping since chart is empty, but adding dummy columns\n', sep='');
		chart_new = cbind(chart[, 1:(geneCol-1)], chart[, (geneCol+1):ncol(chart)], Genes=chart[, geneCol], matrix(nrow=nrow(chart), ncol=numColToAdd));
		names(chart_new)[14:(ncol(chart)+numColToAdd)] = newNames;
		return(chart_new);
	} else {
		modInfo = exn.viewMod(modGeneInfoList, module);
		modkME = exn.getModulekME(module, colors, kMEtable, kMEpvals=NULL);
		ind = match(module, names(modFNCList));
		modkIN = modFNCList[[ind]]$Connectivity;															#print(modkIN)
		modClustCoef = modFNCList[[ind]]$ClusterCoef;
		modMAR = modFNCList[[ind]]$MAR;
		if (all(names(modFNCList)==names(modGeneInfoList$genes))) {
			modDATA = DATA[, names(DATA) %in% modGeneInfoList$genes[[ind]]];
		} else {
			stop('cannot subset data with fnc ind');
		}
		chart_new = cbind(chart[, 1:(geneCol-1)], chart[, (geneCol+1):ncol(chart)], Genes=chart[, geneCol], matrix(nrow=nrow(chart), ncol=numColToAdd));
		names(chart_new)[14:(ncol(chart)+numColToAdd)] = newNames;
		termsIDs = strsplit(chart$Genes, ', ');
		for (term in 1:length(termsIDs)) {
			tmp = exn.EntrezIDsToGeneSymbols(termsIDs[[term]], IDlist=IDs);									#print(tmp)
			genes = exn.geneParser(modInfo$geneAnnos[match(tmp, modInfo$geneSyms), 1]);						#print(genes)
			chart_new$Gene.Syms[term] = paste(tmp[match(termsIDs[[term]], names(tmp))], collapse=', ');		
			chart_new$kME[term] = mean(modkME[exn.mRNAParser(rownames(modkME)) %in% genes, 1]);
			chart_new$kIN[term] = mean(modkIN[exn.mRNAParser(names(modkIN)) %in% genes]);
			chart_new$clustCoef[term] = mean(modClustCoef[exn.mRNAParser(names(modClustCoef)) %in% genes]);
			chart_new$MAR[term] = mean(modMAR[exn.mRNAParser(names(modMAR)) %in% genes]);
			if (length(genes)<3) {
				chart_new$density[term] = NA;
				chart_new$centralization[term] = NA;
				chart_new$heterogeneity[term] = NA;
			} else {
				transcriptFNC = fundamentalNetworkConcepts(adjacency(modDATA[, exn.mRNAParser(names(modDATA)) %in% genes], type=type, power=power));#print(transcriptFNC)
				chart_new$density[term] = transcriptFNC$Density;
				chart_new$centralization[term] = transcriptFNC$Centralization;
				chart_new$heterogeneity[term] = transcriptFNC$Heterogeneity;
			}
			for (gs in 1:length(GSnames)) {
				
			}
		}
	}
	return(chart_new);
}


########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.addGeneInfoToDAVIDChartList = function(modDAVID, IDs, moduleGeneInfoList, colors, kMEtable, DATA, kMEpvals=NULL, reorderColNum=NULL, verbose=T) {
	exn.checkNamesOfModListElements(modDAVID);
	modDAVID_new = list();
	if (verbose) {cat('working on:\n')};
	for (m in 1:length(modDAVID)) {
		#cat('  ', names(modDAVID)[m], '          ', '\r', sep='');flush.console()
		if (verbose) {cat(names(modDAVID)[m], ' ', sep='')};
		modDAVID_new[[m]] = exn.addGeneInfoToDAVIDChart(chart=modDAVID[[m]], moduleGeneInfoList=moduleGeneInfoList, module=names(modDAVID)[m], colors=colors, kMEtable=kMEtable, DATA=DATA, kMEpvals=kMEpvals, reorderColNum=reorderColNum);
	}
	names(modDAVID_new) = names(modDAVID);
	return(modDAVID_new);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.modDAVIDByTerm = function(modDAVID) {
	exn.checkNamesOfModListElements(modDAVID);
	tmp = unlist(modDAVID);
	terms = tmp[grepl('Term', names(tmp))];
	terms = terms[!is.na(terms)];
	termList = list();
	for (tt in 1:length(unique(terms))) {
		termMods = terms[terms == unique(terms)[tt]];
		tmp2 = unlist(strsplit(names(termMods), '\\.'));
		mods = tmp2[!grepl('Term', tmp2)];
		termList[[tt]] = mods;
		names(termList)[tt] = unique(terms)[tt];
	}
	out = list(modsByTerm=termList, uniqueToMod=termList[unlist(lapply(termList, length)) == 1]);
	return(out);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
### 
### 
#
## Value
###

exn.getTermInfoAcrossMods = function(term, modDAVID, modDAVIDByTerm.out=NULL) {
	if (is.null(modDAVIDByTerm.out)) {
		cat('...using exn.modDAVIDByTerm() to get distribution of term "', term, '" across modules\n', sep='');
		modDAVIDByTerm.out = exn.modDAVIDByTerm(modDAVID);
	} else if (any(names(modDAVIDByTerm.out) != c('modsByTerm', 'uniqueToMod'))) {
		stop('modDAVIDByTerm.out must be output from exn.modDAVIDByTerm()');
	}
	termInd = match(term, names(modDAVIDByTerm.out$modsByTerm));
	termMods = modDAVIDByTerm.out$modsByTerm[[termInd]];
	tmpDAVID = modDAVID[match(termMods, names(modDAVID))];
	out = vector(mode='list', length=length(tmpDAVID));
	names(out) = names(tmpDAVID);
	for (m in 1:length(tmpDAVID)) {
		out[[m]] = tmpDAVID[[m]][tmpDAVID[[m]]$Term==term, ];
	}
	return(out);
}

exn.getAllTermsWithGene = function (gene, modDAVIDmod) {
	return(modDAVIDmod[apply(modDAVIDmod, 1, function(f) grepl(gene, f[match('Genes', names(f))])), ]);
}

exn.getAllTermsWithGenes = function (genes, modDAVIDmod) {
	if (sum(genes %in% unlist(strsplit(modDAVIDmod$Genes,', ')))==0) {
		return('No genes in any term');
	}
	mat = exn.getAllTermsWithGene(genes[1], modDAVIDmod);
	for (g in 2:length(genes)) {
		mat = rbind(mat, exn.getAllTermsWithGene(genes[g], modDAVIDmod));
	}
	mat = mat[!duplicated(mat), ];
	mat = cbind(mat, qgenes=rep(NA, nrow(mat)));
	orderBy = c();
	for (tt in 1:nrow(mat)) {
		tgenes = strsplit(mat[tt, match('Genes', colnames(mat))], ', ')[[1]];
		g = genes[genes %in% tgenes];
		mat$qgenes[tt] = paste(g, collapse=', ');
		orderBy = c(orderBy, length(g));
	}
	mat = mat[order(orderBy, decreasing=T), ];
	return(mat);
}

exn.getAllTermsWithGenesAcrossModules = function (genes, modDAVID) {
	outlist = list();
	for (mod in 1:length(modDAVID)) {
		#print(names(modDAVID)[mod]);
		if (nrow(modDAVID[[mod]])==0 | all(is.na(modDAVID[[mod]]))) { 
			#print('YYYYYYY')
			outlist[[mod]] = 'No terms'; 
		} else {
			#print('XXXXXXXX')
			outlist[[mod]] = exn.getAllTermsWithGenes(genes, modDAVID[[mod]]);
		}
		#print(outlist[[mod]])
	}
	names(outlist) = names(modDAVID);
	return(outlist);
}

exn.getAllTermsWithGenesFromTerm = function (term, modDAVIDmod) {
	tcol = which(colnames(modDAVIDmod)=='Term');
	gcol = which(colnames(modDAVIDmod)=='Genes');
	row = which(modDAVIDmod[, tcol]==term); 
	genes = strsplit(modDAVIDmod[row, gcol], ', ')[[1]]    
	return(exn.getAllTermsWithGenes(genes, modDAVIDmod));
}


exn.getDAVIDchartOverlapByTerm = function (chart1, chart2, wc_filter=T) {
	terms1 = chart1[, which(colnames(chart1)=='Term')];
	terms2 = chart2[, which(colnames(chart2)=='Term')];
	overlap = intersect(terms1, terms2);
	chart1only = terms1[!(terms1 %in% terms2)];
	chart2only = terms2[!(terms2 %in% terms1)];
	wc_overlap = sort(table(tolower(unlist(strsplit(overlap, '\\ |~')))), decreasing=T);
	wc_chart1only = sort(table(tolower(unlist(strsplit(chart1only, '\\ |~')))), decreasing=T);
	wc_chart2only = sort(table(tolower(unlist(strsplit(chart2only, '\\ |~')))), decreasing=T);
	if (wc_filter) {
		remove = c('to', 'of', 'and', 'or', 'in');
		wc_overlap = wc_overlap[!(names(wc_overlap) %in% remove)];
		wc_overlap = wc_overlap[wc_overlap>1];
		wc_chart1only = wc_chart1only[!(names(wc_chart1only) %in% remove)];
		wc_chart1only = wc_chart1only[wc_chart1only>1];
		wc_chart2only = wc_chart2only[!(names(wc_chart2only) %in% remove)];
		wc_chart2only = wc_chart2only[wc_chart2only>1];
	}
	wc = list(overlap=wc_overlap, only1=wc_chart1only, only2=wc_chart2only);
	return(list(common=overlap, only1=chart1only, only2=chart2only, wc=wc));
}



########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.checkNamesOfModListElements = function(modList) {
	x = lapply(modList, names);
	x_lengths = unlist(lapply(x, length));
	if ( (!(all(x_lengths[1]==x_lengths))) | length(unique(x_lengths))!=1 ) {
		stop('List elements have different numbers of columns');
	}
	for (element in 2:length(x)) {
		if (!(all(x[[1]]==x[[element]]))) {
			stop('Names in element ', element, ' do not match');
		}
	}
	return(TRUE);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.verboseScatterplotAcrossModDAVIDList = function(modDAVID, cols=c(15, 19), abline=T, abline.color='red', abline.lty='dashed', frame.plot=F, cex=3, pch=20, ...) {
	exn.checkNamesOfModListElements(modDAVID);
	nrows = unlist(lapply(modDAVID, nrow));
	goodChart = modDAVID[nrows>1][[1]];
	if (is.numeric(cols)) {
		if (any(!is.numeric(c(goodChart[, cols[1]], goodChart[, cols[2]])))) {
			stop('At least one of the columns specified does not contain numeric data');
		}
		col1 = cols[1];
		col2 = cols[2];
		xlab = names(goodChart)[col1];
		ylab = names(goodChart)[col2];
	} else if (is.character(cols)) {
		if (!(all(cols %in% names(goodChart)))) {
			stop('cols must match names in list elements');
		}
		col1 = match(cols[1], names(goodChart));
		col2 = match(cols[2], names(goodChart));
		xlab = cols[1];
		ylab = cols[2];
	} else {
		stop('cols must be a 2 element numeric or character vector');
	}
	x = unlist(lapply(modDAVID, function(chart){mean(chart[, col1])}));
	y = unlist(lapply(modDAVID, function(chart){mean(chart[, col2])}));
	verboseScatterplot(x, y, abline=abline, abline.color=abline.color, abline.lty=abline.lty, xlab=xlab, ylab=ylab, frame.plot=frame.plot, col=names(x), cex=cex, pch=pch, ...);
	return(list(x, y));
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.corDAVIDcolsInMods = function(modDAVID, mfrow, cols=c()) {
	exn.checkNamesOfModListElements(modDAVID);
	nrows = unlist(lapply(modDAVID, nrow));
	goodChart = modDAVID[nrows>1][[1]];
	if (is.numeric(cols)) {
		if (any(!is.numeric(c(goodChart[, cols[1]], goodChart[, cols[2]])))) {
			stop('At least one of the columns specified does not contain numeric data');
		}
		col1 = cols[1];
		col2 = cols[2];
		xlab = names(goodChart)[col1];
		ylab = names(goodChart)[col2];
	} else if (is.character(cols)) {
		if (!(all(cols %in% names(goodChart)))) {
			stop('cols must match names in list elements');
		}
		col1 = match(cols[1], names(goodChart));
		col2 = match(cols[2], names(goodChart));
		xlab = cols[1];
		ylab = cols[2];
	} else {
		stop('cols must be a 2 element numeric or character vector');
	}
	cors = c();
	for (m in 1:length(modDAVID)) {
		thisMod = modDAVID[[m]];
		if (nrow(thisMod)<3) {
			cors = c(cors, NA);
			names(cors)[m] = names(modDAVID)[m];
			cat(names(modDAVID)[m], ' only has ', nrow(thisMod), ' rows\n', sep='');
			next;
		}
		cors = c(cors, cor.test(thisMod[,col1], thisMod[,col2])$estimate);
		names(cors)[m] = names(modDAVID)[m];
	}
	return(cors);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.filterDAVIDChart = function(chart, column, value, lessThan=F, verbose=T) {
	if (!is.data.frame(chart)) {
		stop('Input must be a data frame');
	}
	if (is.character(column)) {
		columnName = column;
		column = match(columnName, names(chart));#print(column)
	}
	if (is.character(value)) {
		if (!(is.character(chart[, column]))) {
			stop('Cannot filter ', mode(chart[, column]), ' column by a character value');
		}
		if (verbose) {cat('Filtering $', names(chart)[column], ' to rows that include ', value, '\n', sep='')};
		chart_filtered = chart[grepl(value, chart[, column]), ];
	} else if (is.numeric(value)) {
		if (!(is.numeric(chart[, column]))) {
			stop('Cannot filter ', mode(chart[, column]), ' column by a numeric value');
		}
		if (lessThan) {
			if (verbose) {cat('Filtering $', names(chart)[column], ' to values < ', value, '\n', sep='')};
			chart_filtered = chart[chart[, column]<value, ];
		} else {
			if (verbose) {cat('Filtering $', names(chart)[column], ' to values > ', value, '\n', sep='')};
			chart_filtered = chart[chart[, column]>value, ];
		}
	} else {
		stop('Parameter "value" must be numeric or character');
	}
	return(chart_filtered);
}

########################################################################
#############   
# Depends on:	
#
## Arguments
###
###
#
## Value
###

exn.filterModDAVIDList = function(modDAVID, column, value, lessThan=F, verbose=T) {
	exn.checkNamesOfModListElements(modDAVID);
	if (is.character(column)) {
		columnName = column;
		column = match(columnName, names(modDAVID[[1]]));
	}
	nrows = unlist(lapply(modDAVID, nrow));
	goodColnames = names(modDAVID[nrows>1][[1]]);
	if (verbose) {
		if (is.character(value)) {
			cat('Filtering $', goodColnames[column], ' to rows that include ', value, '\n', sep='');
		} else if (is.numeric(value)) {
			if (lessThan) {
				cat('Filtering $', goodColnames[column], ' to values < ', value, '\n', sep='');
			} else {
				cat('Filtering $', goodColnames[column], ' to values > ', value, '\n', sep='');
			}
		} else {
			stop('Parameter "value" must be numeric or character');
		}
	} 
	modDAVID_filtered = modDAVID;
	for (m in 1:length(modDAVID_filtered)) {
		if (all(is.na(modDAVID_filtered[[m]]))) {
			modDAVID_filtered[[m]] = modDAVID_filtered[[m]];
			next;
		}
		modDAVID_filtered[[m]] = exn.filterDAVIDChart(chart=modDAVID_filtered[[m]], column=column, value=value, lessThan=lessThan, verbose=F);
	}
	return(modDAVID_filtered);
}

########################################################################
#####
##### compare networks 
#####
########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.plotModuleOverlaps = function(colors1, colors2, fcex=1.00, pcex=.7, fcex1=.7, pcex1=1.00, thresh=T, gamma=1, plot=T, ...) {
	overlap = overlapTable(colors1, colors2);
	
	numMat = -log10(overlap$pTable);
	numMat[numMat > 50] = 50;
	 
	textMat = paste(overlap$countTable, '\n', signif(overlap$pTable, 2));
	dim(textMat) = dim(numMat);
	 
	xlabels = paste('M', sort(unique(colors2)));
	ylabels = paste('M', sort(unique(colors1)));
	xSymbols = paste(sort(unique(colors2)), ': ', table(colors2), sep = '');
	ySymbols = paste(sort(unique(colors1)), ': ', table(colors1), sep = '');
	
	if (thresh) {
		thresh=.05/(ncol(textMat)*nrow(textMat));
		textMat[overlap$pTable>thresh] = ''
	}
	if (plot) {
		sizeGrWindow(7, 7); fp = FALSE
		par(mar = c(6, 7, 2, 1.0));
		labeledHeatmap(Matrix = numMat,
					   xLabels = xlabels, xSymbols = xSymbols,
					   yLabels = ylabels, ySymbols = ySymbols,
					   colorLabels = T, 
					   colors = blueWhiteRed(200,gamma=gamma)[80:200],
					   textMatrix = textMat, cex.text = pcex, setStdMargins = F,
					   cex.lab = fcex1,
					   xColorWidth = 0.08,
					   main = '', ...
					   );
	}
	return(overlap);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.modulePreservation = function(DATA1, DATA2, colors1, colors2, setLabels=c('DATA1', 'DATA2'), ranksOnly=T, nPermutations=NULL, networkType='signed', maxModuleSize=max(table(colors1)), saveStats=F, filePrefix=NULL, verbose=3, ...) {
	setLabels = setLabels;
	multiExpr = list();
	multiExpr[[1]] = list(data=DATA1);
	multiExpr[[2]] = list(data=DATA2);
	names(multiExpr) = setLabels;
	colorList = list();
	colorList[[1]] = colors1;
	colorList[[2]] = colors2;
	names(colorList) = setLabels;
	if (ranksOnly) {
		nPermutations = 0;
	} else {
		if (is.null(nPermutations)) {
			stop('if not computing preservation ranks only, i.e. ranksOnly=F, nPermutations must be specified');
		}
		nPermutations = nPermutations;
	}
	mp = modulePreservation(multiExpr, colorList, networkType=networkType, nPermutations=nPermutations, maxModuleSize=maxModuleSize, verbose=verbose, savePermutedStatistics=F, ...);
	if (saveStats) {
		save(mp, file=paste(filePrefix, '_modulePreservation.RData', sep=''));
	}
	return(mp);
}

########################################################################

exn.plotMETraitCorsAgainstNums = function (traitCors, modNums, mfrow, ylim, xlab=NULL, ...) {
	if (is.null(xlab)) {xlab = deparse(substitute(modNums))}
	par(mfrow=mfrow);
	for (tr in 1:ncol(traitCors)) {
		verboseScatterplot(modNums, 
						   traitCors[match(names(modNums), gsub('ME','',rownames(traitCors))), tr],
						   abline=T,
						   pch=19,
						   col=names(modNums),
						   frame.plot=F,
						   xlab=xlab,
						   ylab=colnames(traitCors)[tr],
						   cex=2,
						   ylim=ylim, ...
						   );
	}
}


########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.simulateReplicatesVec = function(DATAvec, n, disp, prob) {
	df = as.data.frame(matrix(nrow=, ncol=))
	for (i in 1:length(DATAvec)) {
		
	}
	return()
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.computeAndPlotMETraitCors = function(traits, MEs, col=blueWhiteRed(50), cex.text=.5, main='', thresh=T, ...) {
	traitCor = cor(MEs, traits, use='p');
	traitCorP = corPvalueStudent(traitCor, nrow(traits));
	textMat = paste(signif(traitCor, 2), '\n(', signif(traitCorP, 1), ')', sep='');
	dim(textMat) = dim(traitCor);
	if (thresh) {
		thresh=.05/(ncol(traits)*nrow(traits));
		textMat[traitCorP>thresh] = ''
	}
	par( mar = c(8, 9, 3, 3) );
	labeledHeatmap(Matrix=traitCor, xLabels=names(traits), yLabels=names(MEs), ySymbols=names(MEs), colorLabels=F, colors=col, textMatrix=textMat, setStdMargins=F, cex.text=cex.text, zlim=c(-1, 1), main=main, ...); 
	return(list(cor=traitCor, pval=traitCorP));
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
###

exn.plotDendroAndColors = function(dendro, colors, main='', rowText=T, groupLabels='module', addGuide=T, guideHang=.05, hang=.1, block=NULL, blockGenes=NULL, ...) {
	if (is.null(block)) {
		dendro = dendro;
		colors = colors;
	} else if (is.numeric(block)) {
		if (is.null(blockGenes)) {
			stop('blockGenes required');
		}
		dendro = dendro[[block]];
		colors = colors[blockGenes[[block]]];
	} else {
		stop('block must be NULL or numeric');
	}
	if (rowText) {
		rowText = colors;
	} else {
		rowText = NULL;
	}
	plotDendroAndColors(dendro, colors, groupLabels=groupLabels, rowText=rowText, main=main, dendroLabels=F, hang=hang, addGuide=addGuide, guideHang=guideHang);
}

# ########################################################################

exn.computeGS = function(traits, DATA, verbose=T, ...) {
	nrow = ncol(DATA);
	ncol = ncol(traits);
	rownames = names(DATA);
	colnames = colnames(traits);
	GSmat = as.data.frame(matrix(nrow=nrow, ncol=ncol, dimnames=list(rownames, paste('GS.', colnames, sep=''))));
	Pmat = as.data.frame(matrix(nrow=nrow, ncol=ncol, dimnames=list(rownames, paste('p.GS.', colnames, sep=''))));
	Qmat = as.data.frame(matrix(nrow=nrow, ncol=ncol, dimnames=list(rownames, paste('q.GS.', colnames, sep=''))));
	for (t in 1:ncol(traits)) {
		if (verbose) {cat('working on "', names(traits)[t], '"\n', sep='')};
		GSmat[, t] = as.data.frame(cor(DATA, traits[, t], use='p'));
		if (any(is.nan(GSmat[, t]))) {
			warning('Some p-vals are NaN, setting them equal to 1');
			GSmat[is.nan(GSmat[, t]), t] = 0.000001;
		}
		Pmat[, t] = as.data.frame(corPvalueFisher(as.matrix(GSmat[, t]), nrow(DATA)));
		if (any(is.nan(Pmat[, t]))) {
			warning('Some p-vals are NaN, setting them equal to 1');
			Pmat[is.nan(Pmat[, t]), t] = 1;
		}
		Qmat[, t] = as.data.frame(qvalue(Pmat[, t], ...)$qvalues);
		if (verbose) {
			cat('# of pvals < 0.05 ...', sum(Pmat[, t]<0.05), '\n');
			cat('# of qvals < 0.05 ...', sum(Qmat[, t]<0.05), '\n');
		}
	}
	tmp = as.data.frame(matrix(nrow=nrow, ncol=ncol*3), row.names=rownames);
	for (col in 1:ncol) {
		tmpcol = (col*3)-2;
		tmp[, tmpcol:(tmpcol+2)] = as.data.frame(cbind(GSmat[, col], Pmat[, col], Qmat[, col]));
		names(tmp)[tmpcol:(tmpcol+2)] = c(names(GSmat)[col], names(Pmat)[col], names(Qmat)[col]);
	}
	return(tmp);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.fundamentalModuleConcepts = function(colors, DATA, type='signed', power=18, verbose=T) {
	modNames = names(table(colors));
	modFNC = list();
	if (verbose) {cat('Working on...\n')};
	for (mod in 1:length(modNames)) {
		if (verbose) {cat(modNames[mod], ' ', sep='')};
		modFNC[[mod]] = fundamentalNetworkConcepts(adjacency(DATA[colors==modNames[mod]], type=type, power=power));
		names(modFNC)[mod] = modNames[mod];
	}
	return(modFNC);
}

# ########################################################################

exn.plotModuleHeatmap = function(module, colors, DATA, power=18, type='signed', modGeneInfoList=NULL, num=NULL, oma=c(5,0,0,5)) {
	par(oma=oma);
	if (!is.null(modGeneInfoList)) {
		names = names(DATA)[colors==module];
		if (any(names(modGeneInfoList)!=c('genes', 'geneAnnos', 'geneSyms', 'geneIDs'))) {
			labCol = names;
			labRow = names;
			warning('Cannot map gene symbols to ids since modGeneInfoList does not match format of output from exn.buildModuleGeneInfo()');
		}
		ind = match(module, names(modGeneInfoList$genes));
		if (any(names!=modGeneInfoList$genes[[ind]])) {
			labCol = names;
			labRow = names;
			warning('Cannot map gene symbols to ids since order does not match');
		} else {
			names = paste(modGeneInfoList$genes[[ind]], ' (', modGeneInfoList$geneSyms[[ind]], ')', sep='');
			labCol = names;
			labRow = names;
		}
	}
	adj = adjacency(as.matrix(DATA[, colors==module]), power=power, type=type);
	heatmap(adj, symm=T, labCol='', labRow='');
}

# ########################################################################

exn.buildNetworkInfoTable = function(modGeneInfoList, kMEtable, GStable, write=F, filename=NULL) {
	if (any(names(modGeneInfoList)!=c('genes', 'geneAnnos', 'geneSyms', 'geneIDs'))) {
		stop('modGeneInfoList must match format of output from exn.buildModuleGeneInfo()');
	}
	if (any(rownames(kMEtable)!=rownames(GStable))) {
		stop('kME and GS rownames do not match');
	}
	if (nrow(kMEtable)!=sum(unlist(lapply(modGeneInfoList$genes, length)))) {
		stop('number of transcripts in kMEtable does not match number in modGeneInfoList');
	}
	modMat = cbind(rep(names(modGeneInfoList$genes)[1], length(modGeneInfoList$genes[[1]])), modGeneInfoList$genes[[1]], modGeneInfoList$geneAnnos[[1]], modGeneInfoList$geneIDs[[1]]);
	colnames(modMat) = gsub('1', 'm', names(modMat));
	for (m in 2:length(modGeneInfoList$genes)) {
		modMat = rbind(modMat, cbind(rep(names(modGeneInfoList$genes)[m], length(modGeneInfoList$genes[[m]])), modGeneInfoList$genes[[m]], modGeneInfoList$geneAnnos[[m]], modGeneInfoList$geneIDs[[m]]));
	}
	names(modMat)[c(1, 2, 8)] = c('module', 'transcriptID', 'EntrezID');
	modMat = as.data.frame(cbind(modMat[match(rownames(kMEtable), modMat$transcriptID), ], GStable, kMEtable));
	rownames(modMat) = modMat$transcriptID;
	modMat = modMat[, -2];
	# if (write) {
		# write.table(modMat, file=filename, quote=F, row.names=T, col.names=T);
	# }
	return(modMat);
}

# ########################################################################
# add column to tmp GS matrix in order to plot kME against gene metric that only exists for a subset of genes in the network 
# metricVec should be vector of metric with gene names 
exn.makeTmpGS = function(metricVec, metricName, GStable) {
	metricVec = metricVec[match(rownames(GStable), names(metricVec)), ];
	tmpGS = cbind(GStable, metricVec);
	names(tmpGS)[ncol(tmpGS)]=paste('GS.', metricName, sep='');
	return(tmpGS);
}


# ########################################################################

exn.plotModuleGSkME = function(module, colors, trait, GStable, kMEtable, cex=2, main='', xlab=NULL, ylab=NULL, horiz=F, ...) {
	if (any(rownames(GStable)!=rownames(kMEtable))) {
		stop('rownames in GS and kME tables do not match');
	}
	modGenes = colors==module;
	GScol = match(paste('GS.', trait, sep=''), names(GStable)); 
	modGS = GStable[modGenes, GScol];
	ymax = max(abs(modGS),na.rm=T);
	kMEcol = match(paste('kME', module, sep=''), names(kMEtable)); 
	modkME = kMEtable[modGenes, kMEcol];
	if (is.null(ylab)) {
		ylab = names(GStable)[GScol];
	}
	if (is.null(xlab)) {
		xlab = names(kMEtable)[kMEcol];
	}
	verboseScatterplot(modkME, modGS, abline=T, abline.color='red', abline.lty='dashed', col=module, pch=20, frame.plot=F, cex=cex, ylim=c(-ymax-.1, ymax+.1), xlim=c(.5, 1), xlab=xlab, ylab=ylab, main=main, ...);
	if (horiz) {
		abline(h=0, col='darkgrey');
	}
}

# ########################################################################

exn.plotModuleAllGSkME = function(module, colors, GStable, kMEtable, mfrow=c(1,ncol(GStable)/3), cex=2, main='', xlab=NULL, ylab=NULL, horiz=F, ...) {
	if (any(rownames(GStable)!=rownames(kMEtable))) {
		stop('rownames in GS and kME tables do not match');
	}
	modGenes = colors==module;
	kMEcol = match(paste('kME', module, sep=''), names(kMEtable));
	par(mfrow=mfrow);
	for (GScol in seq(1,length(names(GStable)),3)) {
		modGS = GStable[modGenes, GScol];
		#ymax = max(abs(modGS));
		modkME = kMEtable[modGenes, kMEcol];
			ylab = names(GStable)[GScol];
			xlab = names(kMEtable)[kMEcol];
			verboseScatterplot(modkME, modGS, abline=T, abline.color='red', abline.lty='dashed', col=module, pch=20, frame.plot=F, cex=cex, ylim=c(-1, 1), xlim=c(.5, 1), xlab=xlab, ylab=ylab, main=main, ...);
			if (horiz) {
				abline(h=0, col='darkgrey');
			}
	}
}


# ########################################################################

exn.plotAllModsGSkME = function(colors, trait, GStable, kMEtable, mfrow, cex=2, order=F, xlab=NULL, ylab=NULL, subset=NULL, horiz=F, returnCors=T, plot=T, ...) {
	modNames = names(table(colors));
	if (is.numeric(subset)) {
		modNames = modNames[subset];print(modNames)
	}
	if (order | returnCors) {
		cors = c();
		for (m in 1:length(modNames)) {
			modGenes = colors==modNames[m];
			modCol = match(modNames[m], gsub('kME', '', names(kMEtable)));
			GScol = match(trait, gsub('GS.', '', names(GStable)));
			cors = c(cors, cor(kMEtable[modGenes, modCol], GStable[modGenes, GScol]));
		}
		names(cors) = modNames;
	}
	if (order) {
		modNames = modNames[order(abs(cors), decreasing=T)];
		cors = cors[order(abs(cors), decreasing=T)];
	}
	if (plot) {
		par(mfrow=mfrow, oma=c(0, 1, 2, 0));
		for (m in 1:length(modNames)) {
			module = modNames[m];
			exn.plotModuleGSkME(module, colors, trait, GStable, kMEtable, cex=cex, xlab=xlab, ylab=ylab, main=paste(module, '\n'), horiz=horiz, ...);
		}
		title(paste('GS.', trait, ' (y-axis) vs kME (x-axis)', sep=''), outer=T);
	}
	if (returnCors) {
		return(cors)
	}
}

# ########################################################################

exn.plotModkMEAndkscore = function(kME, colors, kscores, module, thresh=0, ...) {
	modkME = exn.getModulekME(module,colors,kME);
	modkscore = kscores[names(kscores) %in% rownames(modkME)];	
	modkscore = modkscore[match(rownames(modkME),names(modkscore))];#print(head(modkscore))
	modkscore[is.na(modkscore)] = 0;
	modkscore = modkscore[modkscore>thresh];#print(head(modkscore))
	verboseScatterplot(modkME[rownames(modkME) %in% names(modkscore), 1], modkscore,frame.plot=F,abline=T,xlab=names(modkME)[1], ylab='kscores', col='black',bg=module,pch=21, ...);
	return(modkscore);
}

# ########################################################################

# getModuleUniqueTerms = function() {
	
# }

# ########################################################################

# getTermsInCommon = function() {
	
# }

# ########################################################################

# ########################################################################

exn.plotPowers = function (data, networkType='signed', blockSize=ncol(data)+1, powerVector=NULL, ...) {
	sft = pickSoftThreshold(data, networkType=networkType, blockSize=blockSize, ...);
	
	par(mfrow=c(1,3))
	set = sft$fitIndices;
	x = set$Power;
	y = set$SFT.R.sq;
	plot(x, y, ylim = c(0, 1), type = 'n', xlab = 'Power', ylab = 'Scale-free fit', main = '', frame.plot=F, ...);
	text(x, y, labels = x, ...);
	abline( h = 0.8, col = 'red', lty = 'dashed', ... );
	
	y = set$mean.k.;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Mean k', main = '', frame.plot=F, ...);
	text(x, y, labels = x, ...);
	abline( h = 50, col = 'red', lty = 'dashed', ... );
	
	y = set$slope;
	plot(x, y, type = 'n', xlab = 'Power', ylab = 'Slope', main = '', frame.plot=F, ...);
	text(x, y, labels = x, ...);
	abline( h = -2, col = 'red', lty = 'dashed', ... );
	abline( h = -1, col = 'red', lty = 'dashed', ... );
	
	return(sft);
}

exn.computeConnectivityAndPlotScaleFreeness = function (data, networkType='signed', power=18, blockSize=ncol(data)+1, ...) {
	k = softConnectivity(data, type=networkType, blockSize=blockSize, power=power, ...);
	par(mfrow=c(1,2));
	scaleFreePlot(k);
	hist(k, col='lightgrey', border='darkgrey');
	return(k)
}

########################################################################
#####
##### BiomaRt functions
#####
########################################################################
#############   collect hits from one or more Marts in same data frame
# Depends on:	biomaRt::getBM()
#
## Arguments
### attributes:	character vector of attributes to retrieve from biomart, valid attributes can be found using biomaRt::listAttributes(mart)
### filters:	character vector of filters to constrain biomart query, valid filters can be found using biomaRt::listFilters(mart)
### values:		vector holding value(s) of filter(s), if multiple filters supplied should be list of vectors
### marts:		Mart object, or list of Mart objects, created via biomaRt::useMart()
#
## Value
### hits:		data frame where each row is a hit and each column is an attribute

exn.getBiomaRtHits = function(attributes, filters, values, marts) {
	hits = matrix(ncol=length(attributes));
	colnames(hits) = attributes;
	for (f in 1:length(marts)) {
		#print(names(marts)[f])
		result = getBM(attributes=attributes, filters=filters, values=values, mart=marts[[f]]);
		if (nrow(result) > 0) {
			hits = rbind(hits, result);
		}
	}
	if (sum(is.na(hits[1, ])) == ncol(hits)) {
		hits = hits[-1, ];
	}
	return(hits);
}
########################################################################

.parseMSigDBgenes = function (filename, ...) {
	rawgenes = read.table(filename, header=T, sep='\t', ...);
	return(rawgenes[3:nrow(rawgenes), 1]);
}

.parseMSigDBgenesFromMultipleFiles = function (filenames) {
	outlist = list();
	for ( file in filenames ) {
		outlist[[file]] = .parseMSigDBgenes(filename=file);
	}
	return(outlist);
}

.analyzeMSigDBgenes = function(filename, modGenesList, modkMEList, modDAVID, thresh=.05, mfrow) {
	geneSyms = .parseMSigDBgenes(filename);
	geneEnsembl = .convertIDsWithBiomaRt(fromID='hgnc_symbol', toID='ensembl_gene_id', ids=geneSyms);
	geneEnsembl = geneEnsembl[grep('^ENS',geneEnsembl$ensembl_gene_id), ];
	genesInNet = geneEnsembl$ensembl_gene_id[geneEnsembl$ensembl_gene_id %in% unlist(modGenesList)];
	genekMEs = .getkMERankAndQuantileForGenes(genesInNet, modkMEList);
	modEnrichments = checkGeneListEnrichmentList(geneEnsembl$ensembl_gene_id, modGenesList, unlist(modGenesList));
	mods = names(table(genekMEs$module));
	termEnrichments = .checkTermGeneEnrichmentAllModules(modDAVID=modDAVID, genes=genesInNet, modGenes=modGenesList, thresh=thresh);
	terms = exn.getAllTermsWithGenesAcrossModules(genesInNet, modDAVID);
	genesDAVID = exn.getDAVIDchartForGeneList(genes=genesInNet, background=unlist(modGenesList));
	#.verboseBoxplotComparekMEOfGeneSetsWithinAllModules(modkMEList, genesInNet, mfrow)
	#verboseBoxplot(as.numeric(genekMEs), genekMEs$module)
	return(list(genes=geneEnsembl, kME=genekMEs, modEnrichments=modEnrichments, termEnrichments=termEnrichments, terms=terms, genesDAVID=genesDAVID));
}


.analyzeGenes = function(genes, modGenesList, modkMEList, modDAVID, thresh=.05, mfrow) {
	genesInNet = genes[genes %in% unlist(modGenesList)];
	genekMEs = .getkMERankAndQuantileForGenes(genesInNet, modkMEList);
	modEnrichments = checkGeneListEnrichmentList(genes, modGenesList, unlist(modGenesList));
	mods = names(table(genekMEs$module));
	termEnrichments = .checkTermGeneEnrichmentAllModules(modDAVID=modDAVID, genes=genesInNet, modGenes=modGenesList, thresh=thresh);
	terms = exn.getAllTermsWithGenesAcrossModules(genesInNet, modDAVID);
	genesDAVID = exn.getDAVIDchartForGeneList(genes=genesInNet, background=unlist(modGenesList));
	#.verboseBoxplotComparekMEOfGeneSetsWithinAllModules(modkMEList, genesInNet, mfrow)
	#verboseBoxplot(as.numeric(genekMEs), genekMEs$module)
	return(list(kME=genekMEs, modEnrichments=modEnrichments, termEnrichments=termEnrichments, terms=terms, genesDAVID=genesDAVID));
}


#boxplot(as.numeric(kME) ~ with(xx$kME, reorder(module, -as.numeric(kME), median)), data=xx$kME,cex.axis=.5,col=as.character(levels(with(xx$kME, reorder(module, -as.numeric(kME), median)))))
########################################################################

# hsa = useMart('ensembl','hsapiens_gene_ensembl');
# library(GO.db);
# GOterms = as.list(GOTERM);

# getModGOfromEntrez = function(entrezIDsList) {
	# GOlist = list();
	# for (m in 1:length(entrezIDsList)) {
		# print(names(entrezIDsList)[m])
		# GOlist[[m]] = getBM(attributes=c('go_id'), filters='entrezgene', values=entrezIDsList[[m]], hsa);
		# names(GOlist)[m] = names(entrezIDsList)[m];
	# }
	
	# return(GOlist);
# }


.convertIDsWithBiomaRt = function (biomart='ensembl', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', ids) {
	mart = useMart(biomart=biomart, dataset=dataset);
	return(getBM(attributes=c(fromID, toID), filters=fromID, values=ids, mart=mart));
}

.convertIDsWithBiomaRtList = function (biomart='ensembl', dataset='hsapiens_gene_ensembl', fromID='ensembl_gene_id', toID='entrezgene', modGenes) {
	IDs = list();
	for (mod in 1:length(modGenes)) {
		cat(mod, ':', names(modGenes)[mod], ' ', sep='');
		IDs[[mod]] = .convertIDsWithBiomaRt(biomart=biomart, dataset=dataset, fromID=fromID, toID=toID, ids=modGenes[[mod]]);
		names(IDs)[mod] = names(modGenes)[mod];
	}
	cat('\n');
	return(IDs)
}

########################################################################  
############# 
# Depends on:		biomaRt::getBM()
#
## Arguments
### entrezIDsList:
### GOterms:
### IDandTermsOnly:	boolean indicating 
#
## Value
### GOlist:

exn.getModGOfromEntrez = function(entrezIDsList, GOterms, IDandTermsOnly=T) {
	GOlist = list();
	for (m in 1:length(entrezIDsList)) {
		print(names(entrezIDsList)[m])
		tmp = getBM(attributes=c('go_id'), filters='entrezgene', values=entrezIDsList[[m]], hsa);
		tmp = GOterms[names(GOterms) %in% tmp];
		if (IDandTermsOnly) {
			GOlist[[m]] = unlist(lapply(tmp, Term));
		} else {
			GOlist[[m]] = tmp;
		}
		names(GOlist)[m] = names(entrezIDsList)[m];
	}
	return(GOlist);
}

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# exn.getModGOfromEntrezByGene = function(entrezIDsVec, GOterms, attributes=c('go_id','name_1006','namespace_1003')) {
	# modGObyGene = list();
	# for (g in 1:length(entrezIDsVec)) {
		# modGObyGene[[g]] = getBM(attributes=attributes, filters='entrezgene', values=entrezIDsVec[g], hsa);
		# names(modGObyGene)[g] = names(entrezIDsVec)[g];
	# }
	# return(modGObyGene);
# }

########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

# entrezIDsVec should be vector of human Entrez IDs, each named with gene symbol
exn.getGOfromEntrezByGene = function(entrezIDsVec, GOterms, attributes=c('go_id','name_1006','namespace_1003'), write=F, verbose=T) {
	allGenesGO = list();
	for (g in 1:length(entrezIDsVec)) {
		if (verbose) {
			print(names(entrezIDsVec)[g]);
		}
		allGenesGO[[g]] = getBM(attributes=attributes, filters='entrezgene', values=entrezIDsVec[g], hsa);
		names(allGenesGO)[g] = names(entrezIDsVec)[g];
		if (write) {
			filename = paste('entrez', entrezIDsVec[g], '_', names(entrezIDsVec)[g], '_GOterms.txt', sep='');
			write.table(allGenesGO[[g]], file=filename, quote=F, row.names=F, sep='\t');
		}
	}
	return(allGenesGO);
}

########################################################################
# depends on WGCNA::verboseBoxplot
exn.plotMEsExprWithTraits = function (MEs, factors, mfrow, horiz=NULL, ...) {
	par(mfrow=mfrow);
	for (i in 1:ncol(MEs)) {
		verboseBoxplot(MEs[,i], as.factor(factors), xlab='', ylab='', col=gsub('ME','',names(MEs)[i]),main=names(MEs)[i], ...);
		if (is.numeric(horiz)) {
			abline(h=horiz, col='red', lty='dashed');
		}
	}
}

exn.getNetworkBasics = function (net, data, suffix='.1') {
	assign(paste('dendro', suffix, sep=''), net$dendrograms, envir=.GlobalEnv);
	assign(paste('blockGenes', suffix, sep=''), net$blockGenes, envir=.GlobalEnv);
	assign(paste('colors', suffix, sep=''), net$colors, envir=.GlobalEnv);
	assign(paste('MEs', suffix, sep=''), net$MEs, envir=.GlobalEnv);
	assign(paste('kME', suffix, sep=''), exn.computekME(data, net$MEs)$all, envir=.GlobalEnv);
	assign(paste('modGenes', suffix, sep=''), exn.getModuleGenes(data, net$colors), envir=.GlobalEnv);
	assign(paste('modkMEs', suffix, sep=''), 
		   .getModGenesRankedBykME(names(table(net$colors)),
		   						   net$colors,
		   						   exn.computekME(data, net$MEs)$all
		   						   ), 
		   	envir=.GlobalEnv);
}


exn.computeModulePreservation = function (data1, data2, colors1, colors2, labels=c('1', '2'), networkType='signed', ranksOnly=T, nPermutations=100, verbose=3) {
	multiExpr = list();
	multiExpr[[1]] = list(data=data1);
	multiExpr[[2]] = list(data=data2);
	names(multiExpr) = labels;
	colorList=list();
	colorList[[1]] = colors1;
	colorList[[2]] = colors2;
	names(colorList) = labels;
	
	if (ranksOnly) {
		mp = modulePreservation(multiExpr, colorList, networkType=networkType, nPermutations=0, maxModuleSize=max(table(colors1)), verbose=verbose, savePermutedStatistics=F);
		#tmp = cbind(get(intra,pos=-1), get(inter));
		tmp = cbind(mp[[1]][[1]][[3]][[2]], mp[[1]][[1]][[4]][[2]])
		tmp = tmp[!(rownames(tmp)=='gold'), ];
		tmp2 = apply(-tmp, 2, rank);
		# con - cor.kIM, cor.kME, cor.cor
		medRankCon = apply( tmp2[ , c(8, 9, 11) ], 1, median);
		# density - propVarExplained, meanSignAwareKME, meanSignAwareCorDat, meanAdj
		medRankDen = apply(tmp2[ , c(1, 2, 4, 5) ], 1, median);
		medRankPres = (medRankCon + medRankDen) / 2;
		medRanks = cbind(medRankPres, medRankCon, medRankDen);
		medRanks = medRanks[order(medRanks[,1]), ];
		return(list(medRanks=medRanks, mp=mp));
	} else {
		mp = modulePreservation(multiExpr, colorList, networkType=networkType, maxModuleSize=max(table(colors1)), nPermutations=nPermutations, verbose=verbose, savePermutedStatistics=F);
		return(mp);
	}
}


exn.getDAVIDchartForGeneList = function (genes, background, idType='ENSEMBL_GENE_ID', login='ahilliar@stanford.edu') {
	david = DAVIDWebService$new(email=login);
	addList(david, inputIds=background, idType=idType, listName='BG', listType='Background');
	addList(david, inputIds=genes, idType=idType, listName='genes', listType='Gene');
	setCurrentBackgroundPosition(david, which(getBackgroundListNames(david) == 'BG'));
	setCurrentGeneListPosition(david, 1);
	chart = getFunctionalAnnotationChart(david);
	checkFactor = which(lapply(chart, class)=='factor');
	if (length(checkFactor)>0) {
		chart[, checkFactor] = as.character(chart[, checkFactor]);
	}
	return(chart);	
}

# right now uploads background everytime, need to fix
exn.getDAVIDchartsForModuleSubsetsInOtherNetwork = function (moduleRef, modulesTest, modGenesRef, modGenesTest, idType='ENSEMBL_GENE_ID', login='ahilliar@stanford.edu') {
	if (any(sort(unlist(modGenesRef)) != sort(unlist(modGenesTest)))) {stop()}
	out = list();
	refgenes = modGenesRef[[match(moduleRef, names(modGenesRef))]];
	bg = unlist(modGenesRef);
	for (m in 1:length(modulesTest)) {
		cat(modulesTest[m], '\n');
		testgenes = modGenesTest[[match(modulesTest[m], names(modGenesTest))]];
		out[[m]] = exn.getDAVIDchartForGeneList(genes=refgenes[refgenes %in% testgenes], background=bg, idType=idType, login=login);
		names(out)[m] = modulesTest[m];
	}
	return(out);
}

# depends on circlePlot.R
exn.circlePlotGenes = function (genes, data, colors, type='signed', power=18, returnAdj=F, linecol.base=.9, linecol.gamma=1, cex.labels=c(.5,1.2), min.cex.points=1, max.cex.points=2, orderBy=NULL, ...) {
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
			   min.cex.points=min.cex.points, max.cex.points=max.cex.points, ...
			   );
	if (returnAdj) {
		return(adj);
	}
}

exn.circlePlotCompareGenesAcrossNetworks =  function (genes, data1, data2, colors1, colors2, type='signed', power=18, linecol.base=.9, linecol.gamma=1, cex.labels=c(.5,1.2), min.cex.points=1, max.cex.points=2, return=T) {
	gdata1 = data1[, match(genes, names(data1))];
	gdata2 = data2[, match(genes, names(data2))];
	adj1 = adjacency(gdata1, type=type, power=power);
	adj2 = adjacency(gdata2, type=type, power=power);
	colors1 = colors1[match(genes, names(data1))];
	colors2 = colors2[match(genes, names(data2))];
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
	if (return) {
		return(list(adj1=adj1, adj2=adj2, cols1=colors1, cols2=colors2));
	}
}


.buildModuleSubsetsFromOverlap = function (overlapPvals, modGenesList1, modGenesList2, thresh=NULL) {
	if (is.null(thresh)) {
		thresh = .05 / ( nrow(overlapPvals)*ncol(overlapPvals) );
	}
	modSubsets = list();
	for (row in 1:nrow(overlapPvals)) {
		mod = rownames(overlapPvals)[row];
		omods = colnames(overlapPvals)[overlapPvals[row, ] < thresh];
		if (length(omods)==0) {
			modSubsets[[row]] = list(NULL);
			names(modSubsets)[row] = mod;
			next;
		}
		modgenes = modGenesList1[match(mod, names(modGenesList1))][[1]];
		modlist = list();
		for (m in 1:length(omods)) {
			omodgenes = modGenesList2[match(omods[m], names(modGenesList2))][[1]];
			modlist[[m]] = intersect(modgenes, omodgenes);
			names(modlist)[m] = omods[m];
		}
		modSubsets[[row]] = modlist;
		names(modSubsets)[row] = mod;
	}
	return(modSubsets);
}

.getNumsForModuleOverlapSubsets = function (modSubsets, nums) {
	modNums = modSubsets;
	for (m in 1:length(modSubsets)) {
		this = modSubsets[[m]]
		if (is.null(this)) {next}
		for (subm in 1:length(this)) {
			thesegenes = this[[subm]];
			tmp = nums[match(thesegenes, names(nums))];
			names(tmp) = thesegenes;
			modNums[[m]][[subm]] = tmp;
		}
	}
	return(modNums);
}

.getOriginalkMEForModuleOverlapSubsets = function (modSubsets, modkMEsList) {
	subkMEs = modSubsets;
	for (m in 1:length(modSubsets)) {
		this = modSubsets[[m]];
		mod = names(subkMEs)[m];	
		modkME = modkMEsList[[match(mod, names(modkMEsList))]];
		for (sm in 1:length(this)) {
			if (length(this)==1 & is.null(this[[sm]])) {next}
			smgenes = this[[sm]];
			subkMEs[[m]][[sm]] = modkME[match(smgenes, rownames(modkME)), 1];
			names(subkMEs[[m]][[sm]]) = rownames(modkME)[match(smgenes, rownames(modkME))];
		}
	}
	return(subkMEs);
}

.verboseBoxplotCompareModuleSubsets = function (modSubsetNums, focalMods, ...) {
	fac = sapply(strsplit(names(unlist(modSubsetNums)), '.', fixed=T), function(f) f[1]);
	fac[!(fac %in% focalMods)] = '_rest_';
	tmp = table(fac);
	for (i in 1:length(tmp)) {
		fac[fac==names(tmp)[i]] = paste(fac[fac==names(tmp)[i]], ' (', tmp[i], ')', sep='');
	}
	verboseBoxplot(unlist(modSubsetNums), as.factor(fac), 
				   notch=F, xlab='', 
				   border='darkgrey', ...);
}

.verboseBoxplotComparekMEOfGeneSetsWithinModule = function (modkME, genes, ...) {
	#listname = deparse(substitute(genes));
	#names = c(paste('not ', listname, sep=''), listname);
	verboseBoxplot(modkME[,1], rownames(modkME) %in% genes, frame.plot=F, ylab='kME', xlab='', ...)
}

.verboseBoxplotComparekMEOfGeneSetsWithinAllModules = function (modkMEs, genes, mfrow, ...) {
	par(mfrow=mfrow);
	for (m in 1:length(modkMEs)) {
		if (sum(genes %in% rownames(modkMEs[[m]]))==0) {next}
		.verboseBoxplotComparekMEOfGeneSetsWithinModule(modkMEs[[m]], genes, col=names(modkMEs)[m], main=names(modkMEs)[m], ...);
	}
}

# num should be named numeric vector
.verboseBoxplotAcrossModules = function (modGenes, num, factorgenes, mfrow, names=NULL, ...) {
	par(mfrow=mfrow)
	if (is.null(names)) {
		names = c(paste('not ', deparse(substitute(factorgenes)), sep=''), deparse(substitute(factorgenes)));
	}
	for (m in 1:length(modGenes)) {
		mgenes = modGenes[[m]];
		mnums = num[match(mgenes, names(num))];
		mfact = mgenes %in% factorgenes;
		if (sum(mfact) %in% c(0, 1, (length(mfact)-1))) {next}
		verboseBoxplot(mnums, mfact, col=names(modGenes)[m], frame.plot=F, ylab=deparse(substitute(num)), xlab='', main=names(modGenes)[m], names=names, ...)
	}
}

.verboseScatterplotkMEAgainstNumAllMods = function (modkMEs, modNum, mfrow, ...) {
	if (any(names(modNum) != names(modkMEs))) {stop()}
	par(mfrow=mfrow);
	for(m in 1:length(modkMEs)) {
		verboseScatterplot(modkMEs[[m]][,1], modNum[[m]][match(rownames(modkMEs[[m]]), names(modNum[[m]]))],
						   abline=T, frame.plot=F, col=names(modkMEs)[m], pch=19, xlab=deparse(substitute(modkMEs)), ylab=deparse(substitute(modNum)),...)
	}
}

.plotComparekMEAndExpressionAcrossNetworks = function (homekMEsVec1, homekMEsVec2, data1, data2, colors1, subsetGenes=NULL, subsetModules=NULL, ...) {
	
	plotkME1 = homekMEsVec1;
	plotkME2 = homekMEsVec2;
	plotDATA1 = data1;
	plotDATA2 = data2;
	plotcolors = colors1;
	if (!is.null(subsetModules)) {
		plotkME1 = plotkME1[!(colors1 %in% subsetModules)];
		plotkME2 = plotkME2[!(colors1 %in% subsetModules)];
		plotDATA1 = plotDATA1[, !(colors1 %in% subsetModules)];
		plotDATA2 = plotDATA2[, !(colors1 %in% subsetModules)];
		plotcolors = colors1[!(colors1 %in% subsetModules)];
	}
	
	par(mfrow=c(1,2));
	verboseScatterplot(plotkME1, plotkME2, abline=T, pch=19, col=plotcolors, cex=1.1, frame.plot=F, ...); 
	verboseScatterplot(apply(plotDATA1, 2, mean), apply(plotDATA2, 2, mean), abline=T, pch=19, col=plotcolors, cex=1.1, frame.plot=F, ...);
}


.computekIN = function (module, data, colors, type='signed', power=18, blockSize=(ncol(data)+1), normalize=T, ...) {
	genes = colors==module;
	names = names(data)[genes];
	kIN = softConnectivity(data[, genes], type=type, power=power, blockSize=blockSize, ...);
	names(kIN) = names;
	if (normalize) {
		kIN = kIN/max(kIN);
	}
	return(kIN);
}

.computekINAllModules = function(data, colors, type='signed', power=18, blockSize=(ncol(data)+1), normalize=T, returnVec=F, ...) {
	kINlist = list();
	mods = names(table(colors));
	for (m in 1:length(mods)) {
		kINlist[[m]] = .computekIN(module=mods[m], data=data, colors=colors, type=type, power=power, blockSize=blockSize, normalize=normalize, ...);
		names(kINlist)[m] = mods[m];
	}
	if (returnVec) {
		tmp = unlist(kINlist);
		tmpNames = unlist(strsplit(names(tmp),'\\.'));
		names(tmp) = tmpNames[seq(2, length(tmpNames), 2)];
		return(tmp[match(names(data), names(tmp))]);
	} else {
		return(kINlist);
	}
}

.computekINdiff = function (dataRef, data2, colorsRef, colors2=NULL, type='signed', power=18, blockSize=(ncol(dataRef)+1), ...) {
	if (any(names(data2) != names(dataRef))) {
		stop('datasets must have same genes in same order');
	}
	cat('computing... \n', sep='');
	cat('  kIN in dataRef... ', sep='');
	kINref = .computekINAllModules(data=dataRef, colors=colorsRef, type=type, power=power, blockSize=blockSize, normalize=T, returnVec=T, ...);
	cat('kIN in data2 using colorsRef... ', sep='');
	kINref.pres = .computekINAllModules(data=data2, colors=colorsRef, type=type, power=power, blockSize=blockSize, normalize=T, returnVec=T, ...);
	kIN.pres = kINref - kINref.pres;
	if (!is.null(colors2)) {
		cat('kIN in data2 using colors2... ', sep='');
		kIN2 = .computekINAllModules(data=data2, colors=colors2, type=type, power=power, blockSize=blockSize, normalize=T, returnVec=T, ...);
		kIN.diff = kINref - kIN2;
		return(as.data.frame(cbind(kIN.ref=kINref, kIN.2.ref=kINref.pres, kIN.2=kIN2, kIN.pres=kIN.pres, kIN.diff=kIN.diff, moduleRef=colorsRef, module2=colors2)));
	} else {
		return(as.data.frame(cbind(kIN.ref=kINref, kIN.2.ref=kINref.pres, kIN.2=kIN2, kIN.pres=kIN.pres, moduleRef=colorsRef)));
	}
}

.getkMERankAndQuantileForGenes = function (genes, modkMEList) {
	mat = as.data.frame(matrix(nrow=length(genes), ncol=4));
	dimnames(mat) = list(genes, c('module','kME','rank','quantile'));
	for (g in 1:length(genes)) {
		gene = genes[g];
		tmp = lapply(modkMEList, function(f) if(any(grepl(gene, rownames(f)))){f});
		tmp = tmp[!sapply(tmp, is.null)][[1]];
		module = gsub('kME', '', names(tmp)[1]);
		row = which(rownames(tmp)==gene);
		kME = tmp[row, 1];
		qqtile = row / nrow(tmp);
		mat[g, ] = c(module, signif(kME,3), row, signif(qqtile,3));
	}
	return(mat);
}
#####################

.geneListIntoModules = function (genes, modGenesList) {
	x = list();
	for (mod in names(modGenesList)) {
		x[[mod]] = genes[genes %in% modGenesList[[mod]]];
	}
	return(x);
}

.numVecIntoModules = function (numVec, modGenesList) {
	x = list();
	for (mod in names(modGenesList)) {
		x[[mod]] = numVec[names(numVec) %in% modGenesList[[mod]]];
	}
	return(x);
}

###############

.verboseScatterplotAllColumnPairs = function (mat, mfrow=c(3,ncol(mat)), col=NULL, ...) {
	if (is.null(col)) {
		col = rownames(mat);
	} else {
		col = col;
	}
	par(mfrow=mfrow);
	for (i in 1:ncol(mat)) {
		if (i==ncol(mat)) {break};
		for (j in (i+1):ncol(mat)) {
			verboseScatterplot(mat[,i], mat[,j], 
							   xlab=colnames(mat)[i], ylab=colnames(mat)[j],
							   frame.plot=F, col=col,
							   abline=T, ...
							   );
		}
	}
}

.verboseScatterplotVecAgainstColumns = function (mat, vec, ...) {
	if (any(names(vec) != rownames(mat))) {
		stop('names(vec) must be identical to rownames(mat)');
	}
	par(mfrow=c(2,ceiling(ncol(mat)/2)));
	for (i in 1:ncol(mat)) {
		verboseScatterplot(vec, mat[,i], 
						   xlab=deparse(substitute(vec)), ylab=colnames(mat)[i],
						   frame.plot=F, abline=T, ...
						   );
	}
}
# num vec must have gene names and must be in same order as colors
.verboseScatterplotNumsAcrossModules = function (num1, num2, colors, mfrow, returnCors=T, onlyLabel1=T, ...) {
	modules = names(table(colors));
	cors = c();
	#cors = vector(mode='numeric', length=length(modules));
	#names(cors) = modules;
	par(mfrow=mfrow);
	for (m in modules) {
		mgenes = colors==m;
		n1 = num1[mgenes];
		n2 = num2[mgenes];
		if (onlyLabel1 && m==modules[1]) {
			xlab=deparse(substitute(num1));
			ylab=deparse(substitute(num2))
		} else {
			xlab='';
			ylab='';
		}
		verboseScatterplot(n1, n2, xlab=xlab, ylab=ylab, main=m, 
						   frame.plot=F, abline=T, abline.col='darkgrey',
						   col=m, pch=19, ...
						   );
		cors[m] = cor(n1, n2);
	}
	if (returnCors) {
		return(cors);
	}
}

.histColumns = function (mat, mfrow=c(2,ceiling(ncol(mat)/2)), ...) {
	par(mfrow=mfrow);
	for (i in 1:ncol(mat)) {
		hist(mat[,i], main=colnames(mat)[i], col='grey', border='darkgrey', xlab='', ...);
		abline(v=median(mat[,i]), col='red');
	}
}

.addSectorsBasedOnColumns = function (mat, columns=c(4,5)) {
	mat = as.data.frame(cbind(mat, sector=rep(NA,nrow(mat))));
	for (row in 1:nrow(mat)) {
		c1 = mat[row, columns[1]] > 0;
		c2 = mat[row, columns[2]] > 0;
		if (!c1) {
			if (c2) {
				mat$sector[row] = 'A';
			} else {
				mat$sector[row] = 'C';
			}
		} else {
			if (c2) {
				mat$sector[row] = 'B';
			} else {
				mat$sector[row] = 'D';
			}
		}
	}
	return(mat);
}

.verboseBoxplotColumns = function (mat, factor, mfrow=NULL, hline=NULL, ...) {
	if (is.null(mfrow)) {
		mfrow = c(2, ceiling(ncol(mat)/2));
	}
	par(mfrow=mfrow);
	for (i in 1:ncol(mat)) {
		verboseBoxplot(mat[,i], as.factor(factor), ylab=colnames(mat)[i], xlab='', ...);
		if (is.numeric(hline)) {
			abline(h=hline);
		}
	}
}

# assumes GS is formatted like output from exn.computeGS
# i.e. in groups of 3 columns for each trait, first cor, then pval, then qval
.GSthresh = function (GS, corCols=grep('^GS',names(GS)), filterCols=NULL, thresh=.05, signed=T) {
	if (is.null(filterCols)) {
		filterCols = corCols + 2;
	} else if (!(is.vector(filterCols, mode='numeric'))) {
		stop('filterCols must be NULL or numeric vector');
	}
	th = GS[, filterCols];
	sigs = which(th <= thresh, arr.ind=T);
	nsigs = which(th > thresh, arr.ind=T);
	th[sigs] = as.integer(1);
	th[nsigs] = as.integer(0);
	if (signed) {
		if (nrow(GS) > 1000) {
			cat('Thresholding gene...');
		}     
		for (row in 1:nrow(th)) {
			if (row %% 1000 == 0) {
				cat(' ', row, sep='');
			}
			th[row, GS[row, corCols] < 0 & th[row, ] == 1] = as.integer(-1);
		}
		cat('\n');
	}
	return(th);
}

# assumes q-val column names in GS have format q.GS.trait
.groupGSthreshOutput = function (GSthreshOutput) {
	gr = c();
	th = GSthreshOutput;
	cat(' Grouping gene...');
	for (row in 1:nrow(th)) {
		if (row %% 1000 == 0) {
			cat(' ', row, sep='');
		}
		if (all(th[row, ] == 0)) {
			toPaste = '';
		} else {
			toPaste = c();
			for (i in 1:ncol(th)) {
				if (th[row, i] == 1) {
					toPaste = c(toPaste, unlist(strsplit(colnames(th)[i], '\\.'))[3]);
				}
				if (th[row, i] == -1) {
					toPaste = c(toPaste, paste('-', unlist(strsplit(colnames(th)[i], '\\.'))[3], sep=''));
				}
			}
			toPaste = paste(toPaste, collapse=':');
		}
		gr = c(gr, toPaste);
	}
	names(gr) = rownames(GSthreshOutput);
	return(gr);
}
## very inefficient now since loops twice through genes
.groupMultiGSgenesByTraits = function (GS, corCols=grep('^GS',names(GS)), filterCols=NULL, thresh=.05, signed=T) {
	th = .GSthresh(GS=GS, corCols=corCols, filterCols=filterCols, thresh=thresh, signed=signed);
	gr = .groupGSthreshOutput(th);
	names(th) = gsub('q.', '', names(th));
	return(list(gr=gr, GSth=th, numsig=apply(th, 1, function(f) sum(abs(f)))));
}

# .groupMultiGSgenesByTraits2 = function (GS, corCols=grep('^GS',names(GS)), filterCols=NULL, thresh=.05, signed=T) {
	# gr = c();
	
# }

########################################################################
## 


########################################################################

.GOFunctionAllOntologies = function (genes, refGenes, ...) {
		cat('\n==========  BP  ==========\n');
		tmpBP = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='BP', ...);print(tmpBP)
		if (!is.null(tmpBP)) {
			tmpBP = as.data.frame(cbind(ontology=rep('BP',nrow(tmpBP)), tmpBP));
		}
		cat('\n==========  CC  ==========\n');
		tmpCC = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='CC', ...);print(tmpCC)
		if (!is.null(tmpCC)) {
			tmpCC = as.data.frame(cbind(ontology=rep('CC',nrow(tmpCC)), tmpCC));
		}
		cat('\n==========  MF  ==========\n');
		tmpMF = GOFunction(interestGenes=genes, refGenes=refGenes, ontology='MF', ...);print(tmpMF)
		if (!is.null(tmpMF)) {
			tmpMF = as.data.frame(cbind(ontology=rep('MF',nrow(tmpMF)), tmpMF));
		}
		tmpAll = as.data.frame(matrix(ncol=7, dimnames=list(1, c('ontology', 'goid', 'name', 'refnum', 'interestnum', 'pvalue', 'adjustp'))))
		for (res in c('tmpBP','tmpCC','tmpMF')) {
			if (!is.null(get(res))) {
				tmpAll = as.data.frame(rbind(tmpAll, get(res)))
			}
		}
		tmpAll = tmpAll[-1, ];
		return(tmpAll[order(tmpAll$adjustp, -tmpAll$interestnum), ]);
}

.GOFunctionModules = function (modGenes, refGenes, ...) {
	out = list();
	for (m in 1:length(modGenes)) {
		cat('\n============================== ', names(modGenes)[m], '  ==============================\n', sep='');
		tmpAll = .GOFunctionAllOntologies(genes=modGenes[[m]], refGenes=refGenes, ...);
		out[[names(modGenes)[m]]] = tmpAll;
	}
	return(out);
}



#org.Hs.egGO2ALLEGSlist = as.list(org.Hs.egGO2ALLEGS);
.getGenesForGOTerms = function (terms, maps=org.Hs.egGO2ALLEGSlist, genes=NULL) {
	glist = list();
	for (t in terms) {
		g = maps[[t]];
		if (!is.null(genes)) {
			glist[[t]] = intersect(genes, g);
		} else {
			glist[[t]] = g;					### may include duplicates
		}
	}
	return(glist);
}

# x3 = .getGenesForGOTerms(terms=gofCC$turquoise$goid, genes=modGenes.hsEntrez$turquoise[,2])
# x3ENS = lapply(x3, function(f) modGenes.hsEntrez$turquoise[modGenes.hsEntrez$turquoise[,2] %in% f, 1]);

# gene ids should be Entrez
.addGenesToGOFunctionResults = function (GOFunctionAllOntologiesOutput, maps=org.Hs.egGO2ALLEGSlist, modulegenes) {
	genelist = .getGenesForGOTerms(GOFunctionAllOntologiesOutput$goid, genes=modulegenes);
	new = as.data.frame(cbind(GOFunctionAllOntologiesOutput, genes=rep(NA, nrow(GOFunctionAllOntologiesOutput))));
	for (row in 1:nrow(new)) {
		tgenes = genelist[[match(new$goid[row], names(genelist))]];
		new$genes[row] = paste(tgenes, collapse=',');
	}
	return(new);
}

#tmp = .addGenesToGOFunctionResultsList(tmp, modGenesList=lapply(modGenes.hsEntrez, function(f) f[,2]))

.addGenesToGOFunctionResultsList = function (GOFunctionModulesOutput, maps=org.Hs.egGO2ALLEGSlist, modGenesList) {
	newlist = GOFunctionModulesOutput;
	for (m in 1:length(newlist)) {
		if (nrow(newlist[[m]]) == 0) {
			cat('skipping ', names(newlist)[m], ' since no terms\n', sep='');
			next;
		}
		newlist[[m]] = .addGenesToGOFunctionResults(newlist[[m]], modulegenes=modGenesList[[m]]);
	}
	return(newlist);
}

.addNumToGOFunctionResults = function (GOFunctionAllOntologiesOutput, nums, entrezToENS=NULL, f='median', colname=NULL) {
	if (any(names(GOFunctionAllOntologiesOutput)=='genes')) {
		new = as.data.frame(cbind(GOFunctionAllOntologiesOutput, rep(NA, nrow(GOFunctionAllOntologiesOutput))));
		for (row in 1:nrow(new)) {
			genes = unlist(strsplit(new$genes[row], ','));
			if (!is.null(entrezToENS)) {
				genes = entrezToENS[match(genes, entrezToENS[,2]), 1];
			}
			new[row, ncol(new)] = .getSomeNumForGenes(genes=genes, nums=nums, f=f);
		}
	} else {
		stop('add genes to results first');
	}
	if (is.character(colname)) {
		names(new)[ncol(new)] = colname;
	} else {
		names(new)[ncol(new)] = paste(deparse(substitute(nums)), '.', f, sep='');
	}
	return(new);
}

# tmp = .addNumToGOFunctionResultsList(tmp, kscoresHuman, modGenes.hsEntrez)
# tmp = .addNumToGOFunctionResultsList(tmp, kIN.hsVec, modGenes.hsEntrez)
# tmp = .addNumToGOFunctionResultsList(tmp, kIN.diff, modGenes.hsEntrez)
# tmp = .addNumToGOFunctionResultsList(tmp, .getkMEsFromAssignedModules(modkMEs.hs, data=DATAhs), modGenes.hsEntrez, colname='kME.hs.median')
# tmp = .addNumToGOFunctionResultsList(tmp, .matColToNamedVec(GS.hs,'GS.br'), modGenes.hsEntrez, colname='GS.br.median')


.addNumToGOFunctionResultsList = function (GOFunctionModulesOutputWithGenes, nums, entrezToENS=NULL, f='median', colname=NULL) {
	newlist = GOFunctionModulesOutputWithGenes;
	if (!is.null(entrezToENS)) {
		if (!is.list(entrezToENS) | any(names(newlist) != names(entrezToENS))) {
			stop('entrezToENS should be list in same order as GOFunctionModulesOutputWithGenes');
		}
	}
	if (is.character(colname)) {
		colname = colname;
	} else {
		colname = paste(deparse(substitute(nums)), '.', f, sep='');
	}
	for (m in 1:length(newlist)) {
		if (nrow(newlist[[m]]) == 0) {
			cat('skipping ', names(newlist)[m], ' since no terms\n', sep='');
			next;
		}
		newlist[[m]] = .addNumToGOFunctionResults(GOFunctionAllOntologiesOutput=newlist[[m]], nums=nums, entrezToENS=entrezToENS[[m]], f=f, colname=colname);
	}
	return(newlist);
}



.convertmodGOFunctionResultsListToDataFrame = function (GOFunctionModulesOutput) {
	tmp = GOFunctionModulesOutput;
	df = as.data.frame(cbind(tmp[[1]], module=rep(names(tmp)[1], nrow(tmp[[1]]))));
	for (m in 2:length(tmp)) {
		df = as.data.frame(rbind(df, cbind(tmp[[m]], module=rep(names(tmp)[m], nrow(tmp[[m]])))))
	}
	return(df)
}


########################################################################
########################################################################
########################################################################
########################################################################
#############   
# Depends on:
#
## Arguments
### 
### 
#
## Value
### 

exn.exploreBlockwiseNetwork = function(net, DATA, traits, annotations, IDs, type='signed', power=18, DAVIDlogin='ahilliar@stanford.edu', idType='ENTREZ_GENE_ID', BGname='netBG', writeInfoTable=F, fileBase=NULL, verbose=T, mfrow=c(6,10), ...) {
	if (any(names(net)!=c('colors', 'unmergedColors', 'MEs', 'goodSamples', 'goodGenes', 'dendrograms', 'TOMFiles', 'blockGenes', 'blocks', 'MEsOK'))) {
		stop('net must be output from blockwiseModules() or blockwiseModulesEnriched()');
	}
	args = match.call();
	if (is.null(fileBase)) {
		fileBase = paste(gsub('()', '', args[2]), '_', paste(strsplit(as.character(Sys.time()),' ')[[1]], collapse='_'), sep='');
	}
	colors = net$colors;
	MEs = net$MEs;
	dendros = net$dendrograms;
	blocks = net$blocks;
	blockGenes = net$blockGenes;
	
	
	if (verbose) {cat('...Gathering module and transcript annotations\n  ', sep='')};
	mod = exn.buildModuleGeneInfo(DATA, colors, annotations, IDs);
	trDist = exn.transcriptsAcrossModules(names(DATA), colors);
	all = exn.getAllGeneSymbolsAndHumanEntrezIDsForNetwork(DATA, annotations, IDs, write=F);
	
	# compute kME
	if (verbose) {cat('\n...Computing kME\n', sep='')};
	kME = exn.computekME(DATA, MEs, modules=F)$all;
	
	# compute GS
	if (verbose) {cat('\n...Computing GS\n', sep='')};
	GS = exn.computeGS(traits, DATA, verbose=verbose);
	
	# compute fundamental network concepts
	if (verbose) {cat('\n...Computing fundamental network concepts for modules\n', sep='')};
	modFNC = exn.fundamentalModuleConcepts(colors, DATA, type=type, power=power, verbose=verbose);
	
	# get DAVID results
	if (verbose) {cat('\n...Sending module EntrezID lists to DAVID webserver\n', sep='')};
	david = exn.startAndUploadBgAndModListsToDAVID(login=DAVIDlogin, moduleGeneInfoList=mod, idType=idType, BG=all$geneIDs, BGname=BGname, setBG=T, verbose=verbose);
	if (verbose) {cat('  ...getting module functional annotation charts\n', sep='')};
	modDAVID0 = exn.getModChartsFromDAVID(david, setBgPosition=NULL, factorsToChars=T, verbose=verbose);
	if (verbose) {cat('\n\n...Adding network information to charts\n', sep='')};
	modDAVID = exn.addGeneInfoToDAVIDChartList(modDAVID0, IDs, mod, colors, kME, DATA, verbose=verbose);
	if (verbose) {cat('  ...finding which terms are unique or shared among modules\n', sep='')};
	termDAVID = exn.modDAVIDByTerm(modDAVID);
	
	# build network info table
	if (verbose) {cat('\n...Building network info table\n', sep='')};
	info = exn.buildNetworkInfoTable(mod, kME, GS, write=writeInfoTable, filename=fileBase);
	
	# plot dendrograms
	if (verbose) {cat('\n...Saving dendrograms as "', paste(fileBase, '_dendro-blockX.jpg', sep=''), '"\n', sep='')};
	for (b in 1:length(blockGenes)) {
		jpeg(file=paste(fileBase, '_dendro-block', b, '.jpg', sep=''), width=15, height=5, units='in', quality=100, type='quartz', res=150);
		exn.plotDendroAndColors(dendros[b], colors[blockGenes[[b]]], main=paste('block ', b, ': ', length(blockGenes[[b]]), ' transcripts', sep=''), rowText=T, block=b, blockGenes=blockGenes);
		dev.off();
	}
	
	# plot ME network
	if (verbose) {cat('...Saving ME dendrogram and heatmap as "', paste(fileBase, '_MEnetwork.jpg', sep=''), '"\n', sep='')};
	jpeg(file=paste(fileBase, '_MEnetwork.jpg', sep=''), width=10, height=15, units='in', quality=100, type='quartz', res=150);
	exn.plotEigengeneNetworks2(MEs);
	dev.off();
	
	# plot ME-trait cor
	if (verbose) {cat('...Saving ME-trait correlation heatmap as "', paste(fileBase, '_ME-trait_cors.jpg', sep=''), '"\n', sep='')};
	jpeg(file=paste(fileBase, '_ME-trait_cors.jpg', sep=''), width=5, height=15, units='in', quality=100, type='quartz', res=150);
	exn.computeAndPlotMETraitCors(traits, MEs, main='ME-trait correlations')
	dev.off()
	
	# plot GS vs kME
	if (verbose) {cat('...Saving GS-kME scatterplots as "', paste(fileBase, '_GS.XXX_vs_kME.jpg', sep=''), '"\n', sep='')};
	for (tr in seq(1, ncol(GS), 3)) {
		jpeg(file=paste(fileBase, '_', names(GS)[tr], '_vs_kME.jpg', sep=''), width=20, height=15, units='in', quality=100, type='quartz', res=150);
		exn.plotAllModsGSkME(colors, gsub('GS.', '', names(GS)[tr]), GS, kME, mfrow=mfrow, order=T, cex.main=1,cex.lab=1,cex.axis=1, horiz=T, cex=1.5);
		dev.off();
	}
	
	OUT = list(modInfo=mod, trDist=trDist, kME=kME, GS=GS, fnc=modFNC, modDAVID=modDAVID, termDAVID=termDAVID, netInfo=info);
	return(OUT);
}