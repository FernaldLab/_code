# functions for iteratively removing background genes from network using WGNCA blockwiseModules() function
# Austin Hilliard, White Lab, UCLA, September 2011
# revised in Fernald lab, Stanford, February 2013

#####
# this function, Sep. 2014
blockwiseModulesEnrichedIterate = function(DATA, networkType='signed', power=14, maxBlockSize=ncol(DATA)+1, 
										   deepSplitVec=c(2,4), 
										   mergeCutHeight=0.1, 
										   minModuleSizeVec=c(10,20,40,80,100), 
										   minKMEtoStayVec=c(0.3, 0.4),
										   minCoreKMEVec=c(0.5, 0.6),
										   densityPermTest=F, skipThresh=300, ...
										   ) {
	
	for (DS in deepSplitVec) {
		for (MM in minModuleSizeVec) {
			for (MkME in minKMEtoStayVec) {
				for (MCkME in minCoreKMEVec) {
					saveFileBase = paste(deparse(substitute(DATA)), '_',networkType,'_p',power,'_ds',DS,'_mm',MM,'_mch',mergeCutHeight,'_mKME',MkME,'_mCoreKME',MCkME,sep='');
					cat('Building network: deepSplit=', DS,
					    ', mergeCutHeight=', mergeCutHeight,
					    ', minModuleSize=', MM, 
					    ', minKMEtoStay=', MkME,
					    ', minCoreKME=', MCkME, 
					    '\n  filebase: ', saveFileBase,
					    '\n',
					    sep=''
					    );
					net = blockwiseModulesEnriched(DATA=DATA, maxBlockSize=maxBlockSize, power=power, networkType=networkType, minModuleSize=MM,
											   	   deepSplit=DS, minKMEtoStay=MkME, minCoreKME=MCkME, mergeCutHeight=mergeCutHeight,
											   	   densityPermTest=densityPermTest, skipThresh=skipThresh, verbose=0,
											  	   saveFileBase=saveFileBase, ...
											  	   );
				}
			}
		}
	}	
}



#####


blockwiseModulesEnriched = function(DATA,
								 maxBlockSize = 6000,
								 fixedBlockSize = T,
								 networkType = 'signed',
								 power = 14,
								 deepSplit = 2,
								 minModuleSize = 10,
								 smartMinModThresh = 2,
								 smartMinModMultiplier = 2, 
								 densityPermTest = T,
								 plotModDensities = T,
								 nPerm = 1000,
								 permTestPvalThresh = 0.01,
								 skipThresh = 500,
								 skipGrey = T,
								 onlyGreyThresh = 3,
								 preTOM = NULL,
								 verbose = 1,
								 saveNets = T,
								 saveDendros = T,
								 saveDATA = T,
								 savePermOut = T,
								 saveFileBase = '',
								 setRun = NULL,
								 graphicsType = 'quartz',
								 res = 150, ...) 
{
							 	
	if (ncol(DATA) < nrow(DATA)) 
	{
		DATA = as.data.frame( t(DATA) );
		cat('\n...transposing input data so genes are in columns\n');
	}
	
	bgGenes = 1;	
	onlyGrey = 0;
	if (is.numeric(setRun)) 
	{ 
		run = setRun - 1;
	} else {
		run = 0;
	}
	
	while (sum(bgGenes) > 0) 
	{
		
		numProbes = ncol(DATA);
		run = run + 1;
		cat('\nRUN ', run, '...\n', sep='');
		
		if (onlyGrey >= onlyGreyThresh) 
		{
			cat('   ....skipping perm tests in remaining runs\n');
			densityPermTest = F;
		}
		
		if (numProbes > maxBlockSize & !fixedBlockSize) 
		{
			adjustedBlockSize = (ceiling((numProbes/1000)) * 1000) / 2;
			cat('   ...maxBlockSize adjusted to ', adjustedBlockSize, '\n', sep='');
		} else {
			adjustedBlockSize = maxBlockSize;
			cat('   ...using maxBlockSize = ', adjustedBlockSize, '\n', sep='');
		}
		
		if ( (numProbes < adjustedBlockSize) & densityPermTest & is.null(preTOM) ) 
		{
			saveTOMs = T;
			saveTOMFileBase = 'tmpTOM';
		} else {
			saveTOMs = F;
			saveTOMFileBase = NULL;
		}
		
		cat('...constructing network\n');
		
		if (saveDATA) 
		{
			save(DATA, file = paste(saveFileBase, 'run', run, 'DATA.RData', sep = ''));
		}
		
		net = blockwiseModules(DATA,
							   maxBlockSize = adjustedBlockSize,
							   networkType = networkType,
							   power = power,
							   deepSplit = deepSplit,
							   minModuleSize = minModuleSize,
							   saveTOMs = saveTOMs,
							   saveTOMFileBase = saveTOMFileBase,
							   verbose = verbose, ...
							   );	
		if (saveNets) 
		{
			save(net, file = paste(saveFileBase, 'run', run, 'NET.RData', sep = ''));
		}
		
		modSizes = sort(table(net$colors), decreasing=T);
		cat('\nmodules:\n', sep='');
		print(modSizes);
		
		# if ( is.numeric(smartMinModThresh) )
		# {
			
		# }
		
		collectGarbage();
		
		if (saveDendros) 
		{
			cat('   ...saving dendrogram(s) as .jpgs\n');
			for (block in 1:length(net$blockGenes)) 
			{
				jpeg(file = paste(saveFileBase, 'run', run, 'dendro-block', block, '.jpg', sep=''), 
					 width = 15, 
					 height = 6, 
					 units = 'in', 
					 quality = 100, 
					 type = graphicsType, 
					 res = res
					 );
				plotDendroAndColors(net$dendrograms[[block]],
					net$colors[ net$blockGenes[[block]] ],
					groupLabels = 'module',
					rowText = net$colors[ net$blockGenes[[block]] ],
					main = paste('block ', block, ': ', length(net$blockGenes[[block]]), ' genes', sep=''),
					dendroLabels = F,
					hang = 0.05,
					addGuide = T,
					guideHang = 0.05
					);
				dev.off();
			}
		}
		
		if (densityPermTest) 
		{
			if (saveTOMs) 
			{
				cat('   ...loading TOM\n');
				load('tmpTOM-block.1.RData');
				TOM = as.matrix(TOM);
			} else if ( !(is.null(preTOM)) ) {
				cat('   ...using preTOM\n');
				TOM = preTOM;			
			} else {
				cat('   ...re-computing TOM\n');
				TOM = TOMsimilarityFromExpr(DATA, networkType = networkType, power = power);
			}
			collectGarbage();
			
			modDensities = getModDensitiesFromTOM(TOM, net$colors, plot = plotModDensities, skipGrey = skipGrey);
			cat('\n   ...working on density perm test\n');
			permTest = modDensityPerm(TOM, net$colors, modDensities, nPerm = nPerm, skipThresh = skipThresh, skipGrey = skipGrey);
			weakTOM = names( table(net$colors)[ permTest$pvals > permTestPvalThresh ] );
			cat('      ...weakTOM: ', weakTOM, '\n', sep='');
			bgGenes = net$colors %in% weakTOM;
			cat('   ......removing ', sum(bgGenes), ' background genes\n', sep='');
			DATA = DATA[, !bgGenes];
			
			if (savePermOut) 
			{
				permOut = list(modDensities=modDensities, permTest=permTest, weakTOM=weakTOM);
				save(permOut, file = paste(saveFileBase, 'run', run, 'PermOut.RData', sep = ''));
			}
			if ( length(weakTOM)==1 ) 
			{
				if (weakTOM=='grey') 
				{
					onlyGrey = onlyGrey + 1;
					cat('      onlyGrey = ', onlyGrey, '\n', sep='');
				}
			} else {
				cat('      ....re-setting onlyGrey...\n');
				onlyGrey = 0;
				cat('      onlyGrey = ', onlyGrey, '\n', sep='');
			}
			collectGarbage();
		} else {
			bgGenes = net$colors == 'grey';
			cat('\n   skipping permutation test...\n');
			cat('    ...removing ', sum(bgGenes), ' grey module genes\n', sep='');
			DATA = DATA[, !bgGenes];
			
			
			
			
			collectGarbage();
		}
		
		if (sum(bgGenes) == 0) 
		{
			cat('\nno more background genes, all done\n');
		}
	}						 	
							 						 	
	return(net);						 	
}

####
####
####


getModDensitiesFromTOM = function( TOM, colors, plot = T, diag = F, skipGrey = T ) {
	
	modTable = table(colors);
	nMods = length(modTable);
	modDensities = c();
	
	for( m in 1:nMods ) {
		if( skipGrey ) {
			if( names(modTable)[m] == 'grey' ) {
				modDensities = c(modDensities, 0);
				next;
			}
		}
		modGenes = colors == names(modTable)[m];	#print(modGenes);
		modTOM = TOM[modGenes, modGenes]; 			#print(is.matrix(modTOM));
		#cat( names(modTable)[m], '\n', dim(modTOM), '\n');
		modDensities = c(modDensities, mean( vectorizeMatrix(modTOM, diag = diag) ) );
		}
	names(modDensities) = names(modTable);
	
	if( plot ) {
		plot( as.vector(modTable), modDensities, type = 'n',
			  xlab = 'module size',
			  ylab = 'module density' 
			  );
		text( as.vector(modTable), modDensities, labels = names(modDensities) );
		}
	
	collectGarbage();
	return(modDensities);
	
	}
	
#####

modDensityPerm = function( TOM, colors, modDensities, nPerm = 1000, diag = F, skipThresh = 5000, skipGrey = T ) {
	
	modTable = table(colors);
	nMods = length(modTable);
	
	test = sum( names(modTable) == names(modDensities) ) == nMods;
	if( !test ) { error('module names don\'t match colors') }
	
	modPseudoDensities = list();
	pvals = c();
	
	for( m in 1:nMods ) {
		mod = names(modTable)[m];
		
		if( modTable[m] > skipThresh ) {
			if( mod == 'grey' ) {
				cat('skipping ', mod, ': ', modTable[m], ' probes\n', sep = '');
				modPseudoDensities[[m]] = rep(NA, nPerm);
				pvals = c(pvals, 1);
				cat(' pval = ', 1, '\n', sep = '');
				#print(mod); print(pvals); print(modPseudoDensities);
				next;
			} else {
				cat('skipping ', mod, ': ', modTable[m], ' probes\n', sep = '');
				modPseudoDensities[[m]] = rep(NA, nPerm);
				pvals = c(pvals, 0);
				cat(' pval = ', 0, '\n', sep = '');
				#print(mod); print(pvals); print(modPseudoDensities);
				next;
			}
		}
		if( skipGrey ) {
			if( mod == 'grey' ) {
				cat('skipping ', mod, ': ', modTable[m], ' probes\n', sep = '');
				modPseudoDensities[[m]] = rep(NA, nPerm);
				pvals = c(pvals, 1);
				cat(' pval = ', 1, '\n', sep = '');
				#print(mod); print(pvals); print(modPseudoDensities);
				next;
			}
		}
		
		cat('working on ', mod, ': ', modTable[m], ' probes\n', sep = '');
		pseudoDensities = c();
		for( p in 1:nPerm ) {
			if( p %% 1000 == 0 ) { cat('  permutation ', p, '\n', sep = '') }
			pseudoGenes = sample( 1:ncol(TOM), modTable[m] );	#print(pseudoGenes);
			pseudoTOM = TOM[pseudoGenes, pseudoGenes];			#print(pseudoTOM);
			pseudoDensities = c(pseudoDensities, mean( vectorizeMatrix(pseudoTOM, diag = diag) ));
		}
		
		modPseudoDensities[[m]] = pseudoDensities;
		p = ( sum( pseudoDensities > modDensities[m] ) ) / nPerm;
		pvals = c(pvals, p);
		cat(' pval = ', p, '\n', sep = '');
		collectGarbage();
		
		#print(mod); print(pvals); print(modPseudoDensities);
		
	}
	names(modPseudoDensities) = names(modTable);
	names(pvals) = names(modTable);
	
	permNA = ( 1:length(pvals) )[ is.na(pvals) ];
	if (length(permNA) > 0) {
		for ( p in 1:length(permNA) ) {
			permNAind = which( is.na( modPseudoDensities[[ permNA[p] ]] ) );
			permMed = median( modPseudoDensities[[ permNA[p] ]], na.rm=T );
			modPseudoDensities[[ permNA[p] ]][permNAind] = permMed;
			greater = sum( modPseudoDensities[[permNA[p]]] > modDensities[permNA[p]] );
			pvals[permNA[p]] = greater / nPerm;
		}
	}
	
	out = list(pvals = pvals, modPseudoDensities = modPseudoDensities);
	return(out);
	
}
	
#####

# # # modDensityPermSkipGrey = function( TOM, colors, modDensities, nPerm = 1000, diag = F, verbose = F ) {
	
	# # # modTable = table(colors);
	# # # nMods = length(modTable);
	
	# # # test = sum( names(modTable) == names(modDensities) ) == nMods;
	# # # if( !test ) { error('module names don\'t match colors')}
	
	# # # modPseudoDensities = list();
	# # # pvals = c();
	
	# # # for( m in 1:nMods ) {
		# # # mod = names(modTable)[m];
		# # # cat('working on ', mod, ': ', modTable[m], ' probes\n', sep = '');
		# # # pseudoDensities = c();
		
		# # # if( mod == 'grey' ) {
			# # # modPseudoDensities[[m]] = NULL;
			# # # pvals = c(pvals, 1);
			# # # } else {
				# # # for( p in 1:nPerm ) {
					# # # if( verbose ) {
						# # # if( p %% 500 == 0 ) { cat('  permutation ', p, '\n', sep = '') }
					# # # }
					# # # pseudoGenes = sample( 1:ncol(TOM), modTable[m] );
					# # # pseudoTOM = TOM[pseudoGenes, pseudoGenes];
					# # # pseudoDensities = c(pseudoDensities, mean( vectorizeMatrix(pseudoTOM, diag = diag) ));
					# # # }
					
				# # # modPseudoDensities[[m]] = pseudoDensities;
				# # # p = ( sum( pseudoDensities > modDensities[m] ) ) / nPerm;
				# # # pvals = c(pvals, p);
				# # # collectGarbage();
				# # # }
		
		# # # }
	# # # names(modPseudoDensities) = names(modTable);
	
	# # # permNA = ( 1:length(pvals) )[ is.na(pvals) ];
	
	# # # if (length(permNA) > 0) {
			# # # for ( p in 1:length(permNA) ) {
				# # # permNAind = which( is.na( modPseudoDensities[[ permNA[p] ]] ) );
				# # # permMed = median( modPseudoDensities[[ permNA[p] ]], na.rm=T );
				# # # modPseudoDensities[[ permNA[p] ]][permNAind] = permMed;
				# # # greater = sum( modPseudoDensities[[permNA[p]]] > modDensities[permNA[p]] );
				# # # pvals[permNA[p]] = greater / nPerm;
			# # # }
		# # # }
	
	# # # out = list(pvals = pvals, modPseudoDensities = modPseudoDensities);
	# # # return(out);
	
	# # # }