> blockwiseModules
function (datExpr, blocks = NULL, maxBlockSize = 5000, randomSeed = 12345, 
    corType = "pearson", power = 6, networkType = "unsigned", 
    TOMType = "signed", TOMDenom = "min", deepSplit = 2, detectCutHeight = 0.995, 
    minModuleSize = min(20, ncol(datExpr)/2), maxCoreScatter = NULL, 
    minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, 
    pamStage = TRUE, pamRespectsDendro = TRUE, minCoreKME = 0.5, 
    minCoreKMESize = minModuleSize/3, minKMEtoStay = 0.3, reassignThreshold = 1e-06, 
    mergeCutHeight = 0.15, impute = TRUE, getTOMs = NULL, saveTOMs = FALSE, 
    saveTOMFileBase = "blockwiseTOM", trapErrors = FALSE, numericLabels = FALSE, 
    checkMissingData = TRUE, maxPOutliers = 1, quickCor = 0, 
    pearsonFallback = "individual", cosineCorrelation = FALSE, 
    nThreads = 0, verbose = 0, indent = 0) 
{
    spaces = indentSpaces(indent)
    if (verbose > 0) 
        printFlush(paste(spaces, "Calculating module eigengenes block-wise from all genes"))
    seedSaved = FALSE
    if (!is.null(randomSeed)) {
        if (exists(".Random.seed")) {
            seedSaved = TRUE
            savedSeed = .Random.seed
        }
        set.seed(randomSeed)
    }
    intCorType = pmatch(corType, .corTypes)
    if (is.na(intCorType)) 
        stop(paste("Invalid 'corType'. Recognized values are", 
            paste(.corTypes, collapse = ", ")))
    intTOMType = pmatch(TOMType, .TOMTypes)
    if (is.na(intTOMType)) 
        stop(paste("Invalid 'TOMType'. Recognized values are", 
            paste(.TOMTypes, collapse = ", ")))
    TOMDenomC = pmatch(TOMDenom, .TOMDenoms) - 1
    if (is.na(TOMDenomC)) 
        stop(paste("Invalid 'TOMDenom'. Recognized values are", 
            paste(.TOMDenoms, collapse = ", ")))
    if ((maxPOutliers < 0) | (maxPOutliers > 1)) 
        stop("maxPOutliers must be between 0 and 1.")
    if (quickCor < 0) 
        stop("quickCor must be positive.")
    if (nThreads < 0) 
        stop("nThreads must be positive.")
    if (is.null(nThreads) || (nThreads == 0)) 
        nThreads = .useNThreads()
    if ((power < 1) | (power > 30)) 
        stop("power must be between 1 and 30.")
    intNetworkType = charmatch(networkType, .networkTypes)
    if (is.na(intNetworkType)) 
        stop(paste("Unrecognized networkType argument.", "Recognized values are (unique abbreviations of)", 
            paste(.networkTypes, collapse = ", ")))
    fallback = pmatch(pearsonFallback, .pearsonFallbacks)
    if (is.na(fallback)) 
        stop(spaste("Unrecognized value '", pearsonFallback, 
            "' of argument 'pearsonFallback'.", "Recognized values are (unique abbreviations of)\n", 
            paste(.pearsonFallbacks, collapse = ", ")))
    dimEx = dim(datExpr)
    if (length(dimEx) != 2) 
        stop("datExpr has incorrect dimensions.")
    nGenes = dimEx[2]
    nSamples = dimEx[1]
    allLabels = rep(0, nGenes)
    AllMEs = NULL
    allLabelIndex = NULL
    if (!is.null(blocks) && (length(blocks) != nGenes)) 
        stop("Input error: the length of 'geneRank' does not equal the number of genes in given 'datExpr'.")
    if (!is.null(getTOMs)) 
        warning("getTOMs is deprecated, please use saveTOMs instead.")
    if (checkMissingData) {
        gsg = goodSamplesGenes(datExpr, verbose = verbose - 1, 
            indent = indent + 1)
        if (!gsg$allOK) 
            datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
        nGGenes = sum(gsg$goodGenes)
        nGSamples = sum(gsg$goodSamples)
    }
    else {
        nGGenes = nGenes
        nGSamples = nSamples
        gsg = list(goodSamples = rep(TRUE, nSamples), goodGenes = rep(TRUE, 
            nGenes), allOK = TRUE)
    }
    if (is.null(blocks)) {
        if (nGGenes > maxBlockSize) {
            if (verbose > 1) 
                printFlush(paste(spaces, "....pre-clustering genes to determine blocks.."))
            clustering = projectiveKMeans(datExpr, preferredSize = maxBlockSize, 
                checkData = FALSE, sizePenaltyPower = 5, verbose = verbose - 
                  2, indent = indent + 1)
            gBlocks = clustering$clusters
            if (verbose > 2) {
                printFlush("Block sizes:")
                print(table(gBlocks))
            }
        }
        else gBlocks = rep(1, nGGenes)
        blocks = rep(NA, nGenes)
        blocks[gsg$goodGenes] = gBlocks
    }
    else {
        gBlocks = blocks[gsg$goodGenes]
    }
    blockLevels = as.numeric(levels(factor(gBlocks)))
    blockSizes = table(gBlocks)
    blockOrder = order(-blockSizes)
    nBlocks = length(blockLevels)
    dendros = list()
    TOMFiles = rep("", nBlocks)
    blockGenes = list()
    blockNo = 1
    maxUsedLabel = 0
    while (blockNo <= nBlocks) {
        if (verbose > 1) 
            printFlush(paste(spaces, "..Working on block", blockNo, 
                "."))
        blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks == 
            blockLevels[blockOrder[blockNo]]]
        block = c(1:nGGenes)[gBlocks == blockLevels[blockOrder[blockNo]]]
        selExpr = as.matrix(datExpr[, block])
        nBlockGenes = length(block)
        dissTom = matrix(0, nBlockGenes, nBlockGenes)
        callVerb = max(0, verbose - 1)
        callInd = indent + 2
        CcorType = intCorType - 1
        CnetworkType = intNetworkType - 1
        CTOMType = intTOMType - 1
        warn = 0
        tomResult = .C("tomSimilarity", as.double(selExpr), as.integer(nSamples), 
            as.integer(nBlockGenes), as.integer(CcorType), as.integer(CnetworkType), 
            as.double(power), as.integer(CTOMType), as.integer(TOMDenomC), 
            as.double(maxPOutliers), as.double(quickCor), as.integer(fallback), 
            as.integer(cosineCorrelation), tom = as.double(dissTom), 
            warn = as.integer(warn), as.integer(nThreads), as.integer(callVerb), 
            as.integer(callInd), NAOK = TRUE, DUP = FALSE, PACKAGE = "WGCNA")
        dim(tomResult$tom) = c(nBlockGenes, nBlockGenes)
        if (saveTOMs) {
            TOM = as.dist(tomResult$tom)
            TOMFiles[blockNo] = paste(saveTOMFileBase, "-block.", 
                blockNo, ".RData", sep = "")
            if (verbose > 2) 
                printFlush(paste(spaces, "  ..saving TOM for block", 
                  blockNo, "into file", TOMFiles[blockNo]))
            save(TOM, file = TOMFiles[blockNo])
            rm(TOM)
            collectGarbage()
        }
        dissTom = 1 - tomResult$tom
        dim(dissTom) = c(nBlockGenes, nBlockGenes)
        rm(tomResult)
        collectGarbage()
        if (verbose > 2) 
            printFlush(paste(spaces, "....clustering.."))
        dendros[[blockNo]] = flashClust(as.dist(dissTom), method = "average")
        if (verbose > 2) 
            printFlush(paste(spaces, "....detecting modules.."))
        blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]], 
            deepSplit = deepSplit, cutHeight = detectCutHeight, 
            minClusterSize = minModuleSize, method = "hybrid", 
            maxCoreScatter = maxCoreScatter, minGap = minGap, 
            maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap, 
            pamStage = pamStage, pamRespectsDendro = pamRespectsDendro, 
            distM = dissTom, verbose = verbose - 3, indent = indent + 
                2), silent = TRUE)
        collectGarbage()
        if (verbose > 8) {
            if (interactive()) 
                plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), 
                  dendroLabels = FALSE, main = paste("Block", 
                    blockNo))
        }
        if (class(blockLabels) == "try-error") {
            if (verbose > 0) {
                printFlush(paste(spaces, "*** cutreeDynamic returned the following error:\n", 
                  spaces, blockLabels, spaces, "Stopping the module detection here."))
            }
            else warning(paste("blockwiseModules: cutreeDynamic returned the following error:\n", 
                "      ", blockLabels, "---> Continuing with next block. "))
            next
        }
        if (sum(blockLabels > 0) == 0) {
            if (verbose > 1) {
                printFlush(paste(spaces, "No modules detected in block", 
                  blockNo))
            }
            blockNo = blockNo + 1
            next
        }
        blockLabels[blockLabels > 0] = blockLabels[blockLabels > 
            0] + maxUsedLabel
        maxUsedLabel = max(blockLabels)
        if (verbose > 2) 
            printFlush(paste(spaces, "....calculating module eigengenes.."))
        MEs = try(moduleEigengenes(selExpr[, blockLabels != 0], 
            blockLabels[blockLabels != 0], impute = impute, verbose = verbose - 
                3, indent = indent + 2), silent = TRUE)
        if (class(MEs) == "try-error") {
            if (trapErrors) {
                if (verbose > 0) {
                  printFlush(paste(spaces, "*** moduleEigengenes failed with the following message:"))
                  printFlush(paste(spaces, "       ", MEs))
                  printFlush(paste(spaces, "    ---> Stopping module detection here."))
                }
                else warning(paste("blockwiseModules: moduleEigengenes failed with the following message:", 
                  "\n     ", MEs, "---> Continuing with next block. "))
                next
            }
            else stop(MEs)
        }
        propMEs = MEs$eigengenes
        blockLabelIndex = as.numeric(substring(names(propMEs), 
            3))
        deleteModules = NULL
        changedModules = NULL
        if (verbose > 2) 
            printFlush(paste(spaces, "....checking modules for statistical meaningfulness.."))
        for (mod in 1:ncol(propMEs)) {
            modGenes = (blockLabels == blockLabelIndex[mod])
            corEval = parse(text = paste(.corFnc[intCorType], 
                "(selExpr[, modGenes], propMEs[, mod]", prepComma(.corOptions[intCorType]), 
                ")"))
            KME = as.vector(eval(corEval))
            if (intNetworkType == 1) 
                KME = abs(KME)
            if (sum(KME > minCoreKME) < minCoreKMESize) {
                blockLabels[modGenes] = 0
                deleteModules = c(deleteModules, mod)
                if (verbose > 3) 
                  printFlush(paste(spaces, "    ..deleting module ", 
                    mod, ": of ", sum(modGenes), " total genes in the module\n       only ", 
                    sum(KME > minCoreKME), " have the requisite high correlation with the eigengene.", 
                    sep = ""))
            }
            else if (sum(KME < minKMEtoStay) > 0) {
                if (verbose > 2) 
                  printFlush(paste(spaces, "    ..removing", 
                    sum(KME < minKMEtoStay), "genes from module", 
                    mod, "because their KME is too low."))
                blockLabels[modGenes][KME < minKMEtoStay] = 0
                if (sum(blockLabels[modGenes] > 0) < minModuleSize) {
                  deleteModules = c(deleteModules, mod)
                  blockLabels[modGenes] = 0
                  if (verbose > 3) 
                    printFlush(paste(spaces, "    ..deleting module ", 
                      blockLabelIndex[mod], ": not enough genes in the module after removal of low KME genes.", 
                      sep = ""))
                }
                else {
                  changedModules = union(changedModules, blockLabelIndex[mod])
                }
            }
        }
        if (!is.null(deleteModules)) {
            propMEs = propMEs[, -deleteModules, drop = FALSE]
            modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
            blockLabels[modGenes] = 0
            modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
            allLabels[modAllGenes] = 0
            blockLabelIndex = blockLabelIndex[-deleteModules]
        }
        if (sum(blockLabels > 0) == 0) {
            if (verbose > 1) {
                printFlush(paste(spaces, "No significant modules detected in block", 
                  blockNo))
            }
            blockNo = blockNo + 1
            next
        }
        if (is.null(AllMEs)) {
            AllMEs = propMEs
        }
        else AllMEs = cbind(AllMEs, propMEs)
        allLabelIndex = c(allLabelIndex, blockLabelIndex)
        assigned = block[blockLabels != 0]
        allLabels[gsg$goodGenes][assigned] = blockLabels[blockLabels != 
            0]
        rm(dissTom)
        collectGarbage()
        blockNo = blockNo + 1
    }
    deleteModules = NULL
    goodLabels = allLabels[gsg$goodGenes]
    if (sum(goodLabels != 0) > 0) {
        propLabels = goodLabels[goodLabels != 0]
        assGenes = c(1:nGenes)[gsg$goodGenes[goodLabels != 0]]
        corEval = parse(text = paste(.corFnc[intCorType], "(datExpr[, goodLabels!=0], AllMEs", 
            prepComma(.corOptions[intCorType]), ")"))
        KME = eval(corEval)
        if (intNetworkType == 1) 
            KME = abs(KME)
        nMods = ncol(AllMEs)
        for (mod in 1:nMods) {
            modGenes = c(1:length(propLabels))[propLabels == 
                allLabelIndex[mod]]
            KMEmodule = KME[modGenes, mod]
            KMEbest = apply(KME[modGenes, , drop = FALSE], 1, 
                max)
            candidates = (KMEmodule < KMEbest)
            candidates[!is.finite(candidates)] = FALSE
            if (sum(candidates) > 0) {
                pModule = corPvalueFisher(KMEmodule[candidates], 
                  nSamples)
                whichBest = apply(KME[modGenes[candidates], , 
                  drop = FALSE], 1, which.max)
                pBest = corPvalueFisher(KMEbest[candidates], 
                  nSamples)
                reassign = ifelse(is.finite(pBest/pModule), (pBest/pModule < 
                  reassignThreshold), FALSE)
                if (sum(reassign) > 0) {
                  if (verbose > 2) 
                    printFlush(paste(spaces, " ..reassigning", 
                      sum(reassign), "genes from module", mod, 
                      "to modules with higher KME."))
                  allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
                  changedModules = union(changedModules, whichBest[reassign])
                  if (sum(modGenes) - sum(reassign) < minModuleSize) {
                    deleteModules = c(deleteModules, mod)
                  }
                  else changedModules = union(changedModules, 
                    mod)
                }
            }
        }
    }
    if (!is.null(deleteModules)) {
        AllMEs = AllMEs[, -deleteModules, drop = FALSE]
        genes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
        allLabels[genes] = 0
        allLabelIndex = allLabelIndex[-deleteModules]
        goodLabels = allLabels[gsg$goodGenes]
    }
    if (verbose > 1) 
        printFlush(paste(spaces, "..merging modules that are too close.."))
    if (numericLabels) {
        colors = allLabels
    }
    else {
        colors = labels2colors(allLabels)
    }
    mergedAllColors = colors
    MEsOK = TRUE
    mergedMods = try(mergeCloseModules(datExpr, colors[gsg$goodGenes], 
        cutHeight = mergeCutHeight, relabel = TRUE, impute = impute, 
        verbose = verbose - 2, indent = indent + 2), silent = TRUE)
    if (class(mergedMods) == "try-error") {
        warning(paste("blockwiseModules: mergeCloseModules failed with the following error message:\n    ", 
            mergedMods, "\n--> returning unmerged colors.\n"))
        MEs = try(moduleEigengenes(datExpr, colors[gsg$goodGenes], 
            impute = impute, verbose = verbose - 3, indent = indent + 
                3), silent = TRUE)
        if (class(MEs) == "try-error") {
            if (!trapErrors) 
                stop(MEs)
            if (verbose > 0) {
                printFlush(paste(spaces, "*** moduleEigengenes failed with the following error message:"))
                printFlush(paste(spaces, "     ", MEs))
                printFlush(paste(spaces, "*** returning no module eigengenes.\n"))
            }
            else warning(paste("blockwiseModules: moduleEigengenes failed with the following error message:\n    ", 
                MEs, "\n--> returning no module eigengenes.\n"))
            allSampleMEs = NULL
            MEsOK = FALSE
        }
        else {
            if (sum(!MEs$validMEs) > 0) {
                colors[gsg$goodGenes] = MEs$validColors
                MEs = MEs$eigengenes[, MEs$validMEs]
            }
            else MEs = MEs$eigengenes
            allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, 
                ncol = ncol(MEs)))
            allSampleMEs[gsg$goodSamples, ] = MEs[, ]
            names(allSampleMEs) = names(MEs)
        }
    }
    else {
        mergedAllColors[gsg$goodGenes] = mergedMods$colors
        allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, 
            ncol = ncol(mergedMods$newMEs)))
        allSampleMEs[gsg$goodSamples, ] = mergedMods$newMEs[, 
            ]
        names(allSampleMEs) = names(mergedMods$newMEs)
    }
    if (seedSaved) 
        .Random.seed <<- savedSeed
    if (!saveTOMs) 
        TOMFiles = NULL
    list(colors = mergedAllColors, unmergedColors = colors, MEs = allSampleMEs, 
        goodSamples = gsg$goodSamples, goodGenes = gsg$goodGenes, 
        dendrograms = dendros, TOMFiles = TOMFiles, blockGenes = blockGenes, 
        blocks = blocks, blockOrder = blockLevels[blockOrder], 
        MEsOK = MEsOK)
}
<environment: namespace:WGCNA>
> 