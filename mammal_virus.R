primates.HIV.network.INTERSECT=primates_hiv1_interactions[,1]%in%colnames(DATA)
mammal.HIV.network.INTERSECT=mammals_hiv1_interactions[,1]%in%colnames(DATA)
mammal.PUBMED.network.INTERSECT=mammals_pubmed_interactions[,1]%in%colnames(DATA)

#Intersection of interaction data with mammals
temp1 = mammals_hiv1_interactions[mammal.HIV.network.INTERSECT,]
temp2 = mammals_pubmed_interactions[mammal.PUBMED.network.INTERSECT,]
temp3 = primates_hiv1_interactions[primates.HIV.network.INTERSECT,]

#Subset intersection for only genes with real interactions
mammal.HIV.interaction.GENES=temp1[temp1[,2]==1,1]
mammal.PUBMED.interaction.GENES=temp2[temp2[,2]==1,1]
primate.HIV.interaction.GENES=temp3[temp3[,2]==1,1]

#All mammal/primate interaction ids in network
mammal.HIV.module.genes=table(colors[colnames(DATA)%in%temp1[,1]])
#Only real interactions in network
mammal.HIV.module.INTERACTION.genes=table(colors[colnames(DATA)%in%temp1[temp1[,2]==1,1],1)

#Get gene names in module
mod.genes=exn.getModuleGenes(DATA, colors)

#Check for enrichments of interaction genes
mammal.HIV.enrichments.1=checkGeneListEnrichmentList(mammal.HIV.interaction.GENES,mod.genes,names(DATA))
mammal.PUBMED.enrichments.1 =checkGeneListEnrichmentList(mammal.PUBMED.interaction.GENES,mod.genes,names(DATA)) 
primate.HIV.enrichments.1 =checkGeneListEnrichmentList(primate.HIV.interaction.GENES,mod.genes,names(DATA)) 

#Module membership of enriched genes
lightgreen=exn.getModulekME("lightgreen", colors, kME)
darkred=exn.getModulekME("darkred", colors, kME)
black=exn.getModulekME("black", colors, kME)
darkgrey=exn.getModulekME("darkgrey", colors, kME)
green=exn.getModulekME("green", colors, kME)

#Find kME of interaction genes in module
interaction.genes.mammHIV=lightgreen[rownames(lightgreen)%in%temp1[temp1[,2]==1,1],]
interaction.genes.mammHIV=darkred[rownames(darkred)%in%temp1[temp1[,2]==1,1],]
interaction.genes.mammHIV=green[rownames(green)%in%temp2[temp2[,2]==1,1],]
interaction.genes.mammHIV=darkgrey[rownames(darkgrey)%in%temp1[temp1[,2]==1,1],]

#Subset module for only interaction genes
lightgreen.interaction.genes=rownames(lightgreen)%in%mammal.HIV.interaction.GENES
darkred.interaction.genes=rownames(darkred)%in%mammal.HIV.interaction.GENES
black.interaction.genes=rownames(black)%in%mammal.HIV.interaction.GENES
green.interaction.genes=rownames(green)%in%mammal.PUBMED.interaction.GENES
darkgrey.interaction.genes=rownames(darkgrey)%in%mammal.PUBMED.interaction.GENES

#Test interaction genes kME/p-value in module vs. non-interaction genes 
verboseBoxplot(as.numeric(lightgreen[,1]), as.factor(lightgreen.interaction.genes))
verboseBoxplot(as.numeric(darkred[,1]), as.factor(darkred.interaction.genes))
verboseBoxplot(as.numeric(black[,1]), as.factor(black.interaction.genes))
verboseBoxplot(as.numeric(green[,1]), as.factor(green.interaction.genes))
verboseBoxplot(as.numeric(darkgrey[,1]), as.factor(darkgrey.interaction.genes))