rm(list=ls());
library(Rsubread);

#Bonobo
basename='Bonobo';
setwd("~/Documents/Rsubread/Bonobo");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Bonobo/Bonobo.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	
		 
#Chicken
basename='Chicken';
setwd("~/Documents/Rsubread/Chicken");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Chicken/Chicken.fasta',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Chimp
basename='Chimp';
setwd("~/Documents/Rsubread/Chimp");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Chimp/Chimp.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Gorilla
Basename='Gorilla';
setwd("~/Documents/Rsubread/Gorilla");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Gorilla/Gorilla.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Human
basename='Human';
setwd("~/Documents/Rsubread/Human");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Human/hg19.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Macaque
basename='Macaque';
setwd("~/Documents/Rsubread/Macaque");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Macaque/Macaque.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Mouse
basename='Mouse';
setwd("~/Documents/Rsubread/Mouse");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Mouse/Mouse.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Opossum
basename='Opossum';
setwd("~/Documents/Rsubread/Opossum");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Opossum/Opossum.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Orangutan
basename='Orangutan';
setwd("~/Documents/Rsubread/Orangutan");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Orangutan/Orangutan.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	

#Platypus
basename='Platypus';
setwd("~/Documents/Rsubread/Platypus");
buildindex(basename=basename,
		 reference='/Volumes/handsfiler$/FishStudies/_mammalian_RNAseq/Platypus/Platypus.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	