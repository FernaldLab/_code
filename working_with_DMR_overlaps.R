rm(list=ls());
options(stringsAsFactors=F);
setwd('~/Documents/_BSeq_data/dmr5x_n5_md.15');
#############################################################################################
##### UTRs
#############################################################################################
# read in DMR-UTR results
utr = read.table("Astatotilapia_burtoni.BROADAB1.UTRs.gff3-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
utr = utr[utr$V13<5000 & utr$V13>-5000, ];

# get rid of lines where there was no hit
utr = utr[utr$V6!=".",];

# split into list where each element is a DMR
utrs = split(utr, paste(utr[,1],":",utr[,2],"-",utr[,3],sep=""));

# get number of hits for each DMR
utrsnrow = sapply(utrs, nrow);

# if DMR has >1 hit its most likely the same interval but for different transcript isoforms
# check that all dists are the same
for (i in utrs) { if(length(unique(i$V13))>1){print(i)} }; rm(i)

# check if any DMRs are close to >1 gene
utrsngene = c();
for (i in utrs) {
     tmp = strsplit(grep("mrna", unlist(strsplit(as.character(i$V12), "Parent=")), value=T), ".", fixed=T);
     genes = c();
     for (j in tmp) { genes = c(genes, paste(j[3], j[4], sep=".")) };
     utrsngene = c(utrsngene, length(unique(genes)));
} 
rm(i,tmp,genes,j);

# assuming check was true
# create vec of overlap distances by taking first value from each list component
###utrdists = c();
###for (i in utrs) { utrdists = c(utrdists, i$V13[1]) }; rm(i);

# make matrix of first lines from each list element
###utrred = utrs[[1]][1, ];
###for (i in 2:length(utrs)) { utrred = rbind(utrred, utrs[[i]][1, ]) }; rm(i);

#############################################################################################
##### CDSs
#############################################################################################
cds = read.table("Astatotilapia_burtoni.BROADAB2fix_CDS.gtf-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
cds = cds[cds$V13<5000 & cds$V13>-5000, ];

# get rid of lines where there was no hit
cds = cds[cds$V6!=".",];

# split into list where each element is a DMR
cdss = split(cds, paste(cds[,1],":",cds[,2],"-",cds[,3],sep=""));

# get number of hits for each DMR
cdssnrow = sapply(cdss, nrow);

# if DMR has >1 hit its most likely the same interval but for different transcript isoforms
# check that all dists are the same
for (i in cdss) { if(length(unique(i$V13))>1){print(i)} }; rm(i)

# check if any DMRs are close to >1 gene
cdssngene = c();
for (i in cdss) {
     tmp = grep("gene", unlist(strsplit(as.character(i$V12), ";")), value=T);
     cdssngene = c(cdssngene, length(unique(tmp)));
} 
rm(i,tmp);

# assuming check was true
# create vec of overlap distances by taking first value from each list component
###cdsdists = c();
###for (i in cdss) { cdsdists = c(cdsdists, i$V13[1]) }; rm(i);

# make matrix of first lines from each list element
###cdsred = cdss[[1]][1, ];
###for (i in 2:length(cdss)) { cdsred = rbind(cdsred, cdss[[i]][1, ]) }; rm(i);

#############################################################################################
##### introns
#############################################################################################
int = read.table("Astatotilapia_burtoni.BROADAB2fix.intron.gtf-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
int = int[int$V13<5000 & int$V13>-5000, ];

# get rid of lines where there was no hit
int = int[int$V6!=".",];

# split into list where each element is a DMR
ints = split(int, paste(int[,1],":",int[,2],"-",int[,3],sep=""));

# get number of hits for each DMR
intsnrow = sapply(ints, nrow);

# if DMR has >1 hit its most likely the same interval but for different transcript isoforms
# check that all dists are the same
for (i in ints) { if(length(unique(i$V13))>1){print(i)} }; rm(i)

# check if any DMRs are close to >1 gene
intsngene = c();
for (i in ints) {
     tmp = grep("gene", unlist(strsplit(as.character(i$V12), ";")), value=T);
     intsngene = c(intsngene, length(unique(tmp)));
} 
rm(i,tmp);

# assuming check was true
# create vec of overlap distances by taking first value from each list component
###intdists = c();
###for (i in ints) { intdists = c(intdists, i$V13[1]) }; rm(i);

# make matrix of first lines from each list element
###intred = ints[[1]][1, ];
###for (i in 2:length(ints)) { intred = rbind(intred, ints[[i]][1, ]) }; rm(i);

#############################################################################################
##### lncRNAs
#############################################################################################
lnc = read.table("abur.lnc.final.gtf-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
lnc = lnc[lnc$V13<5000 & lnc$V13>-5000, ];

# get rid of lines where there was no hit
lnc = lnc[lnc$V6!=".",];

# split into list where each element is a DMR
lncs = split(lnc, paste(lnc[,1],":",lnc[,2],"-",lnc[,3],sep=""));

# only 6 DMR hits for lncRNAs

#############################################################################################
##### miRNAs
#############################################################################################
mir = read.table("abur_miRNAs-130326.fix.bed-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
mir = mir[mir$V13<5000 & mir$V13>-5000, ];

# get rid of lines where there was no hit
mir = mir[mir$V4!=".",];			# note different column

# split into list where each element is a DMR
mirs = split(mir, paste(mir[,1],":",mir[,2],"-",mir[,3],sep=""));

# only 2 DMR hits for miRNAs

#############################################################################################
##### transposons
#############################################################################################
te = read.table("Abur_final_TE.bed-closest_Dref", sep="\t", header=F);

# keep only hits within 5kb
te = te[te$V10<5000 & te$V10>-5000, ];		# note different column

# get rid of lines where there was no hit
te = te[te$V7!=".",];				# note different column

# split into list where each element is a DMR
tes = split(te, paste(te[,1],":",te[,2],"-",te[,3],sep=""));

# get number of hits for each DMR
tesnrow = sapply(tes, nrow);

# some DMRs have >1 te hit, all are unique
tes[tesnrow>1]

#############################################################################################
##### compare utrs, cdss, ints, lncs, mirs, tes
#############################################################################################
library(VennDiagram);library(gridExtra);
###utr_cds = intersect(names(utrs), names(cdss));
###utr_cds_int = intersect(utr_cds, names(ints));

ov = list(utrs = utrs[sapply(utrs, function(f) all(f$V13==0))],
		  cdss = cdss[sapply(cdss, function(f) all(f$V13==0))],
		  ints = ints[sapply(ints, function(f) all(f$V13==0))],
		  lncs = lncs[sapply(lncs, function(f) all(f$V13==0))],
		  mirs = mirs[sapply(mirs, function(f) all(f$V13==0))],
		  tes = tes[sapply(tes, function(f) all(f$V10==0))]
		  );
utrcol = 'blue';
cdscol = 'grey';
intcol = 'yellow';
lnccol = 'red';
tecol = 'green';

ovutrtype = sapply(ov$utrs, function(f) unique(f$V6));

ovcdsnum = c();
for (i in ov$cdss) {
	tmp = unlist(strsplit(i$V12, '; '));
	if (any(tmp == 'exon_number 1')) {
		ovcdsnum = c(ovcdsnum, TRUE);
	} else {
		ovcdsnum = c(ovcdsnum, FALSE);
	}
}; rm(i);

ovintnum = c();
for (i in ov$ints) {
	tmp = unlist(strsplit(i$V12, '; '));
	if (any(tmp == 'intron_number 1;')) {
		ovintnum = c(ovintnum, TRUE);
	} else {
		ovintnum = c(ovintnum, FALSE);
	}
}; rm(i);


A1 = ov$utrs[ovutrtype=='five_prime_utr'];
A2 = ov$cdss[ovcdsnum];
COL = c(utrcol,cdscol);
CAT = c("5'UTR","CDS1");
P1=draw.pairwise.venn(area1=length(A1), area2=length(A2), 
				   cross.area=length(intersect(names(A1),names(A2))),
				   col=COL, fill=COL,
				   fontfamily=rep('sans',3), cat.fontfamily=rep('sans',2),
				   category=CAT
				   );
A1 = ov$utrs[ovutrtype=='five_prime_utr'];
A2 = ov$ints[ovintnum];
COL = c(utrcol,intcol);
CAT = c("5'UTR","INT1");
P2=draw.pairwise.venn(area1=length(A1), area2=length(A2), 
				   cross.area=length(intersect(names(A1),names(A2))),
				   col=COL, fill=COL,
				   fontfamily=rep('sans',3), cat.fontfamily=rep('sans',2),
				   category=CAT
				   );
A1 = ov$cdss[ovcdsnum];
A2 = ov$ints[ovintnum];
COL = c(cdscol,intcol);
CAT = c("CDS1","INT1");
P3=draw.pairwise.venn(area1=length(A1), area2=length(A2), 
				   cross.area=length(intersect(names(A1),names(A2))),
				   col=COL, fill=COL,
				   fontfamily=rep('sans',3), cat.fontfamily=rep('sans',2),
				   category=CAT
				   );
				   


A1 = ov$ints[ovintnum];
A2 = ov$utrs[ovutrtype=='five_prime_utr'];
A3 = ov$cdss[ovcdsnum];			   
COL = c(intcol,utrcol,cdscol);
CAT = c('INT1','5UTR','CDS1');
P4=draw.triple.venn(area1=length(A1), area2=length(A2), area3=length(A3),
				 n12=length(intersect(names(A1), names(A2))), 
				 n23=length(intersect(names(A2), names(A3))), 
				 n13=length(intersect(names(A1), names(A3))),
				 n123=length(intersect(intersect(names(A1), names(A2)), names(A3))),
				 category=CAT,
				 fill=cols, col=cols,
				 fontfamily=rep('sans',7), cat.fontfamily=rep('sans',3),
				# cat.pos=c(-40,40,0),
				# alpha=rep(.7,3), 
				# euler.d=F, scaled=F
				 );
grid.arrange(gTree(children=P1), gTree(children=P2), gTree(children=P3), gTree(children=P4), ncol=2, nrow=2)


alln = sort(unique(c(names(utrs), names(cdss), names(ints), names(lncs), names(mirs), names(tes))));

dmrs = list();
for (d in alln) {
     tmplist = list();
     tmp = as.numeric(unlist(strsplit( unlist(strsplit(d, ":"))[2], "-")));
     tmplist[["width"]] = tmp[2] - tmp[1] + 1;
     if (d %in% names(utrs)) { tmplist[["utrs"]] = utrs[names(utrs) == d][[1]] };
     if (d %in% names(cdss)) { tmplist[["cdss"]] = cdss[names(cdss) == d][[1]] };
     if (d %in% names(ints)) { tmplist[["ints"]] = ints[names(ints) == d][[1]] };
     if (d %in% names(lncs)) { tmplist[["lncs"]] = lncs[names(lncs) == d][[1]] };
     if (d %in% names(mirs)) { tmplist[["mirs"]] = mirs[names(mirs) == d][[1]] };
     if (d %in% names(tes)) { tmplist[["tes"]] = tes[names(tes) == d][[1]] };
     dmrs[[d]] = tmplist;
}; rm(d,tmplist,tmp);

dmrsOverlap = c();
for (d in dmrs) {
     for (i in 1:length(d)) {
     	  if (names(d)[i] == "width") { next };
          if (names(d)[i] == "tes") {  #print(d[[i]])
               if (all(as.numeric(d[[i]]$V10) == 0)) { dmrsOverlap = c(dmrsOverlap, d) };
          } else {#print(d[[i]])
               if (all(as.numeric(d[[i]]$V13) == 0)) { dmrsOverlap = c(dmrsOverlap, d) };
          }
     }
}; rm(d,i);
