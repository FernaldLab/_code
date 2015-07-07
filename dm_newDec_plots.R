# figures in /Volumes/fishstudies/_methylation/_paper/Feb2015_doc_figs

tmp = datC$nucs;
tmp[datC$strand=='-'] = datC$nucsRC[datC$strand=='-'];
.seqLogoFromSeqVec(tmp, ic.scale=F)

.seqLogoFromSeqVec(tmp[datC$higher=='D'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$strand=='+'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$strand=='-'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$diH=='CG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$diH=='CH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$triHH=='CHG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$triHH=='CGH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$triHH=='CHH'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$strand=='+'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$strand=='-'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$strand=='+'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$strand=='-'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CGH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHH'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHG'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CGH'], ic.scale=F);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHH'], ic.scale=F);

.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CG' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CG' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$diH=='CH' & datC$strand=='-'], ic.scale=T);

.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CG' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CG' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$diH=='CH' & datC$strand=='-'], ic.scale=T);

.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHG' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHG' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CGH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CGH' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='D' & datC$triHH=='CHH' & datC$strand=='-'], ic.scale=T);

.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHG' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHG' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CGH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CGH' & datC$strand=='-'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHH' & datC$strand=='+'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$higher=='ND' & datC$triHH=='CHH' & datC$strand=='-'], ic.scale=T);

############

.seqLogoFromSeqVec(tmp[datC$strand=='+' & datC$diH=='CG'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$strand=='-' & datC$diH=='CG'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$strand=='+' & datC$diH=='CH'], ic.scale=T);
.seqLogoFromSeqVec(tmp[datC$strand=='-' & datC$diH=='CH'], ic.scale=T);







############

# figure 1
#
par(mfrow=c(3,3), oma=c(0,1,2,0));
CT = 1.1;
CL = 1.5;
LO = .2;
DM = datC;
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd+rev', ylab='all DM sites', cex.lab=CL);

DM = datC[datC$strand=='+', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd');

DM = datC[datC$strand=='-', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='rev');

DM = datC[datC$higher=='ND', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd+rev', ylab='ND higher', cex.lab=CL);

DM = DM[DM$strand=='+', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd');

DM = datC[datC$higher=='ND', ];
DM = DM[DM$strand=='-', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='rev');

DM = datC[datC$higher=='D', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd+rev', ylab='D higher', cex.lab=CL);

DM = DM[DM$strand=='+', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd');

DM = datC[datC$higher=='D', ];
DM = DM[DM$strand=='-', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='rev');
title('distribution of DM sites across strand, social status, and nucleotides',outer=T)


COL='grey';

# figure 2
# 
DAT = datC;
par(mfrow=c(1,3), oma=c(0,1,0,0));
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='overall\n');
DAT = datC[datC$strand=='+',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='fwd strand, overall\n');
DAT = datC[datC$strand=='-',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='rev strand, overall\n');

# figure 3
# 
DAT = datC;
par(mfrow=c(2,2), oma=c(0,1,0,0));
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='D levels\n', names=c('D-win','D-lose'));
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='ND levels\n', names=c('ND-lose','ND-lwin'));
DAT = datC[datC$higher=='D',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='D-win\n');
DAT = datC[datC$higher=='ND',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='ND-win\n');

# figure 4
# same as half of figure 2 but divided by strand
par(mfrow=c(2,3), oma=c(0,1,0,0));
DAT = datC;
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='D levels\n', names=c('D-win','ND-win'));
DAT = datC[datC$strand=='+',];
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='fwd, D levels\n', names=c('D-win','ND-win'));
DAT = datC[datC$strand=='-',];
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='rev, D levels\n', names=c('D-win','ND-win'));
DAT = datC;
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='ND levels\n', names=c('D-win','ND-win'));
DAT = datC[datC$strand=='+',];
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='fwd, ND levels\n', names=c('D-win','ND-win'));
DAT = datC[datC$strand=='-',];
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='rev, ND levels\n', names=c('D-win','ND-win'));

# figure 5
# 
par(mfrow=c(2,3), oma=c(0,1,0,0));
DAT = datC[datC$higher=='D',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='D-win\n');
DAT = datC[datC$higher=='D' & datC$strand=='+',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='fwd, D-win\n');
DAT = datC[datC$higher=='D' & datC$strand=='-',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='rev, D-win\n');
DAT = datC[datC$higher=='ND',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='ND-win\n');
DAT = datC[datC$higher=='ND' & datC$strand=='+',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='fwd, ND-win\n');
DAT = datC[datC$higher=='ND' & datC$strand=='-',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='rev, ND-win\n');

# figure 2-2
# similar to figure 2 but divided by CG/CH
par(mfrow=c(2,3), oma=c(0,1,0,0));
DAT = datC[datC$diH=='CG',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CG, overall\n');
DAT = datC[datC$strand=='+' datC$diH=='CG',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='fwd, CG, overall\n');
DAT = datC[datC$strand=='-' & datC$diH=='CG',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='rev, CG, overall\n');
DAT = datC[datC$diH=='CH',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CH, overall\n');
DAT = datC[datC$strand=='+' datC$diH=='CH',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='fwd, CH, overall\n');
DAT = datC[datC$strand=='-' & datC$diH=='CH',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='rev, CH, overall\n');

# figure 6
# 
par(mfrow=c(2,4), oma=c(0,1,0,0));
DAT = datC[datC$diH=='CG', ];
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CG, D levels\n', names=c('D-win','D-lose'));
DAT = datC[datC$diH=='CH', ];
verboseBoxplot(DAT$D,DAT$higher, col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CH, D levels\n', names=c('D-win','D-lose'));
DAT = datC[datC$diH=='CG', ];
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='CG, ND levels\n', names=c('ND-lose','ND-lwin'));
DAT = datC[datC$diH=='CH', ];
verboseBoxplot(DAT$ND,DAT$higher, col=COL, frame.plot=F, ylab='', xlab='', main='CH, ND levels\n', names=c('ND-lose','ND-lwin'));
DAT = datC[datC$higher=='D' & datC$diH=='CG',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CG, D-win\n');
DAT = datC[datC$higher=='D' & datC$diH=='CH',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='methylation ratio (C/U)', xlab='', main='CH, D-win\n');
DAT = datC[datC$higher=='ND' & datC$diH=='CG',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='CG, ND-win\n');
DAT = datC[datC$higher=='ND' & datC$diH=='CH',];
verboseBoxplot(c(DAT$ND,DAT$D),c(rep('ND',nrow(DAT)),rep('D',nrow(DAT))), col=COL, frame.plot=F, ylab='', xlab='', main='CH, ND-win\n');

# figure X7
# 
par(mfrow=c(2,5), oma=c(0,1,0,0));COL='grey';
verboseBoxplot(datC$highval,datC$higher,frame.plot=F,xlab='',col=COL, ylab='methylation ratio (C/U)',main='all:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$strand=='-'],datC$higher[datC$strand=='-'],frame.plot=F,xlab='',col=COL,ylab='',main='rev:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$strand=='+'],datC$higher[datC$strand=='+'],frame.plot=F,xlab='',col=COL,ylab='',main='fwd:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CG'],datC$higher[datC$diH=='CG'],frame.plot=F,xlab='',col=COL,ylab='',main='CG:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CH'],datC$higher[datC$diH=='CH'],frame.plot=F,xlab='',col=COL,ylab='',main='CH:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CG' & datC$strand=='+'],datC$higher[datC$diH=='CG' & datC$strand=='+'],frame.plot=F,xlab='',col=COL,ylab='',main='CG,fwd:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CG' & datC$strand=='-'],datC$higher[datC$diH=='CG' & datC$strand=='-'],frame.plot=F,xlab='',col=COL,ylab='',main='CG,rev:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CH' & datC$strand=='+'],datC$higher[datC$diH=='CH' & datC$strand=='+'],frame.plot=F,xlab='',col=COL,ylab='',main='CH,fwd:', names=c('D-win','ND-win'));
verboseBoxplot(datC$highval[datC$diH=='CH' & datC$strand=='-'],datC$higher[datC$diH=='CH' & datC$strand=='-'],frame.plot=F,xlab='',col=COL,ylab='',main='CH,rev:', names=c('D-win','ND-win'));

# figure X8
#
par(mfrow=c(2,7), par(oma=c(0,1,0,0)));COL='grey';
verboseBoxplot(datC$highval,datC$strand,frame.plot=F,xlab='',col=COL, ylab='methylation ratio (C/U)',main='all wins:');
verboseBoxplot(datC$highval[datC$higher=='ND'],datC$strand[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND-win:');
verboseBoxplot(datC$highval[datC$higher=='D'],datC$strand[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D-win:');
verboseBoxplot(datC$highval[datC$higher=='ND' & datC$diH=='CG'],datC$strand[datC$higher=='ND' & datC$diH=='CG'],frame.plot=F,xlab='',col=COL, ylab='',main='CG, ND-win:');
verboseBoxplot(datC$highval[datC$higher=='ND' & datC$diH=='CH'],datC$strand[datC$higher=='ND' & datC$diH=='CH'],frame.plot=F,xlab='',col=COL, ylab='',main='CH, ND-win:');
verboseBoxplot(datC$highval[datC$higher=='D' & datC$diH=='CG'],datC$strand[datC$higher=='D' & datC$diH=='CG'],frame.plot=F,xlab='',col=COL, ylab='',main='CG, D-win:');
verboseBoxplot(datC$highval[datC$higher=='D' & datC$diH=='CH'],datC$strand[datC$higher=='D' & datC$diH=='CH'],frame.plot=F,xlab='',col=COL, ylab='',main='CH, D-win:');
verboseBoxplot(datC$highval,datC$diH,frame.plot=F,xlab='',col=COL, ylab='methylation ratio (C/U)',main='all wins:');
verboseBoxplot(datC$highval[datC$higher=='ND'],datC$diH[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND-win:');
verboseBoxplot(datC$highval[datC$higher=='D'],datC$diH[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D-win:');
verboseBoxplot(datC$highval[datC$higher=='ND' & datC$strand=='+'],datC$diH[datC$higher=='ND' & datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd, ND-win:');
verboseBoxplot(datC$highval[datC$higher=='ND' & datC$strand=='-'],datC$diH[datC$higher=='ND' & datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev, ND-win:');
verboseBoxplot(datC$highval[datC$higher=='D' & datC$strand=='+'],datC$diH[datC$higher=='D' & datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd, D-win:');
verboseBoxplot(datC$highval[datC$higher=='D' & datC$strand=='-'],datC$diH[datC$higher=='D' & datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev, D-win:');

# figure 7
par(mfrow=c(2,2), oma=c(0,1,2,0));COL='grey';
DAT = datC[datC$strand=='+', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd, D+ND:');
DAT = datC[datC$strand=='-', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='rev, D+ND:');
DAT = datC[datC$diH=='CG', ];
verboseBoxplot(DAT$highval,DAT$strand,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CG, D+ND:');
DAT = datC[datC$diH=='CH', ];
verboseBoxplot(DAT$highval,DAT$strand,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CH, D+ND:');

# figure 9
#
par(mfrow=c(3,4), oma=c(0,1,2,0));COL='grey';
verboseBoxplot(datC$highval,datC$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='all:');
DAT = datC[datC$di=='CG', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CG:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CG rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CG fwd:');
DAT = datC[datC$diH=='CH', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CH:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CH rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CH fwd:');
DAT = datC[datC$di=='CA', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CA:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CA rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CA fwd:');
DAT = datC[datC$di=='CC', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CC:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CC rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CC fwd:');










par(mfrow=c(4,3), oma=c(0,0,2,0));COL='grey';
DAT = datC[datC$strand=='-', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$diH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$diH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$di[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$di[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$triHH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$triHH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$triH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$triH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
title('rev strand',outer=T)

par(mfrow=c(4,3), oma=c(0,0,2,0));COL='grey';
DAT = datC[datC$strand=='+', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$diH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$diH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$di[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$di[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$triHH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$triHH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
verboseBoxplot(DAT$highval,DAT$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
verboseBoxplot(DAT$highval[DAT$higher=='ND'],DAT$triH[DAT$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
verboseBoxplot(DAT$highval[DAT$higher=='D'],DAT$triH[DAT$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
title('fwd strand',outer=T)

par(mfrow=c(4,3));COL='grey';
verboseBoxplot(datC$highval,datC$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(datC$highval[datC$strand=='-'],datC$diH[datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(datC$highval[datC$strand=='+'],datC$diH[datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(datC$highval,datC$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(datC$highval[datC$strand=='-'],datC$di[datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(datC$highval[datC$strand=='+'],datC$di[datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(datC$highval,datC$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(datC$highval[datC$strand=='-'],datC$triHH[datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(datC$highval[datC$strand=='+'],datC$triHH[datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(datC$highval,datC$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(datC$highval[datC$strand=='-'],datC$triH[datC$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(datC$highval[datC$strand=='+'],datC$triH[datC$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');

par(mfrow=c(4,3), oma=c(0,0,2,0));COL='grey';
DAT = datC[datC$higher=='ND', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$diH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$diH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$di[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$di[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$triHH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$triHH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$triH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$triH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
title('ND',outer=T)

par(mfrow=c(4,3), oma=c(0,0,2,0));COL='grey';
DAT = datC[datC$higher=='D', ];
verboseBoxplot(DAT$highval,DAT$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$diH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$diH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$di[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$di[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$triHH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$triHH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
verboseBoxplot(DAT$highval,DAT$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='fwd+rev:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$triH[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$triH[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='fwd:');
title('D',outer=T)



par(mfrow=c(4,3), oma=c(0,0,2,0));COL='grey';
DAT = datC[datC$triHH=='CGG', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CGG:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CGG rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CGG fwd:');
DAT = datC[datC$triHH=='CGH', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CGH:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CGH rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CGH fwd:');
DAT = datC[datC$triHH=='CHG', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CHG:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CHG rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CHG fwd:');
DAT = datC[datC$triHH=='CHH', ];
verboseBoxplot(DAT$highval,DAT$higher,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='CHH:');
verboseBoxplot(DAT$highval[DAT$strand=='-'],DAT$higher[DAT$strand=='-'],frame.plot=F,xlab='',col=COL, ylab='',main='CHH rev:');
verboseBoxplot(DAT$highval[DAT$strand=='+'],DAT$higher[DAT$strand=='+'],frame.plot=F,xlab='',col=COL, ylab='',main='CHH fwd:');

###########################

# investigate di/trinucleotide distributions

# investigate types of TEs
par(mfrow=c(1,3));
CX = 1.2
DAT = datCTEhit;
pct1 = sort(table(DAT$V10)) / sum(sort(table(DAT$V10))); 
.barplotStackedFromVecs(list(pct1), axes=F, main='fwd+rev', cex.text=CX);   
DAT = datCTE[datCTEhit$strand=='+', ];
pct1 = sort(table(DAT$V10)) / sum(sort(table(DAT$V10))); 
.barplotStackedFromVecs(list(pct1), axes=F, main='fwd', cex.text=CX);    
DAT = datCTE[datCTEhit$strand=='-', ];
pct1 = sort(table(DAT$V10)) / sum(sort(table(DAT$V10))); 
.barplotStackedFromVecs(list(pct1), axes=F, main='rev', cex.text=CX);

# investigate di/trinucleotide distributions

par(mfrow=c(1,2), oma=c(0,0,2,0));
CT = 1.1;
LO = .2;
for (TYPE in c('DNA','LINE')) {
	DM = datCTE[datCTE$V10==TYPE, ];
	pct1 = sort(table(DM$triH) / length(DM$triH));
	pct2 = sort(table(DM$triHH) / length(DM$triHH));
	pct3 = sort(table(DM$di) / length(DM$di));
	.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main=paste(TYPE,' (n=',nrow(DM),')',sep=''));
}

par(mfrow=c(1,5), oma=c(0,0,2,0));
CT = 1.1;
LO = .2;
for (TYPE in names(table(datCTE$V10))) {
	DM = datCTE[datCTE$V10==TYPE, ];
	pct1 = sort(table(DM$triH) / length(DM$triH));
	pct2 = sort(table(DM$triHH) / length(DM$triHH));
	pct3 = sort(table(DM$di) / length(DM$di));
	.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main=paste(TYPE,' (n=',nrow(DM),')',sep=''));
}


DM = datC[datC$strand=='+', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='fwd');

DM = datC[datC$strand=='-', ];
pct1 = sort(table(DM$triH) / length(DM$triH));
pct2 = sort(table(DM$triHH) / length(DM$triHH));
pct3 = sort(table(DM$di) / length(DM$di));
.barplotStackedFromVecs(list(pct3, pct2, pct1), labOffset=LO, cex.text=CT, axes=F, main='rev');


a1 = d$higher=='D';
a2 = d$strand=='+';
a3 = d$diH=='CG';

draw.triple.venn(area1=sum(a1), area2=sum(a2), area3=sum(a3),
				 n12=sum(a1 & a2), n23=sum(a2 & a3), n13=sum(a1 & a3),
				 n123=sum(a1 & a2 & a3),
				 category=c('D','fwd','CG'),
				 fill=c('red','yellow','green')
				 );
				 
				 
a1 = datC$higher=='ND';
a2 = datC$strand=='+';
a3 = datC$diH=='CG';

draw.triple.venn(area1=sum(a1), area2=sum(a2), area3=sum(a3),
				 n12=sum(a1 & a2), n23=sum(a2 & a3), n13=sum(a1 & a3),
				 n123=sum(a1 & a2 & a3),
				 category=c('ND','fwd','CG'),
				 fill=c('red','yellow','green')
				 );				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 
				 



a1 = datC$higher=='D';
a2 = datC$diH=='CG';
a3 = datC$strand=='-';
a4 = datC$strand=='+';

draw.quad.venn(area1=sum(a1), area2=sum(a2), area3=sum(a3), area4=sum(a4),
			   n12=sum(a1 & a2), n13=sum(a1 & a3), n14=sum(a1 & a4),
			   n23=sum(a2 & a3), n24=sum(a2 & a4), n34=sum(a3 & a4),
			   n123=sum(a1 & a2 & a3), n124=sum(a1 & a2 & a4), n134=sum(a1 & a3 & a4), n234=sum(a2 & a3 & a4),
			   n1234=sum(a1 & a2 & a3 & a4),
			   category=c('D','CG','rev','fwd'),
			   fill=c('grey','blue','yellow','green')
			   );
			   
draw.quintuple.venn(area1=sum(datC$higher=='D'),    
					area2=sum(datC$strand=='+'),
					area3=sum(datC$strand=='-'),
					area4=sum(datC$diH=='CG'),
					area5=sum(datC$diH=='CH'),
					n12=sum(datC$higher=='D' & datC$strand=='+'),
					n13=sum(datC$higher=='D' & datC$strand=='-'),
					n14=sum(datC$higher=='D' & datC$diH=='CG'),
					n15=sum(datC$higher=='D' & datC$diH=='CH'),
					n23=sum(datC$strand=='+' & datC$strand=='-'),
					n24=sum(datC$strand=='+' & datC$diH=='CG'),
					n25=sum(datC$strand=='+' & datC$diH=='CH'),
					n34=sum(datC$strand=='-' & datC$diH=='CG'),
					n35=sum(datC$strand=='-' & datC$diH=='CH'),
					n45=sum(datC$diH=='CG' & datC$diH=='CH'),
					n123=sum(datC$higher=='D' & datC$strand=='+' & datC$strand=='-'),
					n124=sum(datC$higher=='D' & datC$strand=='+' & datC$diH=='CG'),
					n125=sum(datC$higher=='D' & datC$strand=='+' & datC$diH=='CH'),
					n134=sum(datC$higher=='D' & datC$strand=='-' & datC$diH=='CG'),
					n135=sum(datC$higher=='D' & datC$strand=='-' & datC$diH=='CH'),
					n145=sum(datC$higher=='D' & datC$diH=='CG' & datC$diH=='CH'),
					n234=sum(datC$strand=='+' & datC$strand=='-' & datC$diH=='CG'),
					n235=sum(datC$strand=='+' & datC$strand=='-' & datC$diH=='CH'),
					n245=sum(datC$strand=='+' & datC$diH=='CG' & datC$diH=='CH'),
					n345=sum(datC$strand=='-' & datC$diH=='CG' & datC$diH=='CH'),
					n1234=sum(datC$higher=='D' & datC$strand=='+' & datC$strand=='-' & datC$diH=='CG'),
					n1235=sum(datC$higher=='D' & datC$strand=='+' & datC$strand=='-' & datC$diH=='CH'),
					n1245=sum(datC$higher=='D' & datC$strand=='+' & datC$diH=='CG' & datC$diH=='CH'),
					n1345=sum(datC$higher=='D' & datC$strand=='-' & datC$diH=='CG' & datC$diH=='CH'),
					n2345=sum(datC$strand=='+' & datC$strand=='-' & datC$diH=='CG' & datC$diH=='CH'),
					n12345=sum(datC$higher=='D' & datC$strand=='+' & datC$strand=='-' & datC$diH=='CG' & datC$diH=='CH'),
					category=c('D','+','-','CG','CH'),
					fill=c('grey','blue','yellow','red','green')
					)			   


# # # par(mfrow=c(4,3));COL='grey';
# # # verboseBoxplot(datC$highval,datC$diH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='ND'],datC$diH[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='D'],datC$diH[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
# # # verboseBoxplot(datC$highval,datC$di,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='ND'],datC$di[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='D'],datC$di[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
# # # verboseBoxplot(datC$highval,datC$triHH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='ND'],datC$triHH[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='D'],datC$triHH[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');
# # # verboseBoxplot(datC$highval,datC$triH,frame.plot=F,xlab='',col=COL, ylab='higher avg.meth.val',main='D+ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='ND'],datC$triH[datC$higher=='ND'],frame.plot=F,xlab='',col=COL, ylab='',main='ND:');
# # # verboseBoxplot(datC$highval[datC$higher=='D'],datC$triH[datC$higher=='D'],frame.plot=F,xlab='',col=COL, ylab='',main='D:');





################################################
######################################
####################################################

par(mfrow=c(2,1));
verboseBoxplot(dnd$highval[dnd$higher=='D'], dnd$closest[dnd$higher=='D'],main='D winner'); 
verboseBoxplot(dnd$highval[dnd$higher=='ND'], dnd$closest[dnd$higher=='ND'], main='ND winner');``