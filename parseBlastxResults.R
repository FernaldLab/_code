file = '~/Documents/_Burtoni_annotations/H_burtoni_rna_blastx_FISH_ENS_top1-part';
x = read.table(file, header=F, comment.char='', fill=T);

.parseBlastResultsBestHits = function (x, blastName='BLASTX') {
	
	for (row in 1:nrow(x)) {
		if (x[row, 1]=='#') {
			lineType = 'comment';
		} else {
			lineType = 'hit';
		}
		if (lineType=='comment') {
			
		} else if (lineType=='hit') {
			
		} else {
			stop();
		}
	}
}
