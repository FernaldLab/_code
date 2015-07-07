rm(list=ls()); options(stringsAsFactors=F);
setwd('/Volumes/fishstudies-1/_behaviorScoring/_testing');

d0 = list();
for (f in list.files()) {
	tmp = read.csv(f);
	rownames(tmp)= tmp[,1];
	d0[[f]] = tmp[,-1];
}; rm(f,tmp);

f30bs = d0$'first30_Baseline Days_data.csv';
f30ex = d0$'first30_Experimental Days_data.csv';
f30l30bs = d0$'first30andlast30_Baseline Days_data.csv';
f30l30ex = d0$'first30andlast30_Experimental Days_data.csv';

diff_f = f30ex - f30bs;
diff_fl = f30l30ex - f30l30bs;

diff_f_c = diff_f[grep('Count',rownames(diff_f)), ];
diff_fl_c = diff_fl[grep('Count',rownames(diff_fl)), ];

cRows = grep('Count',rownames(f30ex));

.pcChange = function (num1, num2) {
	return( (num2-num1)/num1*100 );
}