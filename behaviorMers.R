# start with vector of behaviors 'b'
# e.g. 
# load('/Volumes/fishstudies/beh.RData')


# want to examine single behavior 'f'
# get positions
bpos = which(b=='f');
numpreceding = 2;
# build matrix of preceding behaviors
m = as.data.frame(matrix(rep(NA, length(bpos)*numpreceding), ncol=numpreceding+1))
for (i in 1:length(bpos)) {
	m[i, ] = b[(bpos[i]-numpreceding):bpos[i]];
}; rm(i)

# get probabilities of preceding sequences
cc = table(apply(m[,1:2], 1, function(f) paste(f,collapse='')));
ccprob = cc / nrow(m);
#         cc         cg         fc         gc         hd 
# 0.08333333 0.41666667 0.25000000 0.16666667 0.08333333

# how often does 'cg' occur overall?
cgc = .countBehaviorSequences(b,'cg')$count;
cgprob = cgc / length(b);

# does it occur more than chance?
cgover = .computeOverrepPval(b,'cg')
# what about 'cgf'?
cgfover = .computeOverrepPval(b,'cgf')

# when 'cg' is observed, how often does it precede 'f'?
cc[names(cc)=='cg'] / cgc;

# when 'f' is observed, how often was it preceded by 'cg'?
ccprob[names(ccprob)=='cg']



# see what follows 'cg'
cgpos = .countBehaviorSequences(b,'cg')$pos
tmp = b[sort(c(cgpos, cgpos+1, cgpos+2))];
cgfollowcounts = table(tmp[seq(3,length(tmp),3)]);