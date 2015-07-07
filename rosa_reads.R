# read in the file and keep only rows that have counts
dat = as.matrix(read.table('~/Downloads/test2.reads',header=T,row.names=1,sep='\t'));
dat = dat[1:299, ];

# inspect first few rows
head(dat);
# 						 male31 male54 sup22 sup62
# mature-non-collapsed-1     19     46     0     7
# mature-non-collapsed-2    147    951  3417   618
# mature-non-collapsed-3      4     30     0     3
# mature-non-collapsed-4      0      8    19     7
# mature-non-collapsed-5      6     61     8   101
# mature-non-collapsed-6      0      0     0     0

# summary stats for each animal
#  male54 has much higher count levels than the other 3,
#  maybe this is disrupting DESeq?
summary(dat);
#     male31             male54           sup22            sup62       
# Min.   :     0.0   Min.   :     0   Min.   :     0   Min.   :     0  
# 1st Qu.:     4.0   1st Qu.:    30   1st Qu.:     9   1st Qu.:    18  
# Median :    38.0   Median :   387   Median :   123   Median :   237  
# Mean   :  5315.8   Mean   : 11075   Mean   :  7241   Mean   :  7243  
# 3rd Qu.:   475.5   3rd Qu.:  3176   3rd Qu.:  1226   3rd Qu.:  1526  
# Max.   :286573.0   Max.   :601540   Max.   :275120   Max.   :284612

# looked at correlations between animals
#  male54 has lowest correlation values
cor(dat);
#           male31    male54     sup22     sup62
# male31 1.0000000 0.6109084 0.8560500 0.6126752
# male54 0.6109084 1.0000000 0.5730014 0.6459010
# sup22  0.8560500 0.5730014 1.0000000 0.6739083
# sup62  0.6126752 0.6459010 0.6739083 1.0000000

# try quantile normalizing data
#  normalizing function is in the preprocessCore library
library(preprocessCore);
dat2 = normalize.quantiles(dat);
dimnames(dat2) = dimnames(dat);
head(dat2);
#                        male31  male54   sup22  sup62
# mature-non-collapsed-1  95.25  24.750    0.00   5.00
# mature-non-collapsed-2 649.00 524.500 5221.25 490.50
# mature-non-collapsed-3  14.75  15.250    0.00   3.00
# mature-non-collapsed-4   0.00   3.875   36.75   5.00
# mature-non-collapsed-5  19.75  34.000   14.00  90.75
# mature-non-collapsed-6   0.00   0.000    0.00   0.00

summary(dat2);
#     male31             male54             sup22              sup62         
# Min.   :     0.0   Min.   :     0.0   Min.   :     0.0   Min.   :     0.0  
# 1st Qu.:    14.8   1st Qu.:    15.2   1st Qu.:    15.2   1st Qu.:    15.2  
# Median :   196.2   Median :   196.2   Median :   196.2   Median :   200.8  
# Mean   :  7718.6   Mean   :  7718.6   Mean   :  7718.6   Mean   :  7718.6  
# 3rd Qu.:  1600.6   3rd Qu.:  1600.6   3rd Qu.:  1600.6   3rd Qu.:  1600.6  
# Max.   :361961.2   Max.   :361961.2   Max.   :361961.2   Max.   :361961.2 

# normalizing didn't change correalations very much
#  they're slightly higher overall and now male54 and sup62 have similar avg 
cor(dat2);
#           male31    male54     sup22     sup62
# male31 1.0000000 0.6204001 0.9016000 0.5989989
# male54 0.6204001 1.0000000 0.5887634 0.7142533
# sup22  0.9016000 0.5887634 1.0000000 0.6410671
# sup62  0.5989989 0.7142533 0.6410671 1.0000000