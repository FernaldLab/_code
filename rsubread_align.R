rm(list=ls());
setwd("~/Documents/Rsubread");
library(Rsubread);
basename='H_burtoni_v1.assembly';

buildindex(basename=basename,
		 reference='H_burtoni_v1.assembly.fa',
		 colorspace=F,
		 memory=4000,
		 TH_subread=24
		 );	 
        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //=========================== indexBuilder setting ===========================\\
# ||                                                                            ||
# ||                Index name : H_burtoni_v1.assembly                          ||
# ||               Index space : base-space                                     ||
# ||                    Memory : 4000 Mbytes                                    ||
# ||          Repeat threshold : 24 repeats                                     ||
# ||                                                                            ||
# ||               Input files : 1 file in total                                ||
# ||                             o H_burtoni_v1.assembly.fa                     ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //================================= Running ==================================\\
# ||                                                                            ||
# || Check the integrity of provided reference sequences ...                    ||
# || No format issues were found                                                ||
# || Scan uninformative subreads in reference sequences ...                     ||
# ||    8%,   0 mins elapsed, rate=7634.9k bps/s, total=830m                    ||
# ||   16%,   0 mins elapsed, rate=8013.9k bps/s, total=830m                    ||
# ||   24%,   0 mins elapsed, rate=8156.2k bps/s, total=830m                    ||
# ||   33%,   1 mins elapsed, rate=8217.0k bps/s, total=830m                    ||
# ||   41%,   1 mins elapsed, rate=8262.5k bps/s, total=830m                    ||
# ||   49%,   1 mins elapsed, rate=8297.0k bps/s, total=830m                    ||
# ||   58%,   1 mins elapsed, rate=8308.6k bps/s, total=830m                    ||
# ||   66%,   1 mins elapsed, rate=8265.7k bps/s, total=830m                    ||
# ||   74%,   1 mins elapsed, rate=8280.8k bps/s, total=830m                    ||
# ||   83%,   1 mins elapsed, rate=8278.3k bps/s, total=830m                    ||
# ||   91%,   2 mins elapsed, rate=8299.4k bps/s, total=830m                    ||
# ||   99%,   2 mins elapsed, rate=8296.8k bps/s, total=830m                    ||
# || 115476 uninformative subreads were found.                                  ||
# || These subreads were excluded from index building.                          ||
# || Build the index...                                                         ||
# ||    8%,   2 mins elapsed, rate=5188.4k bps/s, total=830m                    ||
# ||   16%,   2 mins elapsed, rate=5294.3k bps/s, total=830m                    ||
# ||   24%,   2 mins elapsed, rate=5305.2k bps/s, total=830m                    ||
# ||   33%,   3 mins elapsed, rate=5347.7k bps/s, total=830m                    ||
# ||   41%,   3 mins elapsed, rate=5386.3k bps/s, total=830m                    ||
# ||   49%,   3 mins elapsed, rate=5408.3k bps/s, total=830m                    ||
# ||   58%,   3 mins elapsed, rate=5424.7k bps/s, total=830m                    ||
# ||   66%,   3 mins elapsed, rate=5436.7k bps/s, total=830m                    ||
# ||   74%,   4 mins elapsed, rate=5472.3k bps/s, total=830m                    ||
# ||   83%,   4 mins elapsed, rate=5491.5k bps/s, total=830m                    ||
# ||   91%,   4 mins elapsed, rate=5530.4k bps/s, total=830m                    ||
# ||   99%,   4 mins elapsed, rate=5438.9k bps/s, total=830m                    ||
# || Save current index block...                                                ||
# ||  [ 0.0% finished ]                                                         ||
# ||  [ 10.0% finished ]                                                        ||
# ||  [ 20.0% finished ]                                                        ||
# ||  [ 30.0% finished ]                                                        ||
# ||  [ 40.0% finished ]                                                        ||
# ||  [ 50.0% finished ]                                                        ||
# ||  [ 60.0% finished ]                                                        ||
# ||  [ 70.0% finished ]                                                        ||
# ||  [ 80.0% finished ]                                                        ||
# ||  [ 90.0% finished ]                                                        ||
# ||  [ 100.0% finished ]                                                       ||
# ||                                                                            ||
# ||                      Total running time: 5.6 minutes.                      ||
# ||            Index H_burtoni_v1.assembly was successfully built!             ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//
rm(list=ls());
setwd("~/Documents/Rsubread");
library(Rsubread);
basename='H_burtoni_v1.assembly';
firstreads=list.files('../_LYNLEY_RNAseq/')[grepl('_1_pf', list.files('../_LYNLEY_RNAseq/'))];
firstreads=paste('../_LYNLEY_RNAseq/', firstreads, sep='');
secondreads=list.files('../_LYNLEY_RNAseq/')[grepl('_2_pf', list.files('../_LYNLEY_RNAseq/'))];
secondreads=paste('../_LYNLEY_RNAseq/', secondreads, sep='');

align(index=basename,
	 readfile1=firstreads, readfile2=secondreads, 
	 input_format='gzFASTQ', output_format='SAM', 
	 output_file=gsub('1_pf.fastq.gz', 'Rsubread.sam', firstreads),
	 nsubreads=10, TH1=3, TH2=1, nthreads=4, indels=5, phredOffset=33,
	 tieBreakQS=F, tieBreakHamming=T, unique=T, nBestLocations=1, 
	 minFragLength=50, maxFragLength=600, PE_orientation='f',
	 nTrim5=0, nTrim3=0, readGroupID=NULL, readGroup=NULL, color2base=F,
	 DP_GapOpenPenalty=-1, DP_GapExtPenalty=0, DP_MismatchPenalty=0, DP_MatchScore=2,
	 reportFusions=F
	 );
	 
	 
	 
	 
	 
########################
###################
###############
###########
########
#####
###
#
        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== subread-align setting ===========================\\
# ||                                                                            ||
# ||           Function : Read alignment                                        ||
# ||            Threads : 4                                                     ||
# ||       Input file 1 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||       Input file 2 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||        Output file : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||         Index name : H_burtoni_v1.assembly                                 ||
# ||       Phred offset : 33                                                    ||
# ||                                                                            ||
# ||    Min read1 votes : 3                                                     ||
# ||    Min read2 votes : 1                                                     ||
# ||  Max fragment size : 600                                                   ||
# ||  Min fragment size : 50                                                    ||
# ||                                                                            ||
# ||         Max indels : 5                                                     ||
# ||  # of Best mapping : 1                                                     ||
# ||     Unique mapping : yes                                                   ||
# ||   Hamming distance : yes                                                   ||
# ||     Quality scores : no                                                    ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //====================== Running (14-Jan-2014 16:46:00) ======================\\
# ||                                                                            ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_ATCACG_1 ... ||
# || The input file contains base space reads.                                  ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_ATCACG_2 ... ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||    0% completed,  20 mins elapsed, total=38873k frags, rate=19.1k/s        ||
# ||    6% completed,  22 mins elapsed, total=38740k frags, rate=20.1k/s        ||
# || Detect indels...                                                           ||
# ||   11% completed,  24 mins elapsed, total=38740k frags, rate=3.0k/s         ||
# ||   11% completed,  24 mins elapsed, total=38740k frags, rate=3.0k/s         ||
# ||   11% completed,  25 mins elapsed, total=38740k frags, rate=3.1k/s         ||
# ||   12% completed,  25 mins elapsed, total=38740k frags, rate=3.1k/s         ||
# ||   12% completed,  25 mins elapsed, total=38740k frags, rate=3.1k/s         ||
# ||   12% completed,  25 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   12% completed,  25 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   12% completed,  25 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   13% completed,  25 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# || Realign fragments...                                                       ||
# ||   13% completed,  26 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# ||   13% completed,  27 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# ||   13% completed,  27 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   14% completed,  28 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   14% completed,  28 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   14% completed,  29 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   14% completed,  29 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   14% completed,  30 mins elapsed, total=38740k frags, rate=3.1k/s         ||
# ||   14% completed,  30 mins elapsed, total=38740k frags, rate=3.1k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   15% completed,  31 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   15% completed,  32 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   16% completed,  32 mins elapsed, total=38740k frags, rate=3.2k/s         ||
# ||   16% completed,  32 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# ||   17% completed,  33 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# ||   17% completed,  33 mins elapsed, total=38740k frags, rate=3.3k/s         ||
# ||   17% completed,  33 mins elapsed, total=38740k frags, rate=3.4k/s         ||
# ||   18% completed,  34 mins elapsed, total=38740k frags, rate=3.4k/s         ||
# ||   18% completed,  34 mins elapsed, total=38740k frags, rate=3.5k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   19% completed,  37 mins elapsed, total=38739k frags, rate=7.2k/s         ||
# ||   25% completed,  39 mins elapsed, total=38739k frags, rate=8.5k/s         ||
# || Detect indels...                                                           ||
# ||   30% completed,  41 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   30% completed,  41 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   30% completed,  41 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   31% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   31% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   31% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   31% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   31% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   32% completed,  42 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# || Realign fragments...                                                       ||
# ||   32% completed,  43 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   32% completed,  44 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   32% completed,  44 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   32% completed,  45 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   33% completed,  46 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   33% completed,  46 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# ||   33% completed,  47 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# ||   33% completed,  47 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   33% completed,  48 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   34% completed,  49 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   34% completed,  49 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   35% completed,  50 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   35% completed,  50 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   35% completed,  51 mins elapsed, total=38739k frags, rate=4.5k/s         ||
# ||   36% completed,  51 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# ||   36% completed,  51 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# ||   37% completed,  52 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# ||   37% completed,  52 mins elapsed, total=38739k frags, rate=4.6k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   37% completed,  55 mins elapsed, total=38739k frags, rate=7.0k/s         ||
# ||   43% completed,  58 mins elapsed, total=38739k frags, rate=7.5k/s         ||
# || Detect indels...                                                           ||
# ||   49% completed,  61 mins elapsed, total=38739k frags, rate=5.2k/s         ||
# ||   49% completed,  61 mins elapsed, total=38739k frags, rate=5.2k/s         ||
# ||   49% completed,  62 mins elapsed, total=38739k frags, rate=5.2k/s         ||
# ||   50% completed,  62 mins elapsed, total=38739k frags, rate=5.2k/s         ||
# ||   50% completed,  62 mins elapsed, total=38739k frags, rate=5.2k/s         ||
# ||   50% completed,  63 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   50% completed,  63 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   50% completed,  63 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   50% completed,  64 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# || Realign fragments...                                                       ||
# ||   51% completed,  65 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   51% completed,  66 mins elapsed, total=38739k frags, rate=5.0k/s         ||
# ||   51% completed,  67 mins elapsed, total=38739k frags, rate=5.0k/s         ||
# ||   51% completed,  67 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   52% completed,  68 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   52% completed,  69 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   52% completed,  70 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   52% completed,  71 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   52% completed,  71 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   53% completed,  73 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   53% completed,  73 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   54% completed,  73 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   54% completed,  74 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   54% completed,  74 mins elapsed, total=38739k frags, rate=4.7k/s         ||
# ||   55% completed,  75 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   55% completed,  75 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   56% completed,  76 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   56% completed,  76 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   56% completed,  79 mins elapsed, total=38738k frags, rate=6.2k/s         ||
# ||   62% completed,  82 mins elapsed, total=38738k frags, rate=6.6k/s         ||
# || Detect indels...                                                           ||
# ||   68% completed,  85 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   68% completed,  85 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   68% completed,  86 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   68% completed,  86 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   69% completed,  86 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   69% completed,  86 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   69% completed,  87 mins elapsed, total=38738k frags, rate=5.2k/s         ||
# ||   69% completed,  87 mins elapsed, total=38738k frags, rate=5.1k/s         ||
# ||   69% completed,  87 mins elapsed, total=38738k frags, rate=5.1k/s         ||
# || Realign fragments...                                                       ||
# ||   70% completed,  89 mins elapsed, total=38738k frags, rate=5.1k/s         ||
# ||   70% completed,  90 mins elapsed, total=38738k frags, rate=5.0k/s         ||
# ||   70% completed,  91 mins elapsed, total=38738k frags, rate=5.0k/s         ||
# ||   70% completed,  91 mins elapsed, total=38738k frags, rate=5.0k/s         ||
# ||   71% completed,  92 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   71% completed,  93 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   71% completed,  94 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   71% completed,  95 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   71% completed,  96 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   72% completed,  97 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   72% completed,  98 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   73% completed,  98 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   73% completed,  98 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   73% completed,  99 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   74% completed,  99 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   74% completed, 100 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   75% completed, 100 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   75% completed, 100 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   75% completed, 104 mins elapsed, total=38738k frags, rate=5.9k/s         ||
# ||   81% completed, 107 mins elapsed, total=38739k frags, rate=6.1k/s         ||
# || Detect indels...                                                           ||
# ||   87% completed, 109 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   87% completed, 110 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   87% completed, 110 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   87% completed, 110 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   88% completed, 111 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   88% completed, 111 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   88% completed, 111 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   88% completed, 112 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   88% completed, 112 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# || Realign fragments...                                                       ||
# ||   89% completed, 113 mins elapsed, total=38739k frags, rate=5.1k/s         ||
# ||   89% completed, 114 mins elapsed, total=38739k frags, rate=5.0k/s         ||
# ||   89% completed, 115 mins elapsed, total=38739k frags, rate=5.0k/s         ||
# ||   89% completed, 116 mins elapsed, total=38739k frags, rate=5.0k/s         ||
# ||   89% completed, 117 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   90% completed, 118 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   90% completed, 119 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   90% completed, 120 mins elapsed, total=38739k frags, rate=4.9k/s         ||
# ||   90% completed, 121 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   91% completed, 122 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   91% completed, 122 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   92% completed, 123 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   92% completed, 123 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   92% completed, 124 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   93% completed, 124 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   93% completed, 124 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   93% completed, 125 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# ||   94% completed, 125 mins elapsed, total=38739k frags, rate=4.8k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   94% completed, 127 mins elapsed, total=38738k frags, rate=5.8k/s         ||
# || Detect indels...                                                           ||
# ||   97% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   97% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 129 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# || Realign fragments...                                                       ||
# ||   98% completed, 130 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 130 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 130 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 130 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 131 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 131 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 131 mins elapsed, total=38738k frags, rate=4.9k/s         ||
# ||   98% completed, 131 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   98% completed, 131 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# || 2038685 fragments were processed. Save the mapping results for them...     ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 132 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 133 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||   99% completed, 133 mins elapsed, total=38738k frags, rate=4.8k/s         ||
# ||                                                                            ||
# ||                          Completed successfully.                           ||
# ||                                                                            ||
# \\============================================================================//

# //================================= Summary ==================================\\
# ||                                                                            ||
# ||          Processed : 38738845 fragments                                    ||
# ||             Mapped : 33988761 fragments (87.7%)                            ||
# ||   Correctly paired : 31156285 fragments                                    ||
# ||             Indels : 217298                                                ||
# ||                                                                            ||
# ||       Running time : 133.3 minutes                                         ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# #         ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== subread-align setting ===========================\\
# ||                                                                            ||
# ||           Function : Read alignment                                        ||
# ||            Threads : 4                                                     ||
# ||       Input file 1 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||       Input file 2 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||        Output file : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||         Index name : H_burtoni_v1.assembly                                 ||
# ||       Phred offset : 33                                                    ||
# ||                                                                            ||
# ||    Min read1 votes : 3                                                     ||
# ||    Min read2 votes : 1                                                     ||
# ||  Max fragment size : 600                                                   ||
# ||  Min fragment size : 50                                                    ||
# ||                                                                            ||
# ||         Max indels : 5                                                     ||
# ||  # of Best mapping : 1                                                     ||
# ||     Unique mapping : yes                                                   ||
# ||   Hamming distance : yes                                                   ||
# ||     Quality scores : no                                                    ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //====================== Running (15-Jan-2014 14:33:13) ======================\\
# ||                                                                            ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_CGATGT_1 ... ||
# || The input file contains base space reads.                                  ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_CGATGT_2 ... ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||    0% completed,  24 mins elapsed, total=41874k frags, rate=19.7k/s        ||
# ||    6% completed,  26 mins elapsed, total=41728k frags, rate=19.2k/s        ||
# || Detect indels...                                                           ||
# ||   10% completed,  28 mins elapsed, total=41728k frags, rate=2.6k/s         ||
# ||   10% completed,  28 mins elapsed, total=41728k frags, rate=2.6k/s         ||
# ||   11% completed,  28 mins elapsed, total=41728k frags, rate=2.7k/s         ||
# ||   11% completed,  28 mins elapsed, total=41728k frags, rate=2.7k/s         ||
# ||   11% completed,  29 mins elapsed, total=41728k frags, rate=2.7k/s         ||
# ||   11% completed,  29 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   11% completed,  29 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   11% completed,  29 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   12% completed,  29 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# || Realign fragments...                                                       ||
# ||   12% completed,  30 mins elapsed, total=41728k frags, rate=2.9k/s         ||
# ||   12% completed,  30 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   12% completed,  31 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  32 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  32 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  33 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  33 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  34 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   13% completed,  35 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   14% completed,  35 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   14% completed,  36 mins elapsed, total=41728k frags, rate=2.8k/s         ||
# ||   15% completed,  36 mins elapsed, total=41728k frags, rate=2.9k/s         ||
# ||   15% completed,  37 mins elapsed, total=41728k frags, rate=2.9k/s         ||
# ||   15% completed,  37 mins elapsed, total=41728k frags, rate=2.9k/s         ||
# ||   16% completed,  37 mins elapsed, total=41728k frags, rate=3.0k/s         ||
# ||   16% completed,  38 mins elapsed, total=41728k frags, rate=3.0k/s         ||
# ||   16% completed,  38 mins elapsed, total=41728k frags, rate=3.0k/s         ||
# ||   17% completed,  38 mins elapsed, total=41728k frags, rate=3.1k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   17% completed,  42 mins elapsed, total=41727k frags, rate=7.0k/s         ||
# ||   23% completed,  44 mins elapsed, total=41727k frags, rate=8.4k/s         ||
# || Detect indels...                                                           ||
# ||   28% completed,  45 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   28% completed,  46 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   28% completed,  46 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   28% completed,  46 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   29% completed,  46 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   29% completed,  46 mins elapsed, total=41727k frags, rate=4.4k/s         ||
# ||   29% completed,  46 mins elapsed, total=41727k frags, rate=4.4k/s         ||
# ||   29% completed,  46 mins elapsed, total=41727k frags, rate=4.4k/s         ||
# ||   29% completed,  47 mins elapsed, total=41727k frags, rate=4.4k/s         ||
# || Realign fragments...                                                       ||
# ||   30% completed,  47 mins elapsed, total=41727k frags, rate=4.4k/s         ||
# ||   30% completed,  48 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   30% completed,  49 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   30% completed,  49 mins elapsed, total=41727k frags, rate=4.3k/s         ||
# ||   30% completed,  50 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   30% completed,  51 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   31% completed,  51 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   31% completed,  52 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   31% completed,  53 mins elapsed, total=41727k frags, rate=4.1k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   32% completed,  54 mins elapsed, total=41727k frags, rate=4.1k/s         ||
# ||   32% completed,  54 mins elapsed, total=41727k frags, rate=4.1k/s         ||
# ||   32% completed,  54 mins elapsed, total=41727k frags, rate=4.1k/s         ||
# ||   33% completed,  55 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   33% completed,  55 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   33% completed,  56 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   34% completed,  56 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   34% completed,  56 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# ||   34% completed,  57 mins elapsed, total=41727k frags, rate=4.2k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   35% completed,  60 mins elapsed, total=41727k frags, rate=6.8k/s         ||
# ||   41% completed,  63 mins elapsed, total=41726k frags, rate=7.4k/s         ||
# || Detect indels...                                                           ||
# ||   45% completed,  66 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  66 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  66 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  66 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  66 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  67 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   46% completed,  67 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   47% completed,  67 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   47% completed,  67 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# || Realign fragments...                                                       ||
# ||   47% completed,  68 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   47% completed,  69 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   48% completed,  69 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   48% completed,  70 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   48% completed,  71 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   48% completed,  72 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   48% completed,  72 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   48% completed,  73 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   49% completed,  74 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   49% completed,  75 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   49% completed,  75 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   50% completed,  76 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   50% completed,  76 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   51% completed,  77 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   51% completed,  77 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   51% completed,  77 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   52% completed,  78 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   52% completed,  78 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   52% completed,  82 mins elapsed, total=41726k frags, rate=6.4k/s         ||
# ||   58% completed,  85 mins elapsed, total=41726k frags, rate=6.7k/s         ||
# || Detect indels...                                                           ||
# ||   63% completed,  88 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   63% completed,  88 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   63% completed,  88 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  89 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  89 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  89 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  90 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  90 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# ||   64% completed,  90 mins elapsed, total=41726k frags, rate=5.0k/s         ||
# || Realign fragments...                                                       ||
# ||   65% completed,  92 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   65% completed,  93 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   65% completed,  94 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   65% completed,  95 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   65% completed,  96 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   66% completed,  97 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   66% completed,  97 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   66% completed,  98 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   66% completed,  99 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   67% completed, 101 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   67% completed, 101 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   67% completed, 102 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   68% completed, 102 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   68% completed, 102 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   68% completed, 103 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   69% completed, 103 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   69% completed, 104 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   70% completed, 104 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   70% completed, 107 mins elapsed, total=41726k frags, rate=5.9k/s         ||
# ||   76% completed, 111 mins elapsed, total=41726k frags, rate=6.1k/s         ||
# || Detect indels...                                                           ||
# ||   81% completed, 114 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   81% completed, 114 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   81% completed, 114 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   81% completed, 115 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   81% completed, 115 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   81% completed, 115 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   82% completed, 116 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   82% completed, 116 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   82% completed, 117 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# || Realign fragments...                                                       ||
# ||   82% completed, 118 mins elapsed, total=41726k frags, rate=4.9k/s         ||
# ||   83% completed, 119 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   83% completed, 120 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   83% completed, 121 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   83% completed, 122 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   83% completed, 123 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   83% completed, 124 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   84% completed, 125 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   84% completed, 126 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   84% completed, 127 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   85% completed, 128 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   85% completed, 128 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   85% completed, 128 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   86% completed, 129 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   86% completed, 129 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   86% completed, 130 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   87% completed, 130 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   87% completed, 130 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   88% completed, 133 mins elapsed, total=41726k frags, rate=5.6k/s         ||
# ||   94% completed, 136 mins elapsed, total=41726k frags, rate=5.8k/s         ||
# || Detect indels...                                                           ||
# ||   95% completed, 137 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   95% completed, 137 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   95% completed, 138 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   95% completed, 138 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   95% completed, 138 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   95% completed, 139 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   96% completed, 139 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   96% completed, 139 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   96% completed, 139 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# || Realign fragments...                                                       ||
# ||   96% completed, 140 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   96% completed, 141 mins elapsed, total=41726k frags, rate=4.8k/s         ||
# ||   96% completed, 141 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   96% completed, 142 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   96% completed, 143 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   97% completed, 143 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   97% completed, 144 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   97% completed, 145 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||   97% completed, 145 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# || 5026116 fragments were processed. Save the mapping results for them...     ||
# ||   97% completed, 146 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   98% completed, 147 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   98% completed, 147 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   98% completed, 147 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   98% completed, 148 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   99% completed, 148 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   99% completed, 148 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   99% completed, 148 mins elapsed, total=41726k frags, rate=4.6k/s         ||
# ||   99% completed, 149 mins elapsed, total=41726k frags, rate=4.7k/s         ||
# ||                                                                            ||
# ||                          Completed successfully.                           ||
# ||                                                                            ||
# \\============================================================================//

# //================================= Summary ==================================\\
# ||                                                                            ||
# ||          Processed : 41726276 fragments                                    ||
# ||             Mapped : 35297834 fragments (84.6%)                            ||
# ||   Correctly paired : 32686187 fragments                                    ||
# ||             Indels : 229765                                                ||
# ||                                                                            ||
# ||       Running time : 149.4 minutes                                         ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//


        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== subread-align setting ===========================\\
# ||                                                                            ||
# ||           Function : Read alignment                                        ||
# ||            Threads : 4                                                     ||
# ||       Input file 1 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||       Input file 2 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||        Output file : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||         Index name : H_burtoni_v1.assembly                                 ||
# ||       Phred offset : 33                                                    ||
# ||                                                                            ||
# ||    Min read1 votes : 3                                                     ||
# ||    Min read2 votes : 1                                                     ||
# ||  Max fragment size : 600                                                   ||
# ||  Min fragment size : 50                                                    ||
# ||                                                                            ||
# ||         Max indels : 5                                                     ||
# ||  # of Best mapping : 1                                                     ||
# ||     Unique mapping : yes                                                   ||
# ||   Hamming distance : yes                                                   ||
# ||     Quality scores : no                                                    ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //====================== Running (15-Jan-2014 17:02:39) ======================\\
# ||                                                                            ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TGACCA_1 ... ||
# || The input file contains base space reads.                                  ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TGACCA_2 ... ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||    0% completed,  21 mins elapsed, total=38604k frags, rate=11.1k/s        ||
# ||    6% completed,  24 mins elapsed, total=38470k frags, rate=11.2k/s        ||
# || Detect indels...                                                           ||
# ||   11% completed,  28 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   11% completed,  28 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  29 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  29 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  29 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  30 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  30 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   12% completed,  31 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# ||   13% completed,  31 mins elapsed, total=38470k frags, rate=2.7k/s         ||
# || Realign fragments...                                                       ||
# ||   13% completed,  33 mins elapsed, total=38470k frags, rate=2.6k/s         ||
# ||   13% completed,  34 mins elapsed, total=38470k frags, rate=2.6k/s         ||
# ||   13% completed,  35 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   14% completed,  36 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   14% completed,  37 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   14% completed,  38 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   14% completed,  39 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   14% completed,  40 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   15% completed,  41 mins elapsed, total=38470k frags, rate=2.3k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   15% completed,  42 mins elapsed, total=38470k frags, rate=2.3k/s         ||
# ||   16% completed,  43 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   16% completed,  43 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   16% completed,  44 mins elapsed, total=38470k frags, rate=2.4k/s         ||
# ||   17% completed,  44 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   17% completed,  44 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   17% completed,  45 mins elapsed, total=38470k frags, rate=2.5k/s         ||
# ||   18% completed,  45 mins elapsed, total=38470k frags, rate=2.6k/s         ||
# ||   18% completed,  46 mins elapsed, total=38470k frags, rate=2.6k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   19% completed,  49 mins elapsed, total=38469k frags, rate=4.4k/s         ||
# ||   25% completed,  52 mins elapsed, total=38470k frags, rate=5.1k/s         ||
# || Detect indels...                                                           ||
# ||   30% completed,  55 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   30% completed,  56 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   31% completed,  56 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   31% completed,  57 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   31% completed,  57 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   31% completed,  58 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   31% completed,  58 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   32% completed,  59 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# ||   32% completed,  59 mins elapsed, total=38470k frags, rate=3.5k/s         ||
# || Realign fragments...                                                       ||
# ||   32% completed,  61 mins elapsed, total=38470k frags, rate=3.4k/s         ||
# ||   32% completed,  62 mins elapsed, total=38470k frags, rate=3.4k/s         ||
# ||   33% completed,  63 mins elapsed, total=38470k frags, rate=3.3k/s         ||
# ||   33% completed,  64 mins elapsed, total=38470k frags, rate=3.3k/s         ||
# ||   33% completed,  65 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   33% completed,  67 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   33% completed,  68 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   33% completed,  69 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# ||   34% completed,  70 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   34% completed,  72 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# ||   35% completed,  72 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# ||   35% completed,  72 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# ||   35% completed,  73 mins elapsed, total=38470k frags, rate=3.1k/s         ||
# ||   36% completed,  73 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   36% completed,  74 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   37% completed,  74 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   37% completed,  74 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# ||   37% completed,  75 mins elapsed, total=38470k frags, rate=3.2k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   38% completed,  78 mins elapsed, total=38469k frags, rate=4.3k/s         ||
# ||   44% completed,  81 mins elapsed, total=38469k frags, rate=4.7k/s         ||
# || Detect indels...                                                           ||
# ||   49% completed,  85 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   49% completed,  85 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   50% completed,  85 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   50% completed,  86 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   50% completed,  86 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   50% completed,  87 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   50% completed,  87 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   51% completed,  88 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   51% completed,  88 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# || Realign fragments...                                                       ||
# ||   51% completed,  90 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   51% completed,  91 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   52% completed,  92 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   52% completed,  94 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   52% completed,  95 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# ||   52% completed,  96 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# ||   52% completed,  97 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# ||   53% completed,  99 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   53% completed, 100 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   53% completed, 101 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   54% completed, 102 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   54% completed, 102 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   54% completed, 103 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   55% completed, 103 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   55% completed, 103 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   56% completed, 104 mins elapsed, total=38469k frags, rate=3.4k/s         ||
# ||   56% completed, 104 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# ||   56% completed, 105 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   57% completed, 108 mins elapsed, total=38469k frags, rate=4.2k/s         ||
# ||   63% completed, 111 mins elapsed, total=38469k frags, rate=4.5k/s         ||
# || Detect indels...                                                           ||
# ||   68% completed, 114 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   69% completed, 115 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   69% completed, 115 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   69% completed, 115 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   69% completed, 116 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   69% completed, 116 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   70% completed, 117 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   70% completed, 117 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   70% completed, 118 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# || Realign fragments...                                                       ||
# ||   70% completed, 119 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   70% completed, 121 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   71% completed, 122 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   71% completed, 123 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   71% completed, 125 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   71% completed, 126 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   71% completed, 127 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   72% completed, 128 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   72% completed, 130 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   72% completed, 131 mins elapsed, total=38469k frags, rate=3.5k/s         ||
# ||   73% completed, 132 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   73% completed, 132 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   74% completed, 133 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   74% completed, 133 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   74% completed, 133 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   75% completed, 134 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   75% completed, 134 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   75% completed, 134 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   76% completed, 138 mins elapsed, total=38469k frags, rate=4.2k/s         ||
# ||   82% completed, 141 mins elapsed, total=38469k frags, rate=4.4k/s         ||
# || Detect indels...                                                           ||
# ||   87% completed, 144 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   88% completed, 145 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   88% completed, 145 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   88% completed, 146 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   88% completed, 146 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   88% completed, 146 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   89% completed, 147 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   89% completed, 147 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# ||   89% completed, 148 mins elapsed, total=38469k frags, rate=3.9k/s         ||
# || Realign fragments...                                                       ||
# ||   89% completed, 149 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   90% completed, 151 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   90% completed, 152 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   90% completed, 153 mins elapsed, total=38469k frags, rate=3.8k/s         ||
# ||   90% completed, 155 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   90% completed, 156 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   91% completed, 157 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   91% completed, 159 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   91% completed, 160 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   91% completed, 162 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   92% completed, 162 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   92% completed, 163 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   93% completed, 163 mins elapsed, total=38469k frags, rate=3.6k/s         ||
# ||   93% completed, 164 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   93% completed, 164 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   94% completed, 164 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   94% completed, 165 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   95% completed, 165 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   95% completed, 167 mins elapsed, total=38469k frags, rate=4.2k/s         ||
# || Detect indels...                                                           ||
# ||   98% completed, 168 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 168 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 168 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 168 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# || Realign fragments...                                                       ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 169 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 170 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 170 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 170 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 171 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 171 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   98% completed, 171 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 171 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# || 1768859 fragments were processed. Save the mapping results for them...     ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 172 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||   99% completed, 173 mins elapsed, total=38469k frags, rate=3.7k/s         ||
# ||                                                                            ||
# ||                          Completed successfully.                           ||
# ||                                                                            ||
# \\============================================================================//

# //================================= Summary ==================================\\
# ||                                                                            ||
# ||          Processed : 38469019 fragments                                    ||
# ||             Mapped : 33082622 fragments (86.0%)                            ||
# ||   Correctly paired : 30722517 fragments                                    ||
# ||             Indels : 165387                                                ||
# ||                                                                            ||
# ||       Running time : 173.1 minutes                                         ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//


        # ==========     _____ _    _ ____  _____  ______          _____  
        # =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          # =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            # ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              # ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        # ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
       # Rsubread 1.12.6

# //========================== subread-align setting ===========================\\
# ||                                                                            ||
# ||           Function : Read alignment                                        ||
# ||            Threads : 4                                                     ||
# ||       Input file 1 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||       Input file 2 : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||        Output file : ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L ... ||
# ||         Index name : H_burtoni_v1.assembly                                 ||
# ||       Phred offset : 33                                                    ||
# ||                                                                            ||
# ||    Min read1 votes : 3                                                     ||
# ||    Min read2 votes : 1                                                     ||
# ||  Max fragment size : 600                                                   ||
# ||  Min fragment size : 50                                                    ||
# ||                                                                            ||
# ||         Max indels : 5                                                     ||
# ||  # of Best mapping : 1                                                     ||
# ||     Unique mapping : yes                                                   ||
# ||   Hamming distance : yes                                                   ||
# ||     Quality scores : no                                                    ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//

# //====================== Running (15-Jan-2014 19:55:48) ======================\\
# ||                                                                            ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TTAGGC_1 ... ||
# || The input file contains base space reads.                                  ||
# || Decompress ../_LYNLEY_RNAseq/130913_LYNLEY_0364_AC2HRPACXX_L6_TTAGGC_2 ... ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||    0% completed,  20 mins elapsed, total=37727k frags, rate=11.0k/s        ||
# ||    6% completed,  23 mins elapsed, total=37596k frags, rate=12.4k/s        ||
# || Detect indels...                                                           ||
# ||   11% completed,  27 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   12% completed,  27 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   12% completed,  27 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   12% completed,  28 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   12% completed,  28 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   12% completed,  29 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   13% completed,  29 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   13% completed,  30 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   13% completed,  30 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# || Realign fragments...                                                       ||
# ||   13% completed,  31 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   14% completed,  32 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   14% completed,  33 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   14% completed,  33 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   14% completed,  34 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   14% completed,  35 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   15% completed,  36 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   15% completed,  36 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   15% completed,  37 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   16% completed,  38 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   16% completed,  39 mins elapsed, total=37596k frags, rate=2.6k/s         ||
# ||   16% completed,  39 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   17% completed,  40 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   17% completed,  40 mins elapsed, total=37596k frags, rate=2.7k/s         ||
# ||   17% completed,  40 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   18% completed,  41 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   18% completed,  41 mins elapsed, total=37596k frags, rate=2.8k/s         ||
# ||   19% completed,  42 mins elapsed, total=37596k frags, rate=2.9k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   19% completed,  45 mins elapsed, total=37595k frags, rate=5.0k/s         ||
# ||   25% completed,  48 mins elapsed, total=37595k frags, rate=5.9k/s         ||
# || Detect indels...                                                           ||
# ||   31% completed,  51 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   31% completed,  51 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   31% completed,  52 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  52 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  53 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  53 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  53 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  54 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   32% completed,  54 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# || Realign fragments...                                                       ||
# ||   33% completed,  56 mins elapsed, total=37595k frags, rate=3.7k/s         ||
# ||   33% completed,  56 mins elapsed, total=37595k frags, rate=3.7k/s         ||
# ||   33% completed,  57 mins elapsed, total=37595k frags, rate=3.7k/s         ||
# ||   33% completed,  58 mins elapsed, total=37595k frags, rate=3.6k/s         ||
# ||   34% completed,  59 mins elapsed, total=37595k frags, rate=3.6k/s         ||
# ||   34% completed,  60 mins elapsed, total=37595k frags, rate=3.6k/s         ||
# ||   34% completed,  61 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   34% completed,  62 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   34% completed,  63 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   35% completed,  64 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   35% completed,  64 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   36% completed,  65 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   36% completed,  65 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   37% completed,  66 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   37% completed,  66 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   37% completed,  66 mins elapsed, total=37595k frags, rate=3.5k/s         ||
# ||   38% completed,  67 mins elapsed, total=37595k frags, rate=3.6k/s         ||
# ||   38% completed,  67 mins elapsed, total=37595k frags, rate=3.6k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   39% completed,  71 mins elapsed, total=37595k frags, rate=4.9k/s         ||
# ||   45% completed,  74 mins elapsed, total=37595k frags, rate=5.3k/s         ||
# || Detect indels...                                                           ||
# ||   50% completed,  77 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   51% completed,  78 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   51% completed,  78 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   51% completed,  78 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   51% completed,  79 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   51% completed,  79 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   52% completed,  80 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   52% completed,  80 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   52% completed,  81 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# || Realign fragments...                                                       ||
# ||   52% completed,  82 mins elapsed, total=37595k frags, rate=4.0k/s         ||
# ||   53% completed,  83 mins elapsed, total=37595k frags, rate=4.0k/s         ||
# ||   53% completed,  84 mins elapsed, total=37595k frags, rate=4.0k/s         ||
# ||   53% completed,  85 mins elapsed, total=37595k frags, rate=3.9k/s         ||
# ||   53% completed,  86 mins elapsed, total=37595k frags, rate=3.9k/s         ||
# ||   53% completed,  87 mins elapsed, total=37595k frags, rate=3.9k/s         ||
# ||   54% completed,  88 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   54% completed,  89 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   54% completed,  89 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   55% completed,  91 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   55% completed,  91 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   55% completed,  92 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   56% completed,  92 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   56% completed,  93 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   57% completed,  93 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   57% completed,  93 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   57% completed,  94 mins elapsed, total=37595k frags, rate=3.8k/s         ||
# ||   58% completed,  94 mins elapsed, total=37595k frags, rate=3.9k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   58% completed,  97 mins elapsed, total=37594k frags, rate=4.8k/s         ||
# ||   64% completed, 100 mins elapsed, total=37594k frags, rate=5.1k/s         ||
# || Detect indels...                                                           ||
# ||   70% completed, 104 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   70% completed, 104 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   70% completed, 104 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   71% completed, 105 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   71% completed, 105 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   71% completed, 106 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   71% completed, 106 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   71% completed, 106 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   72% completed, 107 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# || Realign fragments...                                                       ||
# ||   72% completed, 108 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   72% completed, 109 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   72% completed, 110 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   73% completed, 111 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   73% completed, 112 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   73% completed, 113 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   73% completed, 114 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   73% completed, 115 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   73% completed, 116 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   74% completed, 117 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   74% completed, 118 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   75% completed, 118 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   75% completed, 119 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   76% completed, 119 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   76% completed, 119 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   76% completed, 120 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   77% completed, 120 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# ||   77% completed, 120 mins elapsed, total=37594k frags, rate=4.0k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   78% completed, 124 mins elapsed, total=37595k frags, rate=4.7k/s         ||
# ||   84% completed, 127 mins elapsed, total=37595k frags, rate=4.9k/s         ||
# || Detect indels...                                                           ||
# ||   90% completed, 130 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   90% completed, 130 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   90% completed, 131 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   90% completed, 131 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   90% completed, 132 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   90% completed, 132 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   91% completed, 132 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   91% completed, 133 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   91% completed, 133 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# || Realign fragments...                                                       ||
# ||   91% completed, 135 mins elapsed, total=37595k frags, rate=4.3k/s         ||
# ||   92% completed, 136 mins elapsed, total=37595k frags, rate=4.2k/s         ||
# ||   92% completed, 137 mins elapsed, total=37595k frags, rate=4.2k/s         ||
# ||   92% completed, 138 mins elapsed, total=37595k frags, rate=4.2k/s         ||
# ||   92% completed, 139 mins elapsed, total=37595k frags, rate=4.2k/s         ||
# ||   92% completed, 140 mins elapsed, total=37595k frags, rate=4.2k/s         ||
# ||   93% completed, 141 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   93% completed, 141 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   93% completed, 142 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# || 7340032 fragments were processed. Save the mapping results for them...     ||
# ||   94% completed, 144 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   94% completed, 144 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   94% completed, 145 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   95% completed, 145 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   95% completed, 146 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   96% completed, 146 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   96% completed, 146 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   96% completed, 147 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# ||   97% completed, 147 mins elapsed, total=37595k frags, rate=4.1k/s         ||
# || Scan read files for multi-threaded alignment...                            ||
# || Load the 1-th index block...                                               ||
# || Map fragments...                                                           ||
# ||   97% completed, 148 mins elapsed, total=37594k frags, rate=4.8k/s         ||
# || Detect indels...                                                           ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.2k/s         ||
# || Realign fragments...                                                       ||
# ||   99% completed, 149 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# || 894588 fragments were processed. Save the mapping results for them...      ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 150 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||   99% completed, 151 mins elapsed, total=37594k frags, rate=4.1k/s         ||
# ||                                                                            ||
# ||                          Completed successfully.                           ||
# ||                                                                            ||
# \\============================================================================//

# //================================= Summary ==================================\\
# ||                                                                            ||
# ||          Processed : 37594748 fragments                                    ||
# ||             Mapped : 32395235 fragments (86.2%)                            ||
# ||   Correctly paired : 29870391 fragments                                    ||
# ||             Indels : 197225                                                ||
# ||                                                                            ||
# ||       Running time : 151.4 minutes                                         ||
# ||                                                                            ||
# \\===================== http://subread.sourceforge.net/ ======================//