# REDitools

---


## Guidelines to run

1. Based on a paper: 

REDItools implements several positional filters to minimize biases due to sequencing and mapping errors. Below, we report a description of filters and default values, but users should tune them according to the characteristics of input RNAseq or WGS reads. In brackets, we also report command line options to active each filter in REDItoolsDnaRna.py script.

| Quality score | positions with a Phred score <25 are excluded (-q); |
|---|---|
| Mapping quality | reads with a mapping quality score <30 are removed (-m in combination with -u for RNAseq and -U for WGS/WXS). This score is aligner-dependent and should be changed accordingly; |
| Per base coverage | sites not supported by ≥10 reads are filtered out (-c); |
| Bases supporting variation | sites with at least three reads supporting the variation are excluded (-v); |
| Homopolymeric regions | positions in homopolymeric stretches of at least five bases are removed (-O in combination with -l for RNAseq and -L for WGS/WXS); |
| Splice sites | positions near (i.e., within 4 bases) known splice sites are excluded (-V in combination with -w). This filter requires an input table with the list of known spice sites; |
| Mis-mapping correction | positions in reads mapping to multiple genome locations are excluded (-b for RNAseq or -B for WGS/WXS). This filter requires an input list of mis-mapping reads from Blat or similar tools; |
| Unique reads | only reads uniquely mapping are taken into account (-e for RNAseq and -E for WGS/WXS); |
| Paired-end reads | only concordant reads are used (-p for RNAseq and -P for WGS/WXS). This filter works with paired-end reads only; |
| PCR duplicates | reads marked as PCR duplicates are removed (-d for RNAseq and -D for WGS/WXS). This filter requires the preprocessing of the input BAM file with tools able to mark duplicate reads; |
| Trimming | positions located a few bases upstream and/or downstream of single reads are excluded (-a for RNAseq and -A WGS/WXS); |
| Infer strand | positions are handled according to the read strand (-s in combination with -g). This filter requires reads sequenced from a strand-oriented library. For reads that are not strand-oriented, the strand can be inferred using gene annotations provided in gtf format (-G); |
| Strand correction | positions in opposite orientation to the detected strand are removed (-S in combination with -s); |
| Editing level | sites showing editing levels lower than the background value of 0.1 are removed (-n); |
| Multiple changes | positions showing multiple changes are excluded (-z for RNAseq and -Z for WGS/WXS); |
| WGS/WXS support | positions not supported by WGS/WXS reads are discarded (-V); |
| Invariant sites | positions without a reference mismatch are removed (-R); |
| Selected changes | positions not showing a specific substitution type are removed (-W); |
| Specific sites | positions not included in the user-provided list of sites are excluded (-T) or included (-K); |
| Specific genomic region | positions not included in the user-provided genomic region are removed (-Y followed by the genomic region in the format chr |start-end). |


## Scripts:


Reditools2 : https://github.com/hato-lab/A-to-I-edit/blob/2a09982e6b61014774aacadff838af55f14ff5a4/Alignment_REDItools2_scripts.sh#L245
