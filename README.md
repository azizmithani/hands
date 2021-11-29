# hands
HANDS: a tool for genome-wide discovery of subgenome-specific base-identity in polyploids


HANDS is depreciated. Please use HANDS2.
HANDS v1.1
Usage: java -jar hands.jar <input parameters>
Input Parameters
	-h or -help  :   Display this help
	-s <str>     :   Polyploid SAM file
	-g <str>     :   GFF file containing gene coordinates
	-hsp <str>   :   Polyploid HSP file
	-snp<n> <str>:   Diploid # n SNP file
	-bd <str>    :   Polyploid Base Distribution File (Optional)
	-bd<1> <str> :   Diploid # n Base Distribution File (Optional)
	-out<n> <str>:   Sub-Genome # n Output File
	-sp <double> :   SNP Pair Proportion Threshold (Default: 0.05)
	-pm <double> :   SNP Pattern Matching Threshold (Default: 0.5)
	-r <boolean> :   Rectify Assignment using Reference Genome (Default: FALSE)
	-d <int>     :   Use Genome <int> as Distant Genome (Default: <null>)
Note: At most one Diploid SNP file can be missing. Use "" for the missing file.
      HANDS supports up to 10 genomes.
