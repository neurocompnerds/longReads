#!/bin/bash

#SBATCH -J ngmlr
#SBATCH -o /fast/users/%u/launch/ngmlr.slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load HTSlib/1.3.1-foss-2016b
module load SAMtools/1.3.1-foss-2016b
module load sambamba/0.6.6-foss-2016b

# run the executable
# A script to map PacBio data with SV detection optimised. Designed for the Phoenix supercomputer
# Set common paths
exeDir=/data/neurogenetics/executables/ngmlr/ngmlr-0.2.6/

usage()
{
echo "# A script to map PacBio data with SV detection optimised. Designed for the Phoenix supercomputer
# Requires: samtools, sambamba, ngmlr
# This script assumes your sequence files are fastq converted using the smrtlink bam2fastq tool.  
#
# Usage sbatch $0 -p Prefix [ -s /path/to/sequence/files -o /path/to/output -g /path/to/genome.fa.gz ] | [ -h | --help ]
#
# Options
# -p	REQUIRED. Prefix for file name. Can be any string of text without spaces or reserved special characters.
# -s	OPTIONAL. /path/to/sequence/files . Default $FASTDIR/fastq
# -g	OPTIONAL. Genome to map to. Default /data/neurogenetics/RefSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/PacBio/\$outPrefix is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you might be reading this too!
#
# 
# Original:  Mark Corbett, 04/01/2018, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 
"
}
## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-p )			shift
					outPrefix=$1
					;;
		-s )			shift
					seqDir=$1
					;;
		-g )			shift
					genome=$1
					shortGenome=$(basename $genome)
					;;
		-o )			shift
					workDir=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done

if [ -z "$outPrefix" ]; then # If no output prefix is specified
	usage
	echo "#ERROR: You need to specify a prefix for your files this can be any text without spaces or reserved special characters.
	# eg. -p myFileName"
	exit 1
fi
if [ -z "$seqDir" ]; then # If path to sequences not specified then do not proceed
	seqDir=$FASTDIR/fastq
	echo "#INFO: Using $seqDir to locate files"
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$FASTDIR/PacBio/$outPrefix
	echo "#INFO: Using $workDir as the output directory"
fi
if [ -z "$genome" ]; then # If genome not specified then use hg38
	genome=/data/neurogenetics/RefSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	shortGenome=GRCh38.dna.primary_assembly.fa
	echo "#INFO: Using $genome as the default genome build"
fi

tmpDir=$FASTDIR/tmp/$outPrefix # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d $tmpDir ]; then
	mkdir -p $tmpDir
fi
if [ ! -d $workDir ]; then
	mkdir -p $workDir
fi

cat $seqDir/*.fastq.gz > $tmpDir/$outPrefix.reads.fastq.gz

$exeDir/ngmlr --bam-fix --no-progress -t 16 -x pacbio -r $genome -q $tmpDir/$outPrefix.reads.fastq.gz |\
samtools view -@ 16 -bT $genome - |\
sambamba sort -l 5 -m 4G -t 16 --tmpdir=$tmpDir -o $workDir/$outPrefix.sort.ngmlr.$shortGenome.bam /dev/stdin
sambamba index $workDir/$outPrefix.sort.ngmlr.$shortGenome.bam
# Clean up
rm -r $tmpDir

