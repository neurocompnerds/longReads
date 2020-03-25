#!/bin/bash

#SBATCH -J sniffles
#SBATCH -o /fast/users/%u/launch/sniffles.slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=08:00:00
#SBATCH --mem=16GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# run the executable
# A script to structural variants from pacbio aligned with ngmlr. Designed for the Phoenix supercomputer
# Set common paths and defaults
exeDir=/data/neurogenetics/executables/sniffles/Sniffles-master/bin/sniffles-core-1.0.7
supportReads=10
usage()
{
echo "# A script to structural variants from pacbio aligned with ngmlr. Designed for the Phoenix supercomputer
# Requires: sniffles
# This script assumes your sequences were mapped with ngmlr or BWA-MEM with the -M option (see sniffles docs https://github.com/fritzsedlazeck/Sniffles/wiki).  
#
# Usage sbatch $0 -b /path/to/bamFile.bam [ -s 10 -o /path/to/output ] | [ -h | --help ]
#
# Options
# -b	REQUIRED. Path to your .bam file
# -s	OPTIONAL. number of supporting reads (default 10)
# -o	OPTIONAL. Path to where you want to find your file output (if not specified the bam folder is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you might be reading this too!
#
# 
# Original:  Mark Corbett, 04/01/2018, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 
"
}

## Load modules
module load HTSlib/1.3.1-foss-2016b

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-b )			shift
					bamFile=$1
					bamDir=$(dirname $bamFile)
					outPrefix=$(basename -s .bam $bamFile)
					;;
		-s )			shift
					supportReads=$1 #Lazy overwrite of the default if this is specified
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

if [ -z "$bamFile" ]; then # If no output prefix is specified
	usage
	echo "#ERROR: You need to tell me which bam file you want to run this on.
	# eg. -b /path/to/bamFile.bam"
	exit 1
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$bamDir
	echo "#INFO: Using $workDir as the output directory"
fi
if [ -z "$genome" ]; then # If genome not specified then use hg38
	genome=/data/neurogenetics/RefSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	shortGenome=GRCh38.dna.primary_assembly.fa
	echo "#INFO: Using $genome as the default genome build"
fi

echo "#INFO: Using $supportReads supporting reads to call variants"

tmpDir=$FASTDIR/tmp/$outPrefix # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d $tmpDir ]; then
	mkdir -p $tmpDir
fi
if [ ! -d $workDir ]; then
	mkdir -p $workDir
fi

$exeDir/sniffles -m $bamFile -s $supportReads -v $workDir/$outPrefix.vcf -r 600 --tmp_file $tmpDir --cluster
bgzip $workDir/$outPrefix.vcf
tabix $workDir/$outPrefix.vcf.gz

# Clean up
rm -r $tmpDir

