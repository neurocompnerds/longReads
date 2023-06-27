#!/bin/bash

#SBATCH -J bam2fastq
#SBATCH -o /fast/users/%u/launch/slurm-%j.out
#SBATCH -A robinson
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --mem=3GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load smrtlink/5.0.1

# run the executable
# A batch script to convert PacBio unmapped bam files to fastq. Designed for the Phoenix supercomputer

usage()
{
echo "# A batch script to convert PacBio unmapped bam files to fastq. Designed for the Phoenix supercomputer
# Requires: SMRTlink v5.0.1.
# This script assumes your sequence files are unmapped bams.  You need to supply a list of files and paths in a text file.
#
# Usage sbatch --array o-(n-1) $0 -f files.txt -o /path/to/output | [ -h | --help ]
#
# Options
# -f	REQUIRED. List of files in a text file .fofn
# -o	OPTIONAL. /path/to/output default $FASTDIR/fastq
#
# 
# Original:  Mark Corbett, 03/01/2018, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 
"
}
## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-f )			shift
					seqFiles=$1
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

if [ -z "$seqFiles" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the file that lists your unmapped PacBio bam files"
	exit 1
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$FASTDIR/fastq
	echo "Using $workDir as the output directory"
fi

tmpDir=$FASTDIR/tmp # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d $tmpDir ]; then
	mkdir -p $tmpDir
fi
if [ ! -d $workDir ]; then
	mkdir -p $workDir
fi

## Define the array
bamFile=(`cat $seqFiles`)

## Do the conversion
cd $workDir
bam2fastq -o $(basename -s .bam ${bamFile[$SLURM_ARRAY_TASK_ID]}) ${bamFile[$SLURM_ARRAY_TASK_ID]} 

