#!/bin/bash

#SBATCH -J mm2ont-cDNA
#SBATCH -o /fast/users/%u/log/mm2ont-DNA-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=05:00:00
#SBATCH --mem=32GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Modules needed
modSAMtools="SAMtools/1.9-foss-2016b"
modHTSlib="HTSlib/1.9-foss-2016b"

# Hard coded paths and variables
minimapProg="/data/neurogenetics/executables/minimap2-2.17_x64-linux/minimap2"
genomeBuild="/data/neurogenetics/RefSeq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
cores=8 # Set the same as above for -n

usage()
{
echo "# Script for mapping Oxford Nanopore reads to the human genome.
# This script bypasses the mm2.ONT.DNA.wdl workflow and just submits the job directly.  The script sets up all of the required inputs using the information 
# submitted via the flags below or default options provided.
# REQUIREMENTS: As a minimum you need the fastq_pass folder and the final_summary_xxx.txt file from your nanopore run.
#
# Usage sbatch --array 0-(\"n - 1 barcodes\") $0 -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -s	REQUIRED. Path to the folder containing the fastq_pass folder.  Your final_summary_xxx.txt must be in this folder.
# -S	OPTIONAL (with caveats). Sample name which will go into the BAM header. If not specified, then it will be fetched 
#       from the final_summary_xxx.txt file.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/ONT/DNA/\$sampleName is used)
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line, eg. MySeqProject20200202-PalindromicDatesRule). 
#                 Default \"SQK-LSK109_\$protocol_group_id\"
# -I	OPTIONAL. Unique ID for the sequence (to go into the @RG line). If not specified the script will make one up.
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Written by Mark Corbett, 28/04/2020
# Modified: (Date; Name; Description)
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-s )			shift
					seqPath=$1
					;;
		-S )			shift
					sampleName=$1
					;;
		-o )			shift
					workDir=$1
					;;
		-L )			shift
					LB=$1
					;;
		-I )			shift
					ID=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done
if [ -z "$seqPath" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "## ERROR: You need to specify the path to the folder containing your fastq_pass folder. Don't include fastq_path in this name."
	exit 1
fi
if [ ! -d "$seqPath/fastq_pass" ]; then # If the fastq_pass directory does not exist then do not proceed
    usage
    echo "## ERROR: The fastq_pass directory needs to be in $seqPath. Don't include fastq_path in this name."
	exit 1
fi

# Assume final_summary file exists for now
finalSummaryFile=$(find $seqPath/final_summary_*)

if [ -z "$sampleName" ]; then # If sample name not specified then look for the final_summary_xxx.txt file or die
	if [ -f "$finalSummaryFile" ]; then
		sampleName=$(grep sample_id $finalSummaryFile | cut -f2 -d"=")
		echo "## INFO: Using sample name $sampleName from $finalSummaryFile"
	else
	    usage
	    echo "## ERROR: No sample name supplied and final_summary_*.txt file is not available. I need at least one of these to proceed"
		exit 1
	fi
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$FASTDIR/ONT/DNA/$sampleName
	echo "## INFO: Using $workDir as the output directory"
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi
if [ -z "$LB" ]; then # If library not specified try to make a specfic one or use "SQK-DCS109" as default
    if [ -f "$finalSummaryFile" ]; then
		protocol_group_id=$(grep protocol_group_id $finalSummaryFile | cut -f2 -d"=") 
		LB="SQK-LSK109-$protocol_group_id"
	else 
	    LB="SQK-LSK109"
	fi
fi
echo "## INFO: Using $LB for library name"

# Collate sequence files if not already done.
# You could just scatter each file as an array job to an alignment but potentially this will be slower by loading up the queue
if [ ! -f "$seqPath/$sampleName.fastq.gz" ]; then
    cd $seqPath/fastq_pass
	fileType=$(ls | head -n1 | rev | cut -d"." -f1 | rev)
	    case $fileType in
            gz)    cat *.gz > ../$sampleName.fastq.gz
			       ;;
		    fastq) cat *.fastq > ../$sampleName.fastq
	               cd ../
			       gzip $sampleName.fastq
				   ;;
			fq)    cat *.fastq > ../$sampleName.fastq
	               cd ../
			       gzip $sampleName.fastq
				   ;;
			*)     usage
			       echo "## ERROR: There doesn't appear to be fastq sequence files in $seqPath/fastq_pass. This script checks for .gz, .fastq and .fq file types.  The first file type found was $fileType"
				   exit 1
		esac
else
    echo "## WARN: A fastq file $seqPath/$sampleName.fastq.gz already exists so I'm going to use it.  If this isn't what you wanted you'll need to remove or move this file before you run this workflow again."
fi

seqFile=$seqPath/$sampleName.fastq.gz

if [ -z "$ID" ]; then # If no ID then fetch from the .fastq file
	ID=$(zcat $seqFile | head -n 1 | tr " " "\n" | grep runid | cut -f2 -d"=")
fi
echo "## INFO: Using $ID for the sequence ID"

# Build the input .json file
cd $workDir

## Run the script ##
module load $modSAMtools
module load $modHTSlib
${minimapProg} -ax map-ont \
-R "@RG\\tID:${ID}\\tLB:${LB}\\tPL:ONT\\tSM:${sampleName}" \
-t ${cores} ${genomeBuild} ${seqFile} |\
samtools view -bT ${genomeBuild} - |\
samtools sort -l 5 -m 4G -@${cores} -T${sampleName} -o ${workDir}/${sampleName}.sort.bam -
samtools index ${workDir}/${sampleName}.sort.bam
