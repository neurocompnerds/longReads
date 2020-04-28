#!/bin/bash

#SBATCH -J mm2ont-cDNA-wdl-sub
#SBATCH -o /fast/users/%u/log/mm2ont-cDNA-wdl-sub-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=05:05:00
#SBATCH --mem=3GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Modules needed
modJava="Java/9.0.4"
modSAMtools="SAMtools/1.9-foss-2016b"
modHTSlib="HTSlib/1.9-foss-2016b"

# Hard coded paths and variables
cromwellPath="/data/neurogenetics/executables/cromwell"
cromwellJar="cromwell-48.jar"
minimapProg="/data/neurogenetics/executables/minimap2-2.17_x64-linux/minimap2"
genomeBuild="/data/neurogenetics/RefSeq/GATK/hg38/Homo_sapiens_assembly38.fasta"
scriptDir="/data/neurogenetics/git/PhoenixScripts/mark/longReads/ONT"
cores=8

usage()
{
echo "# Script for mapping Oxford Nanopore cDNA reads to the human genome.
# This script coordinates submission of the mm2.ONT.cDNA.wdl workflow to phoenix.  The script sets up all of the required inputs using the information 
# submitted via the flags below or default options provided.
# REQUIREMENTS: As a minimum you need the fastq_pass folder and the final_summary_xxx.txt file from your nanopore run.
#
# Usage sbatch $0 -s /path/to/sequences -o /path/to/output -c /path/to/config.cfg -S SAMPLE -L LIBRARY -I readGroupID] | [ - h | --help ]
#
# Options
# -s	REQUIRED. Path to the folder containing the fastq_pass folder.  Your final_summary_xxx.txt must be in this folder.
# -S	OPTIONAL (with caveats). Sample name which will go into the BAM header. If not specified, then it will be fetched 
#       from the final_summary_xxx.txt file.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/ONT/cDNA/\$sampleName is used)
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line, eg. MySeqProject20200202-PalindromicDatesRule). 
#                 Default \"SQK-DCS109_\$protocol_group_id\"
# -I	OPTIONAL. Unique ID for the sequence (to go into the @RG line). If not specified the script will make one up.
# -h or --help	Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Written by Mark Corbett, 21/02/2020
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
					outputDir=$1
					;;
		-L )			shift
					LB=$1
					;;
		-I )			shift
					readGroupID=$1
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
	echo "## ERROR: You need to specify the path to your fastq_pass folder"
	exit 1
fi
if [ ! -d "$seqPath/fastq_pass" ]; then # If the fastq_pass directory does not exist then do not proceed
    usage
    echo "## ERROR: The fastq_pass directory needs to be in $seqPath"
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
if [ -z "$outputDir" ]; then # If no output directory then use default directory
	outputDir=$FASTDIR/ONT/cDNA/$sampleName
	echo "## INFO: Using $FASTDIR/ONT/cDNA/$sampleName as the output directory"
fi
if [ ! -d "$outputDir" ]; then
	mkdir -p $outputDir
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

if [ -z "$readGroupID" ]; then # If no readGroupID then fetch from the .fastq file
	readGroupID=$(zcat $seqFile | head -n 1 | tr " " "\n" | grep runid | cut -f2 -d"=")
fi
echo "## INFO: Using $readGroupID for the sequence ID"

# Build the input .json file
cd $outputDir
echo "{
  \"minimap2_ONT_cDNA.mimimap2.samtools\": \"$modSAMtools\",
  \"minimap2_ONT_cDNA.mimimap2.platform\": \"ONT\",
  \"minimap2_ONT_cDNA.mimimap2.htslib\": \"$modHTSlib\",
  \"minimap2_ONT_cDNA.mimimap2.seqFile\": \"$seqFile\",
  \"minimap2_ONT_cDNA.mimimap2.program\": \"$minimapProg\",
  \"minimap2_ONT_cDNA.mimimap2.cores\": \"$cores\",
  \"minimap2_ONT_cDNA.mimimap2.LB\": \"$LB\",
  \"minimap2_ONT_cDNA.mimimap2.sampleName\": \"$sampleName\",
  \"minimap2_ONT_cDNA.mimimap2.genomeBuild\": \"$genomeBuild\",
  \"minimap2_ONT_cDNA.mimimap2.outputDir\": \"$outputDir\",
  \"minimap2_ONT_cDNA.mimimap2.readGroupID\": \"$readGroupID\"
}
" > $outputDir/$sampleName.inputs.json

## Submit the workflow to the queue ##
module load $modJava
java -Dconfig.file=$scriptDir/cromwell_slurm.conf \
-jar $cromwellPath/$cromwellJar \
run $scriptDir/mm2.ONT.DNA.wdl \
--inputs $outputDir/$sampleName.inputs.json
