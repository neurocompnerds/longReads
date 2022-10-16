#!/bin/bash

#SBATCH -J LAST-gDNA
#SBATCH -o /hpcfs/users/%u/log/LAST-gDNA-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=10:00:00
#SBATCH --mem=48GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Modules needed
modList=("arch/haswell" "SAMtools/1.9-foss-2016b" "HTSlib/1.9-foss-2016b" "R" "Python/3.7.0")

# Hard coded paths and variables
userDir="/hpcfs/users/${USER}"
lastProgDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/last/bin/"
tandemGenotypesProgDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/tandem-genotypes"
cores=8 # Set the same as above for -n

# Genome list (alter case statement below to add new options)
set_genome_build() {
case "${buildID}" in
    GRCh38 )    genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/LAST/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                ;;
    hs37d5 )    genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/LAST/hs37d5/hs37d5.fa.gz"
                ;;
    GRCm38 | mm10 ) buildID="GRCm38"   
                    genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GRCm38_68.fa"
                ;;
    T2T_CHM13v2 )   genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/T2T_CHM13v2.0.ucsc.ebv.fa.gz"
                ;;
    * )         buildID="GRCh38"
                genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/LAST/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                echo "## WARN: Genome build ${buildID} not recognized, the default genome will be used."
                ;;
esac
}

usage()
{
echo "# Script for mapping Oxford Nanopore reads to the human genome.
#
# REQUIREMENTS: As a minimum you need the fastq_pass folder and the final_summary_xxx.txt file from your nanopore run.
#
# Usage sbatch $0 -s /path/to/sequences -o /path/to/output -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
# Usage (with barcodes) sbatch --array 0-(n-1 barcodes) $0 -s /path/to/sequences -b -o /path/to/output -S SAMPLE -L LIBRARY -I ID] | [ - h | --help ]
#
# Options
# -s	REQUIRED. Path to the folder containing the fastq_pass folder.  Your final_summary_xxx.txt must be in this folder.
# -b    DEPENDS.  If you used barcodes set the -b flag.  If you want meaningful sample ID add a file called barcodes.txt to the sequence folder with the 
#                 tab delimited barcode and ID on each line.
# -g    OPTIONAL. Genome build to use, select from either GRCh38, hs37d5, T2T_CHM13v2 or GRCm38. Default is GCA_000001405.15_GRCh38_no_alt_analysis_set
# -S	OPTIONAL.(with caveats). Sample name which will go into the BAM header. If not specified, then it will be fetched 
#                 from the final_summary_xxx.txt file.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $userDir/ONT/DNA/\$sampleName is used)
# -L	OPTIONAL. Identifier for the sequence library (to go into the @RG line, eg. MySeqProject20200202-PalindromicDatesRule). 
#                 Default \"SQK-LSK110_\$protocol_group_id\"
# -I	OPTIONAL. Unique ID for the sequence (to go into the @RG line). If not specified the script will make one up.
# -h or --help	  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Written by Mark Corbett, 01/09/2020
# Modified: (Date; Name; Description)
# 03/02/2022; Mark; Update module loading. Add genome build selection.
# 11/10/2022; Mark Corbett; Add buildID to .bam file name. Add T2T_CHM13v2 to genome list
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-s )			shift
					seqPath=$1
					;;
		-b )			shift
					barcodes=true
					;;
		-S )			shift
					sampleName=$1
					;;
        -g )            shift
                        buildID=$1
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
if [ ! -d "$seqPath/pass" ]; then # If the fastq_pass directory does not exist then do not proceed
    usage
    echo "## ERROR: The fastq_pass directory needs to be in $seqPath. Don't include fastq_pass in this name."
	exit 1
fi

# Assume final_summary file exists for now
finalSummaryFile=$(find $seqPath/final_summary_*)

if "$barcodes"; then
    if [ -f "$seqPath/barcodes.txt" ]; then
        echo "## INFO: Found barcodes.txt file"
		BC=($(cut -f1 $seqPath/barcodes.txt))
		sampleName=($(cut -f2 $seqPath/barcodes.txt))
	else
	    BC=($(ls $seqPath/pass/bar*))
		sampleName=($(ls $seqPath/pass/bar*))
		echo "## INFO: Using generic barcodes as sample names (if you supply a barcodes.txt file you can specify sample names)."
	fi
fi
        
if [ -z "${sampleName[$SLURM_ARRAY_TASK_ID]}" ]; then # If sample name not specified then look for the final_summary_xxx.txt file or die
	if [ -f "$finalSummaryFile" ]; then
		sampleName=$(grep sample_id $finalSummaryFile | cut -f2 -d"=")
		echo "## INFO: Using sample name $sampleName from $finalSummaryFile"
	else
	    usage
	    echo "## ERROR: No sample name supplied and final_summary_*.txt file is not available. I need at least one of these to proceed"
		exit 1
	fi
fi

# Set the genome build using the function defined above.
set_genome_build
echo "## INFO: Using the following genome build: $genomeBuild"

if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$userDir/ONT/DNA/${sampleName[$SLURM_ARRAY_TASK_ID]}
	echo "## INFO: Using $workDir as the output directory"
fi

if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi

if [ -z "$LB" ]; then # If library not specified try to make a specfic one or use "SQK-DCS109" as default
    if [ -f "$finalSummaryFile" ]; then
		protocol_group_id=$(grep protocol_group_id $finalSummaryFile | cut -f2 -d"=") 
		LB="SQK-LSK110-$protocol_group_id"
	else 
	    LB="SQK-LSK110"
	fi
fi
echo "## INFO: Using $LB for library name"

# Collate sequence files if not already done.
# You could just scatter each file as an array job to an alignment but potentially this will be slower by loading up the queue

if [ ! -f  "$seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz" ]; then
    cd $seqPath/pass/${BC[$SLURM_ARRAY_TASK_ID]}
	fileType=$(ls | head -n1 | rev | cut -d"." -f1 | rev)
	    case $fileType in
            gz)    cat *.gz > $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz
			       ;;
		    fastq) cat *.fastq > $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq
	               cd ../
			       gzip $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq
				   ;;
			fq)    cat *.fastq > $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq
	               cd ../
			       gzip $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq
				   ;;
			*)     usage
			       echo "## ERROR: There doesn't appear to be fastq sequence files in $seqPath/fastq_pass/${BC[$SLURM_ARRAY_TASK_ID]}. This script checks for .gz, .fastq and .fq file types.  The first file type found was $fileType"
				   exit 1
		esac
else
    echo "## WARN: A fastq file $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz already exists so I'm going to use it.
	               If this isn't what you wanted you'll need to remove or move this file before you run this workflow again."
fi

if [ -z "$ID" ]; then # If no ID then fetch from the .fastq file
	ID=$(zcat $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz | head -n 1 | tr " " "\n" | grep runid | cut -f2 -d"=")
fi
echo "## INFO: Using $ID for the sequence ID"

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Run the script ##
cd $workDir
${lastProgDir}/last-train -P${cores} -Q0 ${genomeBuild} $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz > $workDir/${sampleName[$SLURM_ARRAY_TASK_ID]}.par
${lastProgDir}/lastal -P${cores} -p $workDir/${sampleName[$SLURM_ARRAY_TASK_ID]}.par ${genomeBuild} $seqPath/${sampleName[$SLURM_ARRAY_TASK_ID]}.fastq.gz | \
${lastProgDir}/last-split -fMAF | gzip > $workDir/${sampleName[$SLURM_ARRAY_TASK_ID]}.maf.gz
