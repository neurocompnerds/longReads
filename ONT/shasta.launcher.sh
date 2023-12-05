#!/bin/bash

# Set fixed variables:
shastaScript="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/mark/longReads/ONT/shasta.sh"
gunzipScript="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call/utilities/bulkUnzip.sh"
gzipScript="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call/utilities/bulkGzip.sh"
userDir="/hpcfs/users/${USER}"
CONFIG=Nanopore-May2022
logDir="/hpcfs/users/${USER}/log"
if [ ! -d "$logDir" ]; then
    mkdir -p $logDir
fi

usage()
{
echo "# Script that coordinates tasks for assembling Oxford Nanopore reads using Shasta.
# This script accepts compressed Nanopore FASTQ files as input but it will temporarily unzip them, this might affect other programs.
# REQUIREMENTS: As a minimum you need the fastq_pass or pass folder, or your reads in a single fastq or fasta file from your nanopore run.
#
# Usage: $0 -s /path/to/sequences | /path/to/fastq_file [-o /path/to/output -S SAMPLE -c config] | [ - h | --help ]
#
# Options
# -s	REQUIRED. Path to the folder containing the fastq_pass or pass folder OR to a fastq or fasta file.  Your final_summary_xxx.txt should be in this folder to automatically specify sample name.
# -S	OPTIONAL (with caveats). Sample name which will go into the file names header. If not specified, then it will be fetched 
#                from the final_summary_xxx.txt file.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $userDir/ONT/DNA/\$sampleName is used)
# -c    OPTIONAL. Config to use run $SHASTA --command listConfigurations to choose or the default is Nanopore-May2022, which may not be best.
# -h or --help	  Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
# 
# Original: Written by Nandini Sandran
# Modified: (Date; Name; Description)
# 25/09/2023; Mark Corbett; Implement pre run tests for options and test for gz vs uncompressed sequence files.
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        seqPath=$1
                                        ;;
                -o )                    shift
                                        outDir=$1
                                        ;;
                -S )                    shift
                                        sampleName=$1
                                        ;;
                -c )                    shift
                                        CONFIG=$1
                                        ;;
                -h | --help )           shift
                                        usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

# Assume $seqPath and the final_summary file exists for now
finalSummaryFile=$(find $seqPath/final_summary_*)
if [ ! -f "$finalSummaryFile" ]; then # If SeqPath was a fastq file then check the folder containing that fastq for the final_summary_* file
    finalSummaryFile=$(find $(dirname ${seqPath})/final_summary_*)
fi

if [ -z "${sampleName}" ]; then # If sample name not specified then look for the final_summary_xxx.txt file or die
    if [ -f "$finalSummaryFile" ]; then
        sampleName=$(grep sample_id $finalSummaryFile | cut -f2 -d"=")
        echo "## INFO: Using sample name $sampleName from $finalSummaryFile"
    else
        usage
        echo "## ERROR: No sample name supplied and final_summary_*.txt file is not available. I need at least one of these to proceed"
        echo "## INFO: The final summary file looked for was: $finalSummaryFile"
        exit 1
    fi
fi

if [ -z "${seqPath}" ]; then # If path to sequences not specified then die complaining about it
    usage
    echo "## ERROR: You need to specify the path to the folder containing your fastq_pass or pass folder or a fastq or fasta file directly."
    exit 1
fi

# Check if $seqPath is a file or a folder
# If it is a folder then assume this is a nanopore run folder and find the fastq_pass or pass folder then concatenate fastq files within

if [ ! -f  "$seqPath" ]; then
    if [ ! -d "$seqPath/fastq_pass" ]; then # If the fastq_pass directory does not exist then see if it is just called pass or do not proceed
        if [ ! -d "$seqPath/pass" ]; then
            usage
            echo "## ERROR: The fastq_pass or pass directory needs to be in $seqPath. Don't include fastq_pass or pass in this path."
            exit 1
        else
            fqDir="pass"
        fi
    else
        fqDir="fastq_pass"
    fi
    if [ ! -f  "${seqPath}/${sampleName}.fastq.gz" ]; then
        cat ${seqPath}/${fqDir}/*.gz > ${seqPath}/${sampleName}.fastq.gz
        seqPath=${seqPath}/${sampleName}.fastq.gz
    else
        seqPath=${seqPath}/${sampleName}.fastq.gz
        echo "## WARN: Hey I already found you have a ${sampleName}.fastq.gz file so I'm going to use it, if that isn't what you wanted then you might need to move or delete that file."
    fi
fi

# Is the file gzipped?
fileType=$(file ${seqPath} | grep gzip)
if [ -z "${fileType}" ]; then
    fileIsZipped=FALSE
else
    fileIsZipped=TRUE
    echo "${seqPath}" > file.to.unzip.for.shasta.txt
fi

# Find or create the output directory
if [ -z "$outDir" ]; then # If no output directory then use default directory
    outDir=$userDir/assemblies/ONT/${sampleName}
    echo "## INFO: Using $outDir as the output directory"
fi

if [ ! -d "$outDir" ]; then
    mkdir -p $outDir
fi

if [ "${fileIsZipped}" = TRUE ]; then
    seqFile=$(basename ${seqPath} .gz)
    seqFolder=$(dirname ${seqPath})
    UnzipJob=`sbatch $gunzipScript -i file.to.unzip.for.shasta.txt`
    UnzipJob=$(echo $UnzipJob | cut -d" " -f4)
    seqPath=${seqFolder}/${seqFile}
    echo "${seqFolder}/${seqFile}" > file.to.gzip.txt
    shastaJob=`sbatch --dependency=afterok:${UnzipJob} ${shastaScript} -s ${seqPath} -o ${outDir} -S ${sampleName} -c ${CONFIG}`
    shastaJob=$(echo $shastaJob | cut -d" " -f4)
    sbatch --dependency=afterok:${shastaJob} $gzipScript -i file.to.gzip.txt
else
    sbatch ${shastaScript} -s ${seqPath} -o ${outDir} -S ${sampleName} -c ${CONFIG}
    echo "## INFO: Please consider storing your sequence files in gzipped format.  This launcher script is written to handle decompressing reads for shasta then compressing again when everything is finished."
fi
