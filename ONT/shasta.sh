#!/bin/sh -l
#SBATCH -J ShastaOrShaman
#SBATCH -o /hpcfs/users/%u/log/shasta.slurm-%j.out
#SBATCH --mem=1024GB
#SBATCH -p highmem
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=10:00:00
#Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# Set fixed variables:
SHASTA="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/shasta-Linux-0.11.1"
gzipScript="/hpcfs/groups/phoenix-hpc-neurogenetics/scripts/git/neurocompnerds/map-n-call-utilities/bulkGzip.sh"
userDir="/hpcfs/users/${USER}"
CONFIG=Nanopore-May2022

usage()
{
echo "# Script for assembling Oxford Nanopore reads.
# This script accepts compressed Nanopore FASTQ files as input but it will temporarily unzip them, this might affect other programs.
# REQUIREMENTS: As a minimum you need the fastq_pass or pass folder or your reads in a single fastq or fasta file from your nanopore run.
#
# Usage: sbatch $0 -s /path/to/sequences | /path/to/fastq_file [-o /path/to/output -S SAMPLE -c config] | [ - h | --help ]
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
                                        OUTDIR=$1
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
if [ ! -f "$finalSummaryFile" ]; then
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
# If it is a folder then find the fastq_pass or pass folder and concatenate fastq files within

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
        echo "## WARN: Hey I already found you have a ${sampleName}.fastq.gz file so I'm going to use it if that isn't what you wanted then you might need to move or delete that file."
    fi
fi

# Is the file gzipped?
fileType=$(file ${seqPath} | grep gzip)
if [ -z "${fileType}" ]; then
    fileIsZipped=FALSE
else
    fileIsZipped=TRUE
    echo "## WARN:  Shasta takes only uncompressed fastq or fasta files as input. It is highly recommended that you use the shasta.launcher.sh to offload compressing and recompressing to other (not highmem) nodes"
fi

# Find or create the output directory
if [ -z "$OUTDIR" ]; then # If no output directory then use default directory
    OUTDIR=$userDir/assemblies/ONT/${sampleName}
    echo "## INFO: Using $OUTDIR as the output directory"
fi

if [ ! -d "$OUTDIR" ]; then
    mkdir -p $OUTDIR
fi

if [ "${fileIsZipped}" = TRUE ]; then
    seqFile=$(basename ${seqPath} .gz)
    seqFolder=$(dirname ${seqPath})
    gunzip ${seqPath}
    seqPath=${seqFolder}/${seqFile}
    echo "${USER} is currently running Shasta on this file which is why it is not gzipped.  Please leave until the run is finished." > ${seqFolder}/${seqFile}.Not.gzipped.txt
fi

${SHASTA} --input ${seqPath} --config ${CONFIG} --assemblyDirectory ${OUTDIR} --Reads.minReadLength 5000 --threads 30

if [ "${fileIsZipped}" = TRUE ]; then
    gzip ${seqPath}
    rm ${seqFolder}/${seqFile}.Not.gzipped.txt
fi
