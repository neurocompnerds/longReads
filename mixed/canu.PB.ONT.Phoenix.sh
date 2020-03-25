#!/bin/bash
# canu.PB.ONT.Phoenix.sh A script to de novo assemble PacBio and Oxford Nanopore sequencing using Canu. Designed for the Phoenix supercomputer

usage()
{
echo "# A script to de novo assemble PacBio and Oxford Nanopore sequencing using Canu. Designed for the Phoenix supercomputer
# Requires: canu
# This script assumes your PacBio sequence files are fastq converted using the smrtlink bam2fastq tool.
# You can have PacBio, Oxford Nanopore or both. Set up your sequence files or links to them in one root directory with subdirectories for PacBio and ONT (case and spelling sensitive)
# ONLY put a PacBio folder in if you have that data and ONLY put an ONT folder in if you have that type data otherwise don't create them. If you have both types of data, make both folders.
# eg. 
# /path/to/my/fastq/PacBio/pacBio.reads.fastq.gz
# /path/to/my/fastq/ONT/ONT.reads.fastq.gz
#
# Usage: 
# $0 -p Prefix [ -s /path/to/sequence/files -o /path/to/output ] | [ -h | --help ]
#
# Options:
# -p	REQUIRED. Prefix ( value = \$outPrefix ) for file names. Can be any string of text without spaces or reserved special characters.
# -s	OPTIONAL. /path/to/sequence/files The top level directory that will contain either PacBio/ or ONT/ folders or both depending on what data you have. Default is $FASTDIR/fastq/Canu/\$outPrefix
#                 If you intend to use the default then you'll need to make that directory manually before running the script.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $FASTDIR/Canu/\$outPrefix is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you might be reading this too!
#
# 
# Original:  Mark Corbett, 01/06/2018, mark dot corbett at adelaide.edu.au
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
					pbDir=$1
					;;
		--ont )			shift
					ontDir=$1
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

# Sort out all the things
if [ -z "$outPrefix" ]; then # If no output prefix is specified then get eaten by a hippo
	usage
	echo "#ERROR: You need to specify a prefix for your files this can be any text without spaces or reserved special characters.
	# eg. -p myFileName"
	exit 1
fi
if [ -z "$seqDir" ]; then # If path to sequences not specified then use default location
	seqDir=$FASTDIR/fastq/Canu/$outPrefix
	echo "#INFO: Using $seqDir to locate files"
fi
if [ -z "$workDir" ]; then # If no output directory then use default location
	workDir=$FASTDIR/Canu/$outPrefix
	echo "#INFO: Using $workDir as the output directory"
fi

if [ ! -d $workDir ]; then
	mkdir -p $workDir
fi

if [ -d $seqDir/PacBio ]; then
	pbReads="-pacbio-raw $seqDir/PacBio/*.fastq.gz"
fi
if [ -d $seqDir/ONT ]; then
	ontReads="-nanopore-raw $seqDir/ONT/*.fastq.gz"
fi
#
#                     O_
#               \-----\/---/
# Launch Canu  ~~\~~~~/~~~/~~
module load canu/1.7-foss-2016b

# Start paddling
canu -p $outPrefix -d $workDir genomeSize=3.6g \
gridOptionsJobName="canu" \
gridOptions=" -A robinson -p batch --time=72:00:00 --mail-type=END --mail-type=FAIL --mail-user=$USER@adelaide.edu.au" \
merylMemory=125 \
gridOptionsBAT=" -o $FASTDIR/launch/canu.slurm-%j.out -A robinson -p highmem --time=72:00:00 --mail-type=END --mail-type=FAIL --mail-user=$USER@adelaide.edu.au" \
$pbReads $ontReads

