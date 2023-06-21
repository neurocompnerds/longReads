#!/bin/bash
#SBATCH -J canu.launch
#SBATCH -o /hpcfs/users/%u/log/canu.launch-slurm-%j.out
#SBATCH -p skylake,icelake,skylakehm,v100cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=3-00:00:00
#SBATCH --mem=1GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# Set common paths and defaults
userDir="/hpcfs/users/${USER}"

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
# sbatch $0 -p Prefix -s /path/to/sequence/files [ -o /path/to/output -g genome_size ] | [ -h | --help ]
#
# Options:
# -p	REQUIRED. Prefix ( value = \$outPrefix ) for file names. Can be any string of text without spaces or reserved special characters.
# -s	REQUIRED. /path/to/sequence/files The top level directory that will contain either PacBio/ or ONT/ folders or both depending on what data you have.
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $userDir/Canu/\$outPrefix is used)
# -g    OPTIONAL. Genome size. Default 3.2g for human (mouse is 2.7g).
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
					seqDir=$1
					;;
		-g )			shift
					genomeSize=$1
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
if [ -z "$seqDir" ]; then # If path to sequences not specified then fail
    usage
	echo "## ERROR: You need to tell me where to locate your fastq files"
	exit 1
fi
if [ -z "$workDir" ]; then # If no output directory then use default location
	workDir=$userDir/Canu/$outPrefix
	echo "#INFO: Using $workDir as the output directory"
fi
if [ -z "$genomeSize" ]; then # If genome not specified assume human
	genomeSize="3.2g"
	echo "#INFO: Using $genomeSize as the genome size which corresponds to human."
fi

if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi

if [ -d "$seqDir/PacBio" ]; then
	pbReads="-pacbio $seqDir/PacBio/*.fastq.gz"
fi
if [ -d "$seqDir/ONT" ]; then
	ontReads="-nanopore $seqDir/ONT/*.fastq.gz"
fi
#
#                     O_
#               \-----\/---/
# Launch Canu  ~~\~~~~/~~~/~~
module purge
module use /apps/skl/modules/all
module load canu/2.1.1

# Start paddling
canu -p $outPrefix -d $workDir genomeSize=$genomeSize \
gridOptionsJobName="canu" \
gridOptions="-N 1 -p skylake,icelake,skylakehm,v100cpu --time=72:00:00 --mail-type=END --mail-type=FAIL --mail-user=$USER@adelaide.edu.au" \
merylMemory=125 \
gridOptionsBAT=" -o $userDir/log/canu.slurm-%j.out -p highmem --time=72:00:00 --mail-type=END --mail-type=FAIL --mail-user=$USER@adelaide.edu.au" \
$pbReads $ontReads

