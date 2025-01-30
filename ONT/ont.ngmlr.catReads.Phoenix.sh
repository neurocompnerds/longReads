#!/bin/bash

#SBATCH -J ngmlr
#SBATCH -o /hpcfs/users/%u/log/ngmlr.ont.slurm-%j.out
#SBATCH -p icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=${USER}@adelaide.edu.au

# run the executable
# A script to map ONT data with SV detection optimised. Designed for the Phoenix supercomputer
# Set common paths
exeDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/ngmlr/ngmlr-0.2.7/"
userDir="/hpcfs/users/${USER}"
sambambaProg="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/sambamba-0.8.2-linux-amd64-static"
# Modules needed
modList=("SAMtools/1.17-GCC-11.2.0" "HTSlib/1.17-GCC-11.2.0")

usage()
{
echo "# A script to map Oxford Nanopore data with SV detection optimised. Designed for the Phoenix supercomputer
# Requires: samtools, sambamba, ngmlr
# This script assumes your sequence files are fastq they may need to be extracted from the ONT fast5 HDF5 archive first.  
#
# Usage sbatch $0 -p Prefix [ -s /path/to/sequence/files -o /path/to/output -g /path/to/genome.fa.gz ] | [ -h | --help ]
#
# Options
# -p	REQUIRED. Prefix for file name. Can be any string of text without spaces or reserved special characters.
# -s	OPTIONAL. /path/to/sequence/files . Default $userDir/sequences/ONT/\$outPrefix/fastq_pass
# -g	OPTIONAL. Genome to map to. Default /hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
# -o	OPTIONAL. Path to where you want to find your file output (if not specified an output directory $userDir/ONT/\$outPrefix is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you might be reading this too!
#
# 
# Original:  Mark Corbett, 04/01/2018, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 26/03/2020; Mark; Update to v0.2.7, Update to use info in the final_report.txt file
# 15/09/2020; Mark Corbett; Update Phoenix paths
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
	echo "## ERROR: You need to specify a prefix for your files this can be any text without spaces or reserved special characters.
	# eg. -p myFileName"
	exit 1
fi
if [ -z "$seqDir" ]; then # If path to sequences not specified then do not proceed
	seqDir=$userDir/sequences/ONT/$outPrefix/fastq_pass
	echo "## WARN: Using $seqDir to locate files, if your files aren't here the script might run but you won't get any alignments"
fi
if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$userDir/ONT/$outPrefix
	echo "## INFO: Using $workDir as the output directory"
fi
if [ -z "$genome" ]; then # If genome not specified then use hg38
	genome="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
	shortGenome="GRCh38_no_alt_analysis_set"
	echo "## INFO: Using $genome as the default genome build"
fi

tmpDir=$userDir/tmp/$SLURM_JOB_ID # Use a tmp directory for the sambamba temp files
if [ ! -d "$tmpDir" ]; then
	mkdir -p $tmpDir
fi
if [ ! -d "$workDir" ]; then
	mkdir -p $workDir
fi

# load modules
for mod in "${modList[@]}"; do
    module load $mod
done

cat $seqDir/*.fastq.gz > $tmpDir/$outPrefix.reads.fastq.gz

$exeDir/ngmlr --bam-fix --no-progress -t 16 -x ont -r $genome -q $tmpDir/$outPrefix.reads.fastq.gz |\
samtools view -@ 16 -bT $genome - |\
$sambambaProg sort -l 5 -m 4G -t 16 --tmpdir=$tmpDir -o $workDir/$outPrefix.sort.ngmlr.ont.$shortGenome.bam /dev/stdin
$sambambaProg index $workDir/$outPrefix.sort.ngmlr.ont.$shortGenome.bam
# Clean up
rm -r $tmpDir
