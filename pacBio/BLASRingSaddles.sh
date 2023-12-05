#!/bin/bash

#SBATCH -J BLASR
#SBATCH -o /hpcfs/users/%u/log/slurm-%j.out
#SBATCH -p skylake,icelake,a100cpu
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=01:00:00
#SBATCH --mem=64GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# load modules
module load blasr/5.0.9-GCC-5.3.0-binutils-2.25

# run the executable
# A script to map PacBio unmapped bam files using BLASR designed for the Phoenix supercomputer

usage()
{
echo "# A script to map PacBio unmapped bam files using BLASR designed for the Phoenix supercomputer
# Requires: BLASRv5.0, samtools, sambamba.
# This script assumes your sequence files are unmapped bams.  You need to supply a list of files in a text file that has extension .fofn
#
# Usage sbatch $0 -p Prefix -f files.fofn -g /data/neurogenetics/RefSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -o /path/to/output | [ -h | --help ]
#
# Options
# -p	REQUIRED. Prefix for file outputs
# -f	REQUIRED. List of files in a text file .fofn
# -g	OPTIONAL. Genome to map to
# -o	OPTIONAL. /path/to/output
#
# 
# Original:  Mark Corbett, 18/12/2017, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 
"
}
## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-p )			shift
					Prefix=$1
					;;
		-f )			shift
					seqFiles=$1
					;;
		-g )			shift
					genome=$1
					;;
		-o )			shift
					WORKDIR=$1
					;;
		-h | --help )		usage
					exit 0
					;;
		* )			usage
					exit 1
	esac
	shift
done

if [ -z "$Prefix" ]; then # If no file prefix specified then do not proceed
	usage
	echo "#ERROR: You need to specify a file prefix (PREFIX) referring to your sequence files eg. PREFIX_R1.fastq.gz."
	exit 1
fi
if [ -z "$seqFiles" ]; then # If path to sequences not specified then do not proceed
	usage
	echo "#ERROR: You need to specify the file that lists your unmapped PacBio bam files"
	exit 1
fi
if [ -z "$genome" ]; then # If sample name not specified then use "OUTPREFIX"
	genome=/data/biohub/Refs/human/hg38_hg20_GRCh38p5/GRCh38.p5.genome.fa
	shortGenome=GRCh38.p5
	echo "Using $genome as the default genome build"
fi
if [ -z "$WORKDIR" ]; then # If no output directory then use default directory
	WORKDIR=$FASTDIR/PacBio/$Prefix
	echo "Using $WORKDIR as the output directory"
fi

tmpDir=$FASTDIR/tmp/$OUTPREFIX # Use a tmp directory for all of the GATK and samtools temp files
if [ ! -d $tmpDir ]; then
	mkdir -p $tmpDir
fi
if [ ! -d $WORKDIR ]; then
	mkdir -p $WORKDIR
fi

## Start of the script ##
cd $WORKDIR
blasr $seqFiles $genome -sam -unaligned $Prefix.unaligned.txt -minMatch 15 -nproc 16 > $WORKDIR/$Prefix.$shortGenome.pacbio.sam

module unload blasr/5.0.9-GCC-5.3.0-binutils-2.25
module load HTSlib/1.3.1-GCC-5.3.0-binutils-2.25 
module load SAMtools/1.3.1-GCC-5.3.0-binutils-2.25

samtools view -bT $genome $WORKDIR/$Prefix.$shortGenome.pacbio.sam |\
samtools sort -m 4G -@16 -T$SAMPLE -o $SAMPLE.samsort.bwa.$BUILD.bam -o $WORKDIR/$Prefix.$shortGenome.pacbio.bam -
samtools index $WORKDIR/$Prefix.$shortGenome.pacbio.bam
