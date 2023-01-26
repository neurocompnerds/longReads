#!/bin/bash

#SBATCH -J sniffles
#SBATCH -o /hpcfs/users/%u/log/sniffles-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --time=00:30:00
#SBATCH --mem=8GB

# Notification configuration 
#SBATCH --mail-type=END                                         
#SBATCH --mail-type=FAIL                                        
#SBATCH --mail-user=%u@adelaide.edu.au

# run the executable
# A script to structural variants from pacbio aligned with ngmlr. Designed for the Phoenix supercomputer
# Set common paths and defaults
userDir="/hpcfs/users/${USER}"
neuroDir="/hpcfs/groups/phoenix-hpc-neurogenetics"
refDir="$neuroDir/RefSeq"
customModDir="$neuroDir/executables/easybuild/modules/all"
modList=("arch/skylake" "SAMtools/1.12" "Anaconda3/2020.07")

supportReads=10
usage()
{
echo "# A script to structural variants from pacbio aligned with ngmlr. Designed for the Phoenix supercomputer
# Requires: sniffles installed into a conda environment called sniffles
# This script assumes your sequences were mapped with ngmlr or BWA-MEM with the -M option (see sniffles docs https://github.com/fritzsedlazeck/Sniffles/wiki).  
#
# Usage sbatch $0 -b /path/to/bamFile.bam[.cram] [ -g /path/to/genome.fa[.gz] -o /path/to/output ] | [ -h | --help ]
#
# Options
# -b <arg>  REQUIRED. Path to your .bam or .cram file
# -g <arg>  OPTIONAL: Path to the original reference that your BAM file was mapped to. The script will try to locate the right genome based on the @SQ lines in the bam header if you don't set this.
# -o <arg>  OPTIONAL. Path to where you want to find your file output (if not specified the bam folder is used)
# -h or --help	Prints this message.  Or if you got one of the options above wrong you might be reading this too!
#
# 
# Original:  Mark Corbett, 04/01/2018, mark dot corbett at adelaide.edu.au
# Modified: (Date; Name; Description)
# 25/01/2023; Mark; Update to Sniffles2
# 
"
}

select_genome_build()
{
case "${genomeSize}" in
    3099922541 )    buildID="GRCh38"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
                    ;;
    3217346917 )    buildID="hs38DH"
                    genomeBuild="$refDir/hs38DH.fa"
                    ;;
    3137454505 )    buildID="hs37d5"
                    genomeBuild="$refDir/hs37d5.fa.gz"
                    ;;
    2730871774 )    buildID="GRCm38"   
                    genomeBuild="$refDir/GRCm38_68.fa"
                    ;;
    3117463893 )    buildID="CHM13v2"
                    genomeBuild="$refDir/T2T_CHM13v2.0.ucsc.ebv.fa.gz"
                    ;;
    3137161264 )    buildID="hg19"
                    genomeBuild="$refDir/ucsc.hg19.fasta"
                    ;;
    3105715063 )    buildID="GRCh38.hs38d1"
                    genomeBuild="$refDir/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz"
                    ;;
    3099750718 )    buildID="GRCh38"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                    ;;
    3031042417 )    buildID="GRCh38.blacklist"
                    genomeBuild="$refDir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.r101.s501.blacklist.gz"
                    ;;
    3101804741 )    buildID="hg19_1stM_unmask_ran_all"
                    genomeBuild="$refDir/hg19_1stM_unmask_ran_all.fa"
                    ;;
    * )         usage
                echo "## ERROR: Genome length $genomeSize was not matched, you may need to specify the genome build directly."
                exit 1
                ;;
esac
}

## Set Variables ##
while [ "$1" != "" ]; do
	case $1 in
		-b )		shift
					bamFile=$1
					bamDir=$(dirname $bamFile)
					outPrefix=$(basename -s .bam $bamFile)
					;;
        -g )        shift
                    genomeBuild=$1
                    ;;
		-o )		shift
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

if [ -z "$bamFile" ]; then # If the bam or cram file isn't specified then give up and be grumpy about it
	usage
	echo "#ERROR: You need to tell me which bam or cram file you want to run this on.
	# eg. -b /path/to/bamFile.bam[.cram]"
	exit 1
fi

fullBamFile=$( basename ${bamFile} )
baseBamFile=${fullBamFile%.*}

## Load modules ##
module use $customModDir
for mod in "${modList[@]}"; do
    module load $mod
done

if [ -z "$genomeBuild" ]; then # If genome not specified then check the bam file to work it out
    genomeSize=$(samtools view -H ${bamFile} | grep @SQ | cut -f3 | cut -f2 -d":" | awk '{s+=$1} END {printf "%.0f\n", s}' -)
    select_genome_build
fi

if [ -z "$workDir" ]; then # If no output directory then use default directory
	workDir=$userDir/variants/SV/
	echo "#INFO: Using $workDir as the output directory"
fi

if [ ! -d $workDir ]; then
	mkdir -p $workDir
fi

conda activate sniffles

sniffles --input $bamFile --vcf $workDir/$baseBamFile.vcf.gz --reference $genomeBuild

conda deactivate
