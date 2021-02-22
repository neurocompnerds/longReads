#!/bin/bash -l

#SBATCH -J megaMeth
#SBATCH -o /hpcfs/users/%u/log/megaMeth-slurm-%j.out
#SBATCH -A robinson
#SBATCH -p v100
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --gres=gpu:2
#SBATCH --time=04:00:00
#SBATCH --mem=150GB

# Notification configuration
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=${USER}@adelaide.edu.au

# Hard coded paths and variables
# Module list
modList=("arch/skylake" "arch/haswell" "Anaconda3/2020.07" "HTSlib/1.10.2-foss-2016b" "SAMtools/1.10-foss-2016b")

# Paths
userDir="/hpcfs/users/${USER}"
guppyProgDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/ont-guppy/bin"
configDir="/hpcfs/groups/phoenix-hpc-neurogenetics/executables/rerio/basecall_models"
config="res_dna_r941_min_modbases_5mC_5hmC_v001.cfg"
cores=30 # Set the same as above for -n

# Genome list (alter case statement below to add new options)
set_genome_build() {
case "${buildID}" in
    GRCh38 )    genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi"
                ;;
    hs37d5 )    genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/hs37d5.fa.gz"
                ;;
    * )         genomeBuild="/hpcfs/groups/phoenix-hpc-neurogenetics/RefSeq/GCA_000001405.15_GRCh38_no_alt_analysis_set.mmi"
                echo "## WARN: Genome build ${buildID} not recognized, the default genome will be used."
                ;;
esac
}

usage()
{
echo "# Script for calling methylation from Oxford Nanopore reads.
#
# REQUIREMENTS: As a minimum you need the fast5_pass folder and the final_summary_xxx.txt file from your nanopore run.
#
# Usage: sbatch $0 -s /path/to/fast5_folders -g [GRCh38 | hs37d5] -B | [ - h | --help ]
#
# Options
# -s    REQUIRED.    Path to the folder containing the fast5 folder or folders (don't include fast5_pass)
# -S    CONDITIONAL. Either a sample name must be supplied or the final_summary.txt file is needed in /path/to/fast5_folders
# -o    OPTIONAL.    Path to output.  Default is /path/to/fast5_folders/megalodon_output
# -g    OPTIONAL.    Genome build to use, select from either GRCh38 or hs37d5. Default is GCA_000001405.15_GRCh38_no_alt_analysis_set
# -B    OPTIONAL.    Flag to make an in silico bisulfite conversion of the 5mC data
# -h or --help       Prints this message.  Or if you got one of the options above wrong you'll be reading this too!
#
# Original: Written by Mark Corbett, 15/02/2021
# Modified: (Date; Name; Description)
#
"
}

## Set Variables ##
while [ "$1" != "" ]; do
        case $1 in
                -s )                    shift
                                        seqPath=$1
                                        ;;
                -S )                    shift
                                        sampleName=$1
                                        ;;
                -o )                    shift
                                        workDir=$1
                                        ;;
                -g )                    shift
                                        buildID=$1
                                        ;;
                -B )                    shift
                                        convertToBS="--mod-map-base-conv C T --mod-map-base-conv m C"
                                        echo "## INFO: In silico bisulfite conversion will be done on 5mC calls."
                                        ;;
                -h | --help )           usage
                                        exit 0
                                        ;;
                * )                     usage
                                        exit 1
        esac
        shift
done

if [ -z "$seqPath" ]; then # If path to sequences not specified then do not proceed
        usage
        echo "## ERROR: You need to specify the path to the folder containing your fast5_pass folder. Don't include fast5_pass in $seqPath."
        exit 1
fi
if [ ! -d "$seqPath/fast5_pass" ]; then # If the fastq_pass directory does not exist then do not proceed
    usage
    echo "## ERROR: The fast5_pass directory needs to be in $seqPath. Don't include fast5_pass in $seqPath."
        exit 1
fi

# Set the genome build using the function defined above.
set_genome_build
echo "## INFO: Using the following genome build: $genomeBuild"

# Assume final_summary file exists for now
finalSummaryFile=$(find $seqPath/final_summary_*)

if [ -z "${sampleName}" ]; then # If sample name not specified then look for the final_summary_xxx.txt file or die
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
        workDir=$seqPath/megalodon_output
        echo "## INFO: Using $workDir as the output directory"
fi

if [ ! -d "$workDir" ]; then
        mkdir -p $workDir
fi

## Load modules ##
for mod in "${modList[@]}"; do
    module load $mod
done

## Run the script ##
conda activate megalodon

megalodon $seqPath/fast5_pass \
--guppy-server-path $guppyProgDir/guppy_basecall_server \
--guppy-params "-d $configDir" \
--guppy-config $config \
--devices cuda:all \
--processes $cores \
--outputs basecalls mappings mod_mappings mods per_read_mods \
--mod-motif hm CG 0 $convertToBS \
--reference $genomeBuild \
--mod-map-emulate-bisulfite \
--output-directory $workDir \
--overwrite

conda deactivate

cd $workDir
## Sort the bam files (for some reason --sort-mappings does not work on phoenix)
for f in mappings.bam mod_mappings.5mC.bam mod_mappings.5hmC.bam; do
    samtools sort -@ $cores -m 4G $f > $sampleName.$f && samtools index $sampleName.$f
done
 
## Compress some of the data to save a bit of space
gzip basecalls.fastq &
gzip log.txt &
gzip mappings.summary.txt &
gzip sequencing_summary.txt
wait

