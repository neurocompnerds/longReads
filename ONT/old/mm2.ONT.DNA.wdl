workflow minimap2_ONT_cDNA {
    call mimimap2
}
task mimimap2 {
    String htslib
    String samtools
    String program
    String genomeBuild
    String buildID
    String seqFile
    String readGroupID
    String sampleName
    String platform
    String LB
    String outputDir
    Int cores
    command {
        module load ${htslib}
        module load ${samtools}
        ${program} -ax map-ont \
        -R "@RG\\tID:${readGroupID}\\tLB:${LB}\\tPL:${platform}\\tSM:${sampleName}" \
        -t ${cores} \
        ${genomeBuild} ${seqFile} |\
        samtools view -bT ${genomeBuild} - |\
        samtools sort -l 5 -m 4G -@${cores} -T${sampleName} -o ${outputDir}/${sampleName}.sort.${buildID}.bam -
        samtools index ${outputDir}/${sampleName}.sort.${buildID}.bam
    }
    output {
        File sortedBAM = "${outputDir}/${sampleName}.sort.${buildID}.bam"
        File sortedBAI = "${outputDir}/${sampleName}.sort.${buildID}.bam.bai"
    }
    runtime {
        job_title: "mm2ont-DNA"
        slurm_out: "/hpcfs/users/%u/log/mm2ont-DNA-slurm-%j.out"
        nodes: "1"
        cores: "9"
        queue: "batch"
        requested_memory_mb_per_core: "4"
        time: "05:00:00"
    }
}