#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Pipeline parameters - MUST be at top level, outside any blocks
 */

params.pod5_dir = 'data/pod5/'
params.model = 'rna004_130bps_sup@v5.1.0'
params.modif = 'pseU'
params.reference = 'data/reference/Homo_sapiens.rRNA.fasta'
params.seqtagger_sif = null

process dorado_model {
    storeDir "${workflow.launchDir}/models"

    input: 
    val model

    output:
    path "${model}"

    script:
    """
    dorado download --model ${model}
    """
}

process dorado_basecall {
    input:
    path pod5_dir
    path reference_file
    val model
    val modif
    
    output:
    path "output.bam"
    
    script:
    """
    dorado basecaller \
        ${model} \
        ${pod5_dir} \
        --modified-bases ${modif} \
        --reference ${reference_file} \
        --mm2-opts "-x map-ont -N 0 -k 13" \
        > output.bam \
        2> dorado_basecall.log
    """
}

process sort_index_bam {
    input:
    path bam_file
    
    output:
    path "*.sorted.bam", emit: bam
    path "*.sorted.bam.bai", emit: bai
    
    script:
    """
    samtools sort -o ${bam_file.baseName}.sorted.bam ${bam_file}
    samtools index ${bam_file.baseName}.sorted.bam
    """
}   

process download_seqtagger_image {
    storeDir "${workflow.launchDir}/data"
    
    output:
    path 'seqtagger.sif'
    
    script:
    """
    apptainer pull seqtagger.sif docker://lpryszcz/seqtagger:2.0a
    """
}

process seqtagger_index {
    storeDir "${workflow.launchDir}/data/pod5/demux"
    
    input:
    path sif_file
    
    output:
    path "*.demux.tsv.gz"
    
    script:
    """
    apptainer run --nv \
        --bind ${workflow.launchDir}/data:/data \
        ${sif_file} \
        mRNA -k /opt/app/models/b96_RNA004 -r \
        -i /data/pod5/ \
        -o /data/pod5/demux
    
    cp /data/pod5/demux/*.demux.tsv.gz .
    """
}

process seqtagger_demultiplex {
    publishDir "${workflow.launchDir}/results/demultiplexed_bams", mode: 'copy'
    
    input:
    path bam_file
    path index
    path sif_file
    
    output:
    path "*.bam"
    
    script:
    """
    apptainer run --nv --bind \$(pwd):/work ${sif_file} \
        bam_split_by_barcode.py \
        -i /work/${index} \
        -f /work/${bam_file} \
        -o /work/
    """
}

process pileup {
    publishDir "${workflow.launchDir}/results/pileups", mode: 'copy'
    
    input:
    path demuxed_bam
    path reference_file
    
    output:
    path "*.bed"
    
    script:
    """
    modkit pileup \
        ${demuxed_bam} \
        ${demuxed_bam.baseName}.bed \
        --ref ${reference_file} \
        --filter-threshold 0.90 \
        --threads ${task.cpus}
    """
}

process Rreport {
    publishDir "${workflow.launchDir}/results", mode: 'copy'
    
    input:
    path pileup_files
    
    output:
    path "report.html"
    
    script:
    """
    Rscript ${projectDir}/scripts/generate_report.R \
        --input ${pileup_files} \
        --output report.html
    """
}

// ===== WORKFLOW =====

workflow {
    pod5_dir_ch = channel.fromPath(params.pod5_dir, type: 'dir')
    reference_ch = channel.fromPath(params.reference)

    // Download model if needed
    model_ch = dorado_model(params.model)
    
    // Use existing .sif OR download it
    if (params.seqtagger_sif) {
        sif_ch = channel.fromPath(params.seqtagger_sif)
    } else {
        sif_ch = download_seqtagger_image()
    }
    
    index_ch = seqtagger_index(sif_ch)
    
    // Basecall all pod5 files at once
    bam_ch = dorado_basecall(
        pod5_dir_ch,
        reference_ch,
        model_ch,
        params.modif
    )
    
    // Sort and index
    sorted_ch = sort_index_bam(bam_ch)
    
    // Demultiplex using the index
    demux_ch = seqtagger_demultiplex(
        sorted_ch.bam,
        index_ch,
        sif_ch
    )
    
    // Pileup per demuxed file
    pileup_ch = pileup(
        demux_ch.flatten(),
        reference_ch
    )
    
    // Generate report
    Rreport(pileup_ch.collect())
}
