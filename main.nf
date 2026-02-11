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

/*params.bam_dir = 'results/demultiplex_bam/'
*params.sample_sheet = 'data/RNA_barcodes_sampleSheet.csv'
*/

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
    storeDir "${workflow.launchDir}/results/demultiplexed_bams"

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
    publishDir "${workflow.launchDir}/results/sorted_bams", mode: 'copy'
    input:
    path bam_file
    
    output:
    tuple path("*.sorted.bam"), path("*.sorted.bam.bai"), emit: sorted
    
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
    """
}

process seqtagger_demultiplex {
    storeDir "${workflow.launchDir}/results/demultiplexed_bams"
    errorStrategy 'ignore'

    input:
    path bam_file
    path index
    path sif_file

    output:
    path "output.bc_*.bam"

    script:
    """
    apptainer run --nv \
        --bind \$(pwd):/work \
        --bind ${workflow.launchDir}/data:/data \
        ${sif_file} \
        bam_split_by_barcode.py \
        -i /data/pod5/demux/${index} \
        -f /work/${bam_file} \
        -o /work/output; exit 0
    """
}

process pileup {
    publishDir "${workflow.launchDir}/results/pileups", mode: 'copy'

    input:
    tuple path(bam), path(bai)
    each path(reference_file)

    output:
    path "*.bed"

    script:
    """
    modkit pileup ${bam} ${bam.baseName}.bed \
        --ref ${reference_file} \
        -n 40000 \
        --filter-threshold 0.90
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
    Rscript ${workflow.projectDir}/scripts/generate_report.R \
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

    // Basecall all pod5 files
    if (file("${workflow.launchDir}/results/demultiplexed_bams/output.bam").exists()) {
        bam_ch = Channel.fromPath("${workflow.launchDir}/results/demultiplexed_bams/output.bam")
    } else {
        bam_ch = dorado_basecall(pod5_dir_ch, reference_ch, model_ch, params.modif)
    }
        
    // Demultiplex using the index
    demux_ch = seqtagger_demultiplex(
        bam_ch,
        index_ch,
        sif_ch
    )
    
        // Sort and index
    sorted_ch = sort_index_bam(demux_ch.flatten())

    // Pileup per demuxed file
    pileup_ch = pileup(
        sorted_ch,
        reference_ch
    )
    
    // Generate report
    Rreport(pileup_ch.collect())
}