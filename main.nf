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
params.sample_sheet = 'data/sample_sheet.csv'
params.region = 'data/reference/region.bed'

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
    def sample = bam_file.baseName
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

process change_name {
    publishDir "${workflow.launchDir}/results/demultiplexed_bams", mode: 'copy'

    input:
    tuple val(barcode), path(bam_file), val(sample_name)

    output:
    path "${sample_name}.bam"

    script:
    """
    cp ${bam_file} ${sample_name}.bam
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

process SingleReads {
    publishDir "${workflow.launchDir}/results/SingleReads/"

    input:
    tuple path(bam), path (bai)
    each path(reference_file)
    each path(single_region)

    output:
    path "*tsv"

    script:
    """
    modkit extract calls \
        --include-bed ${single_region} \
        --ref ${reference_file} \
        --filter-threshold 0.90 \
        ${bam} ${bam.baseName}.tsv
    """
}

process Rreport_all {
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

process Rreport_single {
    publishDir "${workflow.launchDir}/results", mode: 'copy'
    
    input:
    path sinread_files
    
    output:
    path "report.html"
    
    script:
    """
    Rscript ${workflow.projectDir}/scripts/generate_report.R \
        --input ${sinread_files} \
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
    
    // Reading sample sheet
    sample_map = Channel.fromPath(params.sample_sheet)
                        .splitCsv(sep: ',', header: true)
                        .map {row-> [row.barcode, row.sample_name]}

    // Extraction of barcode numbers
    bam_ch = demux_ch.flatten()
                    .map { file ->
                        def matcher = (file.name =~ /output\.bc_(\d+)\.bam/)
                        def barcode = matcher[0][1]
                        [barcode, file]
    }

    // Matching barcodes with samples
    renamed_ch = bam_ch.join(sample_map, by: 0)
                        .map { barcode, bam_file, sample_name -> tuple(barcode, bam_file, sample_name)
        }

    renamed_bams = change_name(renamed_ch)

    // Sort and index
    sorted_ch = sort_index_bam(renamed_bams.flatten())

    // SingleReads per bam
    sinread_ch = SingleReads(
        sorted_ch,
        reference_ch,
        channel.fromPath(params.region)
    )

    // Pileup per bam
    pileup_ch = pileup(
        sorted_ch,
        reference_ch
    )
    
    // Generate report
    Rreport_all(pileup_ch.collect())
    Rreport_single(sinread_ch.collect())
    
}