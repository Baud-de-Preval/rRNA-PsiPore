#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */

params.pod5_dir = 'data/pod5/'
params.model = 'rna004_130bps_sup@v5.3.0'
params.modif = 'pseU_2OmeU,inosine_m6A_2OmeA,m5C_2OmeC,2OmeG'
params.reference = 'data/reference/Homo_sapiens.rRNA.fasta'
params.seqtagger_sif = null
params.sample_sheet = 'data/sample_sheet.csv'
params.region = 'data/reference/region.bed'
params.multiplex = true

/*
* Pipeline processes
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
    storeDir "${workflow.launchDir}/results/demultiplexed_bams/${modif}"

    input:
    tuple val(modif), path(pod5_dir), path(reference_file), val(model)

    output:
    tuple val(modif), path("${modif}.bam")

    script:
    """
    dorado basecaller \
        ${model} \
        ${pod5_dir} \
        --modified-bases ${modif} \
        --reference ${reference_file} \
        --emit-summary \
        --mm2-opts "-x map-ont -N 0 -k 13" \
        > ${modif}.bam
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
    apptainer run --nv --no-home \
        --bind ${workflow.launchDir}/data:/data \
        ${sif_file} \
        mRNA -k /opt/app/models/b96_RNA004 -r \
        -i /data/pod5/ \
        -o /data/pod5/demux    
    """
}

process seqtagger_demultiplex {
    storeDir "${workflow.launchDir}/results/demultiplexed_bams/${modif}"

    input:
    tuple val(modif), path(bam_file), path(index), path(sif_file)

    output:
    tuple val(modif), path("output.bc_*.bam")

    script:
    """
    apptainer run --nv --no-home\
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
    publishDir "${workflow.launchDir}/results/demultiplexed_bams/${modif}", mode: 'copy'

    input:
    tuple val(modif), val(sample_name), path(bam_file)
    
    output:
    tuple val(modif), path("${sample_name}_${modif}.bam")
    
    script:
    """
    mv "${bam_file}" "${sample_name}_${modif}.bam"
    """
}

process sort_index_bam {
    publishDir "${workflow.launchDir}/results/sorted_bams/${modif}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple val(modif), path(bam_file)
    
    output:
    tuple val(modif), path("${bam_file.baseName}.sorted.bam"), path("${bam_file.baseName}.sorted.bam.bai"), emit: sorted
    
    script:
    """
    samtools sort -o ${bam_file.baseName}.sorted.bam ${bam_file}
    samtools index ${bam_file.baseName}.sorted.bam
    """
}

process pileup {
    publishDir "${workflow.launchDir}/results/pileups/${modif}", mode: 'copy'

    input:
    tuple val(modif), path(bam), path(bai)
    path reference_file

    output:
    tuple val(modif), path("*.bed")

    script:
    """
    modkit pileup ${bam} ${bam.baseName}.bed \
        --ref ${reference_file} \
        -n 40000 \
        --filter-threshold 0.90
    """
}

process SingleReads {
    publishDir "${workflow.launchDir}/results/SingleReads/${modif}", mode: 'copy'

    input:
    tuple val(modif), path(bam), path(bai)
    path reference_file
    path single_region

    output:
    tuple val(modif), path("*.tsv")

    script:
    """
    modkit extract calls \
        --include-bed ${single_region} \
        --ref ${reference_file} \
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
    pod5_dir_ch  = Channel.fromPath(params.pod5_dir, type: 'dir')
    reference_ch = Channel.fromPath(params.reference)
    region_ch = Channel.fromPath(params.region)
    model_ch = dorado_model(params.model)

    modif_ch = Channel.of(params.modif.split(','))
        .flatten()
        .map { it.trim() }

    basecall_input = modif_ch
        .combine(pod5_dir_ch)
        .combine(reference_ch)
        .combine(model_ch)

    bam_ch = dorado_basecall(basecall_input)

if (params.multiplex) {
    sif_ch = params.seqtagger_sif ? 
        Channel.fromPath(params.seqtagger_sif) : 
        download_seqtagger_image()

    index_ch = seqtagger_index(sif_ch)
    
    demux_input = bam_ch
        .combine(index_ch)
        .combine(sif_ch)
    
    demux_ch = seqtagger_demultiplex(demux_input)

    // Read sample sheet into map
    def sample_map = [:]
    file(params.sample_sheet).eachLine { line ->
        if (!line.startsWith('barcode')) {  // skip header
            def parts = line.split(',')
            sample_map[parts[0]] = parts[1]
        }
    }

    // Parse BAMs and lookup sample names
    renamed_ch = demux_ch.flatMap { modif, files ->
        files.collect { file ->
            def matcher = (file.name =~ /output\.bc_(\d+)\.bam/)
            if (matcher) {
                def barcode = matcher[0][1]
                def sample_name = sample_map[barcode]
                if (sample_name) {
                    [modif, sample_name, file]
                }
            }
        }
    }

    renamed_bam = change_name(renamed_ch)
    sorted_bam = sort_index_bam(renamed_bam)

} else {
    sorted_bam = sort_index_bam(bam_ch)
}

pileup_result = pileup(sorted_bam.sorted, reference_ch)
sinread_result = SingleReads(sorted_bam.sorted, reference_ch, region_ch)
}