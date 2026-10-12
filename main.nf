#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Pipeline parameters
 */

params.pod5_dir = 'data/pod5/'
params.model = 'rna004_130bps_sup@v5.2.0'
params.modif = 'pseU_2OmeU,inosine_m6A_2OmeA,2OmeG,m5C_2OmeC'
params.reference = 'data/reference/Homo_sapiens.rRNA.fasta'
params.seqtagger_sif = null
params.sample_sheet = 'data/sample_sheet.csv'
params.region = 'data/reference/region_Alban.bed'
params.multiplex = true
/*params.summary_file = 'data/pod5/sequencing_summary_PBK77343_f8441023_27ec0d94.txt'*/

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
    tuple val(modif), path("${modif}.bam"), emit: bam

    script:
    """
    dorado basecaller \
        ${model} \
        ${pod5_dir} \
        --modified-bases ${modif} \
        --reference ${reference_file} \
        --mm2-opts "-x map-ont -N 5 -k 13" \
        > ${modif}.bam
    """
}

process nanostat_sample {
    tag "${sample_name}_${modif}"
    storeDir "${workflow.launchDir}/results/QC/per_sample"
    errorStrategy 'ignore'

    input:
    tuple val(sample_name), val(modif), path(bam), path(bai)

    output:
    path("${sample_name}_${modif}_NanoStats.txt")

    script:
    """
    NanoStat --bam ${bam} \\
        --name ${sample_name}_${modif}_NanoStats.txt \\
        --threads ${task.cpus}
    """
}

process multiqc_report {
    publishDir "${workflow.launchDir}/results/QC", mode: 'copy'

    input:
    path nanostat_files

    output:
    path("combined_QC_report.html")
    path("combined_QC_report_data")

    script:
    """
    multiqc . --force --filename combined_QC_report.html
    """
}

process download_seqtagger_image {
    storeDir "${workflow.launchDir}/data"
    
    output:
    path 'seqtagger.sif'
    
    script:
    """
    apptainer pull seqtagger.sif docker://lpryszcz/seqtagger:1.0d
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
    export MPLCONFIGDIR=/tmp/matplotlib-\$\$
    apptainer run --nv --no-home \
        --bind ${workflow.launchDir}/data:/data \
        ${sif_file} \
        mRNA -k /opt/app/models/b96_RNA004 -r \
        -i /data/pod5/ \
        -o /data/pod5/demux    
    """
}

process bam_tag_merging {
    storeDir "${workflow.launchDir}/results/merged_bams"

    input:
    tuple val(label), path(bams)

    output:
    tuple val(label), path("${label}.merged.bam")

    script:
    """
    sorted=""
    for b in ${bams}; do
        samtools sort -N -@ ${task.cpus} -o nsort.\$b \$b
        sorted="\$sorted nsort.\$b"
    done
    ${projectDir}/script/merge_tag.py -t ${task.cpus} -o ${label}.merged.bam \$sorted
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
    publishDir { "${workflow.launchDir}/results/demultiplexed_bams/${modif}" }, mode: 'copy'

    input:
    tuple val(modif), val(sample_name), path(bam_file)
    
    output:
    tuple val(modif), val(sample_name), path("${sample_name}_${modif}.bam")
    
    script:
    """
    mv "${bam_file}" "${sample_name}_${modif}.bam"
    """
}

process sort_index_bam {
    storeDir "${workflow.launchDir}/results/sorted_bams/${modif}"

    input:
    tuple val(modif), val(sample_name), path(bam_file)
    
    output:
    tuple val(sample_name), val(modif), path("${bam_file.baseName}.sorted.bam"), path("${bam_file.baseName}.sorted.bam.bai"), emit: sorted_bam
    
    script:
    """
    samtools sort -o ${bam_file.baseName}.sorted.bam ${bam_file}
    samtools index ${bam_file.baseName}.sorted.bam
    """
}

process pileup {
    scratch true
    storeDir "${workflow.launchDir}/results/pileups"
    input:
    tuple val(sample_name), val(modif), path(bam), path(bai), path(reference_file)
    
    output:
    tuple val(sample_name), val(modif), path("${sample_name}_${modif}.bed")
    
    script:
    """
    modkit pileup ${bam} ${sample_name}_${modif}.bed \
        --ref ${reference_file} \
        -n 40000 \
        --filter-threshold 0.90
    """
}

process SingleReads {
    scratch true
    storeDir "${workflow.launchDir}/results/SingleReads/"
    input:
    tuple val(sample_name), val(modif), path(bam), path(bai), path(reference_file), path(region_file)
    
    output:
    tuple val(sample_name), val(modif), path("${sample_name}_${modif}.tsv")
    
    script:
    """
    modkit extract calls \
        --include-bed ${region_file} \
        --ref ${reference_file} \
        ${bam} ${sample_name}_${modif}.tsv
    """
}

/*
To compare the fraction of pileup and single read calls, 
we want to label the negative calls in the single read output with the modification name.
*/

process label_negative_calls {
    storeDir "${workflow.launchDir}/results/SingleReads/labeled"
    
    input:
    tuple val(sample_name), val(modif), path(tsv)

    output:
    tuple val(sample_name), path("${sample_name}_${modif}_labeled.tsv")

    script:
    """
    awk 'NR==1 {print; next} \$14 == "-" {OFS="\\t"; \$14 = "-(${modif})"; print; next} {print}' ${tsv} \
        > ${sample_name}_${modif}_labeled.tsv
    """
}

// ===== WORKFLOW =====

workflow {
    pod5_dir_ch  = Channel.fromPath(params.pod5_dir, type: 'dir')
    reference_ch = Channel.fromPath(params.reference)
    region_ch    = Channel.fromPath(params.region)
    model_ch     = dorado_model(params.model)

    // e.g. "pseU,m6A_DRACH" -> label modfications as "pseU_2OmeG" for merged BAM and downstream outputs
    def mods_label = params.modif.split(',').collect { it.trim() }.sort(false).join('_')

    modif_ch = Channel.of(params.modif.split(','))
        .flatten()
        .map { it.trim() }

    basecall_input = modif_ch
        .combine(pod5_dir_ch)
        .combine(reference_ch)
        .combine(model_ch)

    bam_ch = dorado_basecall(basecall_input)

    // One BAM per run carrying all modification tags
    merge_input = bam_ch.bam
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { items -> tuple(mods_label, items.collect { it[1] }) }

    merged_bam = bam_tag_merging(merge_input)

    if (params.multiplex) {
        sif_ch = params.seqtagger_sif ?
            Channel.fromPath(params.seqtagger_sif) :
            download_seqtagger_image()

        index_ch = seqtagger_index(sif_ch)

        demux_input = merged_bam
            .combine(index_ch)
            .combine(sif_ch)

        demux_ch = seqtagger_demultiplex(demux_input)

        def sample_map = [:]
        file(params.sample_sheet).eachLine { line ->
            if (!line.startsWith('barcode')) {
                def parts = line.split(',')
                sample_map[parts[0]] = parts[1]
            }
        }

        renamed_ch = demux_ch.flatMap { modif, files ->
            files.collect { file ->
                def matcher = (file.name =~ /output\.bc_(\d+)\.bam/)
                if (matcher) {
                    def sample_name = sample_map[matcher[0][1]]
                    if (sample_name) {
                        return [modif, sample_name, file]
                    }
                }
                return null
            }
        }.filter { it != null }

        renamed_bam   = change_name(renamed_ch)
        sorted_bam_ch = sort_index_bam(renamed_bam).sorted_bam

    } else {
        // sort_index_bam expects (modif, sample_name, bam)
        single_ch     = merged_bam.map { label, bam -> tuple(label, 'sample', bam) }
        sorted_bam_ch = sort_index_bam(single_ch).sorted_bam
    }

    nanostat_sample(sorted_bam_ch)
    multiqc_report(nanostat_sample.out.collect())

    // Multi-mod BAM -> one pileup / one extract per sample
    pileup(sorted_bam_ch.combine(reference_ch))

    sinread_result = SingleReads(sorted_bam_ch.combine(reference_ch).combine(region_ch))
    label_negative_calls(sinread_result)
}