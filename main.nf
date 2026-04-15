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
    storeDir "${workflow.launchDir}/results/demultiplexed_bams/${modif}"

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

process merge_pileups {
    storeDir "${workflow.launchDir}/results/pileups"
    
    input:
    tuple val(sample_name), path(bed_files)
    
    output:
    tuple val(sample_name), path("${sample_name}_merged.bed")
    
    script:
    """
    cat ${bed_files} > ${sample_name}_merged.bed
    """
}

process SingleReads {
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
    errorStrategy 'ignore'
    
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

process merge_singlereads {
    storeDir "${workflow.launchDir}/results/SingleReads"
    errorStrategy 'ignore'

    input:
    tuple val(sample_name), path(tsv_files)
    
    output:
    tuple val(sample_name), path("${sample_name}_merged.tsv")
    
    script:
    def files = tsv_files instanceof List ? tsv_files : [tsv_files]
    def first  = files[0]
    def all    = files.join(' ')

    """
    # Header from first file
    head -1 ${first} > ${sample_name}_merged.tsv

    # Data (no header) from all files
    for f in ${all}; do
        tail -n +2 \$f >> ${sample_name}_merged.tsv
    done
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
    // If multiplexing, run seqtagger to demultiplex BAMs and rename them with sample names
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
        if (!line.startsWith('barcode')) {
            def parts = line.split(',')
            sample_map[parts[0]] = parts[1]
        }
    }

    // Parse BAMs and lookup sample names in sample_sheet
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
    sorted_bam_ch = sort_index_bam(renamed_bam).sorted_bam

} else { // If no multiplexing, simply sort and index the single BAM
    sorted_bam_ch = sort_index_bam(bam_ch).sorted_bam
}

// Run per modification pileup and single read extraction, then merge results for each sample

pileup_input = sorted_bam_ch.combine(reference_ch)
pileup_result = pileup(pileup_input)

merged_pileups = pileup_result
    .groupTuple(by: 0)  // Group by sample_name
    .map { sample_name, modifs, beds -> [sample_name, beds.flatten()] }

merge_pileup_result = merge_pileups(merged_pileups)

sinread_input = sorted_bam_ch.combine(reference_ch).combine(region_ch)
sinread_result = SingleReads(sinread_input)

sinread_labeled = label_negative_calls(sinread_result)

merged_input_ch = sinread_labeled
        | groupTuple(by: 0)
        | map { sample_name, files ->
            tuple(sample_name, files.flatten().findAll { it != null })
        }
        | filter { sample_name, files ->
            files.size() > 0
        }

    merge_singlereads(merged_input_ch)
}