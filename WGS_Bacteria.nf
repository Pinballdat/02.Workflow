#!/home/datnguyen/nextflow nextflow

// Parameater of reads 
params.reads = "./*_{1,2}.fq.gz"
params.outdir = "./Results"


// fromFilePairs == tuple val(sample_id), path(reads)

// Process for Quality Control of raw reads
process FASTQC {
    conda "/home/datnguyen/anaconda3/envs/fastqc/"
    cpus 60
    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/01.Raw_QC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    tuple val(sample_id), path("*.zip"), emit: zip

    script:
    """
    fastqc -f fastq --threads ${task.cpus} -q ${reads}
    """
}

process TRIMMOMATIC {
    conda "/home/datnguyen/anaconda3/envs/trimmomatic/"
    cpus 60
    tag "TRIMMOMATIC on $sample_id"
    publishDir "${params.outdir}/02.Trimmed_reads", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.paired.fq.gz"), emit: paired
    tuple val(sample_id), path("*.unpaired.fq.gz"), emit: unpaired
    tuple val(sample_id), path("*_summary.txt"), emit: summary
    tuple val(sample_id), path("*.log"), emit:log


    script:
    fq_1_paired = sample_id + "_1.paired.fq.gz"
    fq_1_unpaired = sample_id + "_1.unpaired.fq.gz"
    fq_2_paired = sample_id + "_2.paired.fq.gz"
    fq_2_unpaired = sample_id + "_2.unpaired.fq.gz"

    """
    trimmomatic PE -threads ${task.cpus} -trimlog ${sample_id}.log -summary ${sample_id}_summary.txt ${reads[0]} ${reads[1]} ${fq_1_paired} ${fq_1_unpaired} ${fq_2_paired} ${fq_2_unpaired} SLIDINGWINDOW:4:25 MINLEN:100
    """
}


process FASTQC_2 {
    conda "/home/datnguyen/anaconda3/envs/fastqc/"
    cpus 60
    tag "FASTQC on $sample_id after run TRIMMOMATIC"
    publishDir "${params.outdir}/03.Trimmed_QC", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("*")

    script:
    """
    fastqc -f fastq -q ${reads} -t 60
    """
}

process KRAKEN_READ {
    conda "/home/datnguyen/anaconda3/envs/kraken2/"
    cpus 60
    tag "KRAKEN2 for reads on $sample_id"
    publishDir "${params.outdir}/04.Kraken2_read", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_kraken2.txt")
    tuple val(sample_id), path("${sample_id}_kraken2_report.txt")
    script:
    """
    kraken2 --paired --db /data18tb/datnguyen/ONT_Gam/bacteria_db --threads ${task.cpus} --output ${sample_id}_kraken2.txt --report ${sample_id}_kraken2_report.txt ${reads[0]} ${reads[1]}
    """
}

process SPADES {
    conda "/home/datnguyen/anaconda3/envs/spades/"
    cpus 60
    tag "SPADES on $sample_id"
    publishDir "${params.outdir}/13.careful_Assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // path "${params.outdir}/${sample_id}_assembly/contigs.fasta", emit: assembled_contigs
    tuple val(sample_id), path("contigs.fasta"), emit: contigs
    tuple val(sample_id), path("assembly_graph.fastg"), emit: assembly_graph
    tuple val(sample_id), path("*")

    script:
    """
    spades.py --careful --pe1-1 ${reads[0]} --pe1-2 ${reads[1]} --threads $task.cpus -o ./
    """
}

process BANDAGE {
    conda "/home/datnguyen/anaconda3/envs/bandage"
    tag "BANDAGE on $sample_id"
    publishDir "${params.outdir}/05.1.Bandage", mode: 'copy'

    input:
    tuple val(sample_id), path(graph)

    output:
    // path "${params.outdir}/${sample_id}_assembly/contigs.fasta", emit: assembled_contigs
    tuple val(sample_id), path("bandage_info.tsv"), emit: info
    tuple val(sample_id), path("bandage_image.png"), emit: image

    script:
    """
    Bandage info ${graph} --tsv > bandage_info.tsv
    Bandage image ${graph} bandage_image.png
    """
}

process SEQKIT {
    conda "/home/datnguyen/anaconda3/envs/seqkit/"
    tag "Seqkit seq on $sample_id"
    publishDir "${params.outdir}/05.Assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    seqkit seq -m 500 ${contigs} > filter_contigs.fasta
    """
}

process PLASME {
    conda "/home/datnguyen/anaconda3/envs/plasme/"
    cpus 60
    tag "Plasme on $sample_id"
    publishDir "${params.outdir}/05.2.plasme", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    git clone https://github.com/HubertTang/PLASMe.git
    python PLASMe/PLASMe.py ${contigs} ${sample_id}_plasme -d /data18tb/datnguyen/PLASMe/PLASMe/DB -t ${task.cpus} -m balance
    """
}

process QUAST {
    conda "/home/datnguyen/anaconda3/envs/quast/"
    cpus 60
    tag "QUAST on $sample_id"
    publishDir "${params.outdir}/06.Assembly_QC", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    // tuple val(sample_id), path("${sample_id}"), emit: results
    tuple val(sample_id), path("*")

    script:
    """
    quast.py -o ${sample_id} --threads $task.cpus ${contigs}
    """
}

process BUSCO {
    conda "/home/datnguyen/anaconda3/envs/busco/"
    cpus 60
    tag "BUSCO in $sample_id"
    publishDir "${params.outdir}/07.Assembly_completeness", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    busco -i ${contigs} -m genome -c $task.cpus -o ${sample_id} --out_path ./ --auto-lineage-prok -f
    """
}

process KRAKEN2_CONTIG {
    conda "/home/datnguyen/anaconda3/envs/kraken2/"
    cpus 60
    tag "KRAKEN2 for contig on $sample_id"
    publishDir "${params.outdir}/08.Kraken2_contig", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_kraken2.txt")
    tuple val(sample_id), path("${sample_id}_kraken2_report.txt")
    script:
    """
    kraken2 --db /data18tb/datnguyen/ONT_Gam/bacteria_db --threads ${task.cpus} --output ${sample_id}_kraken2.txt --report ${sample_id}_kraken2_report.txt ${contigs}
    """
}

process PROKKA {
    conda "/home/datnguyen/anaconda3/envs/prokka/"
    cpus 60
    tag "PROKKA on $sample_id"
    publishDir "${params.outdir}/09.Prokka", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*.faa"), emit: FAA
    tuple val(sample_id), path("*")

    script:
    """
    prokka --kingdom Bacteria --prefix ${sample_id}_annotation --cpus ${task.cpus} --force --outdir ./ ${contigs}
    """
}

process ABRICATE {
    conda "/home/datnguyen/anaconda3/envs/abricate/"
    tag "ABRICATE on $sample_id"
    cpus 60
    publishDir "${params.outdir}/10.AMR_Virulence", mode:'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    abricate --threads ${task.cpus} --db argannot --csv ${contigs} > ${sample_id}_argannot.AMR.csv
    abricate --threads ${task.cpus} --db card --csv ${contigs} > ${sample_id}_card.AMR.csv
    abricate --threads ${task.cpus} --db megares --csv ${contigs} > ${sample_id}_megares.AMR.csv
    abricate --threads ${task.cpus} --db ncbi --csv ${contigs} > ${sample_id}_ncbi.AMR.csv
    abricate --threads ${task.cpus} --db resfinder --csv ${contigs} > ${sample_id}_resfinder.AMR.csv
    abricate --summary *.AMR.csv > AMR_summary.tab
    abricate --threads ${task.cpus} --db plasmidfinder --csv ${contigs} > ${sample_id}_plasmidfinder.csv
    abricate --threads ${task.cpus} --db vfdb --csv ${contigs} > ${sample_id}_vfdb.csv
    """
}

process MLST {
    conda "/home/datnguyen/anaconda3/envs/mlst"
    tag "MLST on $sample_id"
    publishDir "${params.outdir}/11.mlst", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.csv
    """

}


process EGGNOGMAPPER {
    conda "/home/datnguyen/anaconda3/envs/eggnog_mapper/"
    cpus 60
    tag "EGGNOG_MAPPER on $sample_id"
    publishDir "${params.outdir}/12.eggnog_mapper", mode: 'copy'

    input:
    tuple val(sample_id), path(faa)

    output:
    tuple val(sample_id), path("*.annotations"), emit:annotations
    tuple val(sample_id), path("*")

    script:
    """
    emapper.py -m diamond --itype proteins -i ${faa} -o ${sample_id}.emapper --cpu ${task.cpus} 
    """
}

process COGCLASSIFIER {
    conda "/home/datnguyen/anaconda3/envs/cogclassifier"
    cpus 60
    tag "COGCLASSIFIER on $sample_id"
    publishDir "${params.outdir}/13.cogclassifier", mode: "copy"

    input:
    tuple val(sample_id), path(faa)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    COGclassifier -i ${faa} -o ./ -t 60 
    """
}

process GOINPUT {
    tag "Analysis GO id on $sample_id from Eggnog_mapper"
    publishDir "${params.outdir}/14.GO", mode:"copy"

    input:
    tuple val(sample_id), path(eggnog_res)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    cat ${eggnog_res} | cut -f1,10 | grep -v -e "-" -e "##" > ${sample_id}_GO_id.txt
    """
}

process GOCLASSIFIER {
    tag "GO ANALYSIS on $sample_id"
    publishDir "${params.outdir}/14.GO", mode:"copy"

    input:
    tuple val(sample_id), path(go_id)

    output:
    tuple val(sample_id), path("*")

    script:
    """
    python3 /data18tb/datnguyen/WGS_MinhIBT_KC2_6/GO_analysis.py -i ${go_id} -o ${sample_id}_GO_analysis.csv
    python3 /data18tb/datnguyen/WGS_MinhIBT_KC2_6/GO_chart.py -i ${sample_id}_GO_analysis.csv -o ${sample_id}_GO_graph.lv1.png
    python3 /data18tb/datnguyen/WGS_MinhIBT_KC2_6/GO_chart.lv2.py -i ${sample_id}_GO_analysis.csv -o ${sample_id}_GO_graph.lv2.png
    """
    
}

workflow {
    Channel.fromFilePairs(params.reads, checkIfExists: true).set { raw_reads_channel }
    fastqc_channel = FASTQC(raw_reads_channel)
    trim_channel = TRIMMOMATIC(raw_reads_channel)

    fastqc_channel_2 = FASTQC_2(TRIMMOMATIC.out.paired)
    kraken_channel_1 = KRAKEN_READ(TRIMMOMATIC.out.paired)
    spades_channel = SPADES(TRIMMOMATIC.out.paired)

    bandage_channel = BANDAGE(SPADES.out.assembly_graph)
    seqkit_channel = SEQKIT(SPADES.out.contigs)
    plasme_channel = PLASME(seqkit_channel)

    quast_channel = QUAST(seqkit_channel)
    busco_channel = BUSCO(seqkit_channel)

    kraken_channel_2 = KRAKEN2_CONTIG(seqkit_channel)
    prokka_channel = PROKKA(seqkit_channel)

    abricate_channel = ABRICATE(seqkit_channel)
    mlst_channel = MLST(seqkit_channel)
    eggnog_channel = EGGNOGMAPPER(PROKKA.out.FAA)

    cog_channel = COGCLASSIFIER(PROKKA.out.FAA)
    go_input_channel = GOINPUT(EGGNOGMAPPER.out.annotations)
    go_classifier_channel = GOCLASSIFIER(go_input_channel)
}



/*
==========================================================================
Workflow Event Handler
*/

workflow.onComplete {

    println(workflow.success ? """
    Pipeline execution summary
    ------------------------------
    Completed at: ${workflow.complete}
    Duration: ${workflow.duration}
    Success: ${workflow.success}
    workDir: ${workflow.workDir}
    exit Status: ${workflow.exitStatus}
    """:"""
    Failed: ${workflow.errorReport}
    exit status: ${workflow.exitStatus}
    
    """
    )
}
