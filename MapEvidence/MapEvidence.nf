nextflow.preview.dsl=2

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

params.genome = "/path/to/genome/assembly.fasta"
params.outdir = "results"

params.protein = "/path/to/protein/protein.fasta"

params.transcript = "/path/to/transcript/transcript.fasta"

log.info """
NBIS
  _   _ ____ _____  _____
 | \\ | |  _ \\_   _|/ ____|
 |  \\| | |_) || | | (___
 | . ` |  _ < | |  \\___ \\
 | |\\  | |_) || |_ ____) |
 |_| \\_|____/_____|_____/  Annotation Service

 Annotation preprocessing workflow
 ===================================

 General parameters
     genome          : ${params.genome}
     outdir          : ${params.outdir}
		 protein         : ${params.protein}
		 transcript      : ${params.transcript}

 """

workflow {

    main:
        Channel.fromPath(params.genome, checkIfExists: true)
            .ifEmpty { exit 1, "Cannot find genome matching ${params.genome}!\n" } \
// CHECK IF PROTEIN OR TRANSCRIPT else exit
        | map_evidence

}

workflow map_evidence {

    take:
        genome_assembly
				protein
				transcript

    main:
        map_transcript(genome_assembly)
				minimap2_to_gff(map_transcript.out)
        map_protein(protein)
}

process map_transcript {

    tag "${fasta_file.baseName} ; min length = ${params.min_length}"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path fasta_file
		path fasta_file

    output:
    path "${fasta_file.baseName}_min${params.min_length}.fasta"

    script:
    """
    seqtk seq -A $fasta_file -L ${params.min_length} > ${fasta_file.baseName}_min${params.min_length}.fasta
    """

}

process assembly_generate_stats {

    tag "${fasta_file.simpleName}"
    publishDir "${params.outdir}/stats", mode: 'copy'
    label 'GAAS'

    input:
    path fasta_file

    output:
    path "${fasta_file.baseName}_assembly_report"

    script:
    """
    gaas_fasta_statistics.pl --infile $fasta_file --output ${fasta_file.baseName}_assembly_report
    """
    // gaas_fasta_statistics.pl can be found in the NBIS GAAS repository
}

process busco {

    tag "$fasta"
    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
    path fasta
    each lineage

    output:
    path out

    script:
    out = "busco_${fasta.baseName}_${lineage}"
    """
    : "\${BUSCO_CONFIG_FILE:=/usr/local/config/config.ini}"
    export BUSCO_CONFIG_FILE
    busco -c ${task.cpus} -i $fasta -l $lineage -m genome --out $out
    """
}

workflow.onComplete {
    log.info ( workflow.success ? "\nAnnotation preprocessing complete!\n" : "Oops .. something went wrong\n" )
}
