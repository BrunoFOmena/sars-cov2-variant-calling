// Nextclade dataset is cached outside the work dir; its version is written to the results.
process NEXTCLADE_DATASET {
    storeDir "${projectDir}/data/nextclade_dataset"

    output:
    path "sars-cov-2", emit: dataset

    script:
    """
    nextclade dataset get --name sars-cov-2 --output-dir sars-cov-2
    """
}

process PANGOLIN {
    publishDir "${params.outdir}/lineage", mode: 'copy'

    input:
    path consensus

    output:
    path "pangolin_report.csv"
    path "pangolin_versions.txt"

    script:
    """
    cat ${consensus} > all_consensus.fasta
    pangolin all_consensus.fasta --outfile pangolin_report.csv --threads ${task.cpus}
    pangolin --all-versions > pangolin_versions.txt 2>&1
    """
}

process NEXTCLADE {
    publishDir "${params.outdir}/lineage", mode: 'copy'

    input:
    path consensus
    path dataset

    output:
    path "nextclade_report.tsv"
    path "nextclade_dataset_version.json"

    script:
    """
    cat ${consensus} > all_consensus.fasta
    nextclade run --input-dataset ${dataset} --output-tsv nextclade_report.tsv all_consensus.fasta
    python3 -c "import json;d=json.load(open('${dataset}/pathogen.json'));print(json.dumps({'name':'sars-cov-2','version':d.get('version')}))" > nextclade_dataset_version.json
    """
}
