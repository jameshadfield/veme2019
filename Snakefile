rule all:
    input:
        auspice_tree = "auspice/veme2019_tree.json",
        auspice_meta = "auspice/veme2019_meta.json"

rule files:
    params:
        input_sequences = ["data/Reference.fas", "data/NewSequences.fas"],
        input_metadata = ["data/Reference.csv", "data/NewSequences.csv"],
        input_latlongs = "config/latlongs.tsv",
        reference = "config/KF955363.gb",
        auspice_config = "config/auspice_config.json",
        colors = "config/colors.tsv"

files = rules.files.params

rule parse:
    message: "Parsing provided metadata & sequences into augur formats"
    input:
        metadata = files.input_metadata,
        sequences = files.input_sequences
    output:
        sequences = "results/sequences.fasta",
        metadata = "results/metadata.tsv",
        latlongs = "results/latlongs_per_strain.tsv"
    shell:
        """
        python scripts/parseMetadata.py \
            --metadataIn {input.metadata} \
            --metadataOut {output.metadata} \
            --latlongs {output.latlongs} \
            --sequencesIn {input.sequences} \
            --sequencesOut {output.sequences}
        """

rule join_latlongs:
    message: "Combining lat long data from different sources"
    input:
        latlongs = [rules.parse.output.latlongs, files.input_latlongs]
    output:
        latlongs = "results/latlongs.tsv"
    shell:
        """
        cat {input.latlongs} > {output.latlongs}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.parse.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """

rule tree:
    message: "Building tree using IQ-TREE"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    message:
        """
        Refining tree to add names to internal nodes and inferring a timetree.
        NOTE: this step can drop samples which are extreme outliers in the root-to-tip analysis
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root="best",
        coalescent = "opt",
        clock_filter_iqd = 4,
        date_inference = "marginal"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd} \
            --root {params.root}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """


rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits.json",
    params:
        columns = "location region state country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        auspice_config = files.auspice_config,
        lat_longs = rules.join_latlongs.output.latlongs,
        colors = files.colors,
        annotation_file = "config/genome_annotation_file.json"
    output:
        auspice_tree = rules.all.input.auspice_tree,
        auspice_meta = rules.all.input.auspice_meta
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.annotation_file} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta} \
            --colors {input.colors}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
