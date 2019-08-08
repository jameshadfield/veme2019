# DENGUE-2 nextstrain build for VEME 2019

Data from [CADDE project](https://twitter.com/CaddeProject) for the Nextstrain tutorial at [VEME 2019, Hong Kong](https://rega.kuleuven.be/cev/veme-workshop/2019).

Please see [this virological post](http://virological.org/t/genomic-monitoring-of-dengue-virus-serotype-2-in-brazil-2019/312) for more information.

## Installing nextstrain locally
```
curl http://data.nextstrain.org/nextstrain.yml --compressed -o nextstrain.yml
conda env create -f nextstrain.yml
conda activate nextstrain
npm install --global auspice
```

## Clone this repo
```
## using https
git clone https://github.com/jameshadfield/veme2019.git
cd veme2019
## alternative method
git clone git@github.com:jameshadfield/veme2019.git
cd veme2019
```


## Running the bioinformatics analyses using Augur
This pipeline contains a series of individual augur steps.
We've collected them all into a "Snakemake" file, which will run them all for you via:

```
snakemake clean
snakemake
```

Alternatively, you can run each step individually:

```bash
## parse the metadata into TSV formats for nextstrain
python scripts/parseMetadata.py --metadataIn data/Reference.csv data/NewSequences.csv --metadataOut results/metadata.tsv --latlongs results/latlongs_per_strain.tsv --sequencesIn data/Reference.fas data/NewSequences.fas --sequencesOut results/sequences.fasta
## align genomes using MAFFT
augur align --sequences results/sequences.fasta --reference-sequence config/KF955363.gb --output results/aligned.fasta --fill-gaps --remove-reference
## build a tree using IQ-TREE
augur tree --alignment results/aligned.fasta --output results/tree_raw.nwk
## refine the tree & infer temporal structure using treetime
##  NOTE: this may prune some tips which are outliers in the root-to-tip plot
##  NOTE: this will root the tree to maximise temporal structure
##  Treetime: Sagulenko et al., Virus Evolution, vol 4, https://academic.oup.com/ve/article/4/1/vex042/4794731
augur refine --tree results/tree_raw.nwk --alignment results/aligned.fasta --metadata results/metadata.tsv --output-tree results/tree.nwk --output-node-data results/branch_lengths.json --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 4 --root best
## Infer "ancestral" mutations across the tree
augur ancestral --tree results/tree.nwk --alignment results/aligned.fasta --output results/nt_muts.json --inference joint
## Translate sequences
augur translate --tree results/tree.nwk --ancestral-sequences results/nt_muts.json --reference-sequence config/KF955363.gb --output results/aa_muts.json   
## Use Discrete Trait Reconstruction for location, region, state & country
augur traits --tree results/tree.nwk --metadata results/metadata.tsv --output results/traits.json --columns location region state country --confidence
## Some more metadata parsing!
cat results/latlongs_per_strain.tsv config/latlongs.tsv > results/latlongs.tsv
## Export for visualisation in Auspice
augur export --tree results/tree.nwk --metadata results/metadata.tsv --node-data results/branch_lengths.json results/traits.json results/nt_muts.json results/aa_muts.json config/genome_annotation_file.json --lat-longs results/latlongs.tsv --auspice-config config/auspice_config.json --output-tree auspice/DENV2_tree.json --output-meta auspice/DENV2_meta.json --colors config/colors.tsv
```


## Viewing the results locally
```
auspice view --datasetDir auspice
```
and then open a browser at [localhost:4000](http://localhost:4000)


## Viewing the results using nextstrain community

See [the nextstrain docs](https://nextstrain.org/docs/contributing/community-builds) for more information about community builds.


Since the exported JSONs (the final step in the pipeline above) have been committed to this repo, they can be automatically accessed via [nextstrain.org/community/jameshadfield/veme2019](https://nextstrain.org/community/jameshadfield/veme2019)



