# hAMRoaster

## **H**armonized **AMR O**utput comp**A**ri**S**on **T**ool **ER**


hAMRoaster is an analysis pipeline that can compare the output of nine different tools for detecting AMR genes and provide metrics of their performance. For most tools, hAMRoaster requires preprocessing with hAMRonization. 

#### Installing hAMRoaster

To install hAMRoaster, run `conda install -c ewissel hamroaster` . (upload to bioconda in progress)

#### Test your install

`hAMRoaster -h` 
*(notice the caps)*

## hAMRoaster Commands

hAMRoaster has several required and optional commands. At a minimum, users **must** provide the following flags. 

* `--ham_out` : the full path to the output file from hAMRonization. If using the data published with this tool, the file is in `study_data/`

* `--name`: a unique identifier for this run. This name will be the name of the output directory created for the hAMRoaster run, and the name will be included in all the output file names. 

* `--AMR_key` : the full path to the file of known AMR phenotypes; the file is expected to be a tsv in the following format: taxa_name, antibiotic tested, result of antibiotic resting, and testing standard. This matches the output for the [NCBI BioSample AntiBioGram](https://www.ncbi.nlm.nih.gov/biosample/?term=antibiogram%5bfilter%5d). An example file is in `study_data/mock_2_key.csv`.

##### Optional Arguments

* `--abx_map` : While hAMRoaster comes with it's own classifications for antibiotics and their drug class, users can provide their own classification scheme. This expects the drug class to in the first column, and the antibiotics or drugs that fit into that drug class in the next two columns. hAMRoaster's default classification is located in `db_files/cleaned_drug_class_key.csv`. 

* `--fargene` : If users want to test the output of fARGene runs, they can point to the directory of the fARGene output. Because fARGene requires multiple runs for all the built in models, hAMRoaster expects users to point to a generic `fargene_output` directory, and for each run to be a subdirectory named after the AMR model analyzed for that run. For example, for a fARGene run with the model flag `class_a`, hAMRoaster expects their to be a directory named `class_a` within the fargene output directory specified with this parameter, and that the `class_a` subdirectory contains the output of that fARGene run. 

* `--shortbred` : If users want to test the output of a shortBRED run, they can use this flag to specifcy the full path to the output file from shortbred (include the file in the path). 

* `--shortbred_map` : This flag can point to a specific file that maps the shortbred IDs to their AMR gene names. hAMRoaster included the mapping file created with the shortbred publication in 2016. It is used by default and located in `db_files/ShortBRED_ABR_Metadata.tab`. 


Example usage:
```
hAMRoaster --ham_out amr-benchmarking/hAMRoaster/study_data/ham_sum.tsv  --name test_conda --AMR_key amr-benchmarking/hAMRoaster/study_data/mock_2_key.csv --db_files amr-benchmarking/hAMRoaster/

## the above is used to test the conda install with locally stored study data (available on this repo)
```

# Replicating hAMRoaster publication analysis

The following will walk through the code for the analysis done as an initial test of hAMRoaster and published with the hAMRoaster pipeline, from creating simulated metagenomic sequences to analyzing the data in hAMRoaster. 

## Create Simulated Mock Community

First, a simulated mock community was created. FASTA files were downloaded for eight taxa selected from NCBI's BioSample for their extensive phenotypic antibiotic resistance. NCBI's tool ART was used to simulate FASTQs from a HiSeq 2500 these FASTA files at three coverage levels - 5, 50, and 100x coverage.

```
conda install -c bioconda/label/cf201901 art 

art_illumina -ss HS25 -sam -i reference.fa -p -l 150 -f 5 -m 200 -s 10 -o paired_dat ## -f flag changed for 5,50,100x coverage
```

Once this was done for each of the eight mock community members, the simulated fastqs (`paired_dat_[1|2].fq`), the eight fastqs were combined to create one fastq file with even coverage at each coverage levels.

`cat 'ls -1 A*1* ... | sort -V ' > combo_1.fq  ## the ... represents each taxa being identified in alphabetical order for the cat`

Once fastqs were combined, contigs were assembled using metaSPAdes. 

```
conda create -y -n spades-test -c conda-forge -c bioconda spades python=3.7.7 

conda activate spades-test 

spades.py --meta  --pe1-1 combo_1.fq --pe1-2 combo_2.fq -o fasta/ 
```

With this, three mock communities were generated, again, for the three coverage levels (5x,50x, and 100x).

## Running Mock Community through AMR Detection Tools

This section will walk through how we set up a conda environment for each of the bioinformatic tools included in the hAMRoaster publication, in no particular order.

### [shortBRED](https://github.com/biobakery/shortbred)

shortBRED13 v0.9.3 uses a set of marker genes to search metagenomic data for protein families of interest. The bioBakery team published an AMR gene marker database built from 849 AR protein families derived from the ARDB20 v1.1 and independent curation alongside shortBRED, which is used in this study.

```
conda create –n shortbred 

conda activate shorbred

conda install -c biobakery shortbred 
  
conda install python=3.7 

## in a dir called shortbred database 

wget https://raw.githubusercontent.com/biobakery/shortbred/master/sample_markers/ShortBRED_ARDB_Evaluation_markers.faa ## get database  

shortbred_quantify.py --markers /full/path/to/shortbred_database/ShortBRED_ARDB_Evaluation_markers.faa --wgs /full/path/to/fastas/contigs.fasta --results combo_shortbred.tsv --tmp combo_quant --usearch /full/path/to/usearch11.0.667_i86linux32 
```

### [fARGene](https://github.com/fannyhb/fargene)

fARGene v.0.1 uses Hidden Markov Models to detect AMR genes from short metagenomic data or long read data. This is different from most other tools which compare the reads directly. fARGene has three pre-built models for detecting resistance to quinolone, tetracycline, and beta lactamases, which we test here.

```
conda activate python2 

conda install -c conda-forge -c bioconda fargene 

fargene -i fastq/combo* --store-peptides --hmm-model class_a -o fargene_out_fa/class_a --meta --force 
```

The following arguments were used with the `--hmm-model` flag to run all the pre-built fARGene models, with the output directory matching the name and spelling of the `--hmm-model` flag: `class_a`, `class_b_1_2`, `class_b_3`, `class_c`, `class_d_1`, `class_d_2`, qnr, `tet_efflux`, `tet_rpg`, and `tet_enzyme`.

### [AMRFinderPlus](https://github.com/ncbi/amr/wiki)

AMR Finder Plus v.3.9.3 uses BLASTX translated searches and hierarchical tree of gene families to detect AMR genes. The database is derived from the Pathogen Detection Reference Gene Catalog and was compiled as part of the National Database of Antibiotic Resistant Organisms (NDARO).

```
conda install -y -c bioconda ncbi-amrfinderplus 

amrfinder -n fasta/contigs.fasta --plus --threads 3 -o AMRFinderPlus_out/combo_out 
```

### [RGI](https://github.com/arpcard/rgi)

RGI v5.1.1 uses protein homology and SNP models to predict ‘resistomes’. It uses CARD’s protein homolog models as a database. RGI predicts open reading frames (ORFs) using Prodigal, detects homologs with DIAMOND, and matches to CARD’s database and model cut off values. 

```
conda create --prefix=/home/ewissel/conda/rgi -c bioconda rgi=5.1.1 -c conda-forge -c defaults  ## note I needed to use the prefix flag for this to work

conda activate /home/ewissel/conda/rgi 

rgi -i fasta/contigs.fasta -o rgi_out -t contig 
```

### [ResFinder](https://cge.cbs.dtu.dk/services/ResFinder/)

ResFinder v4.0 is available both as a web based application or the command line. We use ResFinder 4 for this data, which was specifically designed for detecting genotypic resistance in phenotypically resistant samples. ResFinder aligns reads directly to its own curated database without need for assembly. 

Please note that this was run while ResFinder was being updated to version 4.0 (and the docker version was catching up). I ran into issues with the docker version of ResFinder in May/June, and instead used the web based ResFinder results in the publication study; hyperlink to the web based ResFinder in the heading of this section. 

```
docker build -t resfinder . 

docker run --rm -it -v $(pwd)/db_resfinder/:/usr/src/db_resfinder -v $(pwd)/results/:/usr/src/results resfinder -v $(pwd)/dat/combo_1.fq:/usr/src/combo_1.fq -ifq /usr/src/combo_1.fq  /usr/src/combo_2.fq -acq -db_res /usr/src/db_resfinder -o /usr/src/results 
```

### [ABRicate](https://github.com/tseemann/abricate)

ABRIcate v.1.0.1 takes contig FASTA files as inputs and compared reads against a large database created by compiling several existing database, including NCBI AMRFinder Plus, CARD, ResFinder, ARG-ANNOT, MEGARES, EcOH, PlasmidFinder, VFDB, and Ecoli_VF. ABRIcate reports on acquired AMR genes and not mutational resistance.

```
conda activate python3.6 

conda install -c conda-forge -c bioconda -c defaults abricate 

conda install -c bioconda perl-path-tiny 

abricate --check 
 
## if you get path::tiny error, run 

conda update --all 

conda update abricate 

#### 

abricate --list 

abricate fasta/contigs.fasta 
\abricate assembly.fa.gz 
```

### [sraX](https://github.com/lgpdevtools/sraX#srax)

sraX v.1.5 is built as a one step tool; in a single command, sraX downloads a database and aligns contigs to this database with DIAMOND21. By default, sraX uses CARD, though other options can be specified. As we use default settings for all tools, only CARD is used in this study for sraX.

```
conda create –name sraX 

conda activate sraX 

sraX –v 

sraX -i sim_dat/fasta/ -o sraX_out 
```

### [deepARG](https://github.com/gaarangoa/deeparg2.0)

deepARG v.2.0 uses a supervised deep learning based approach for antibiotic resistance gene annotation of metagenomic sequences. It combines three databases—CARD, ARDB, and UNIPROT—and categorizes them into resistance categories. 

```
conda create -n deeparg_env python=2.7.18 

conda activate deeparg_env 

conda install -c bioconda diamond==0.9.24 

#pip install deeparg==1.0.2 to resolve error

deeparg download_data -o /path/to/local/directory/ 

deeparg predict -i fasta/contigs.fasta -o deeparg_out_combo_fasta --model SS –d full/path/to/downloaded_data/deeparg_data 
 ```
 
 ### [StarAMR](https://github.com/phac-nml/staramr#staramr)
 
starAMR v.0.7.2 uses blast+ to compare contigs against a combined database with data from ResFinder, PointFinder, and PlasmidFinder. 

```
conda install -c bioconda staramr==0.7.2 

cpan Moo 

staramr search fasta/contigs.fasta -o staramr_out 
```

With this, all nine bioinformatic tools were run on the same mock community at three coverage levels and combined using [hAMRonization](https://github.com/pha4ge/hAMRonization). For shortBRED and fARGene, which are not a part of hAMRonization at the time of analysis, options were added to hAMRoaster to examine its output (as decribed above in the section on hAMRoaster flags). 

## Combining the Outputs

hAMRonization is a strong tool which creates a unified format for examining the output of many different AMR analysis tools. Note that all the tools accomidated by the hAMRonization format are not included in this study as they did not meet all the eligibility criteria. For a further explanation, look at the table in the hAMRoaster publication which walks through 16 tools for AMR gene detection and their reasons for inclusion/exclusion.

At the time of analysis, hAMRonization was under development and documentation was not yet complete. As such, the database versions and analysis software versions, which are required flags for hAMRonization, reflect a minimally sufficient string for this and not the true version used in this analysis. This is not a recommended use of this tool or these flags, but we are including this detail for the sake of true replication.

`conda install hAMRonization `

To generate the hAMROnized format for each tool:

* ABRicate: `hamronize abricate  abricate_out/abricate_out.txt --reference_database_version db_v_1 --analysis_software_version tool_v_1 --format tsv` 

* starAMR: `hamronize staramr staramr_out/resfinder.tsv --format tsv --output hamr_out/staramr.tsv --reference_database_version db_v_1 --analysis_software_version tool_v_1 `  

* deepARG: `hamronize deeparg deeparg_out/deeparg_out_combo_fasta.mapping.ARG --reference_database_version db_v_1 --analysis_software_version tool_v_1 --format tsv --output ham_out/deeparg_out.tsv `

* sraX: `hamronize srax sraX_out/Results/Summary_files/sraX_detected_ARGs.tsv --reference_database_version db_v_1 --analysis_software_version tool_v_1 --format tsv --output ham_out/srax_out.tsv  --reference_database_id basic --input_file_name sraX_detected_ARGs.tsv `

* ResFinder (txt file from wb-based processing): `hamronize resfinder4 ResFinder_results_tab_5x.txt  --reference_database_version db_v_1 --analysis_software_version tool_v_1 --output ham_out/resfinder_out.tsv --input_file_name contigs.fasta `

* AMRFinderPlus: `hamronize  amrfinderplus AMRFinderPlus_out/combo_out --format tsv --output ham_out/amrfinderplus.tsv --analysis_software_version tool_v_1 --reference_database_version db_v_1 --input_file_name combo_contig_fasta_out `

* RGI: `hamronize rgi --input_file_name rgi_report --analysis_software_version rgi_v1 --reference_database_version card_v1 rgi_out/rgi_out `
 
 
 To combine the hamronized outputs into one table: 
 
 `hamronize summarize * -t tsv -o ham_sum.tsv`
 
 ## Analyzing result
 
This analysis was motivated by the lack of an open-source pipeline for comparing the results once compiled. hAMRoaster was built so that results can easily be compared using several metrics across tools.
