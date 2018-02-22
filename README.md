# QKbusco
A set of scripts that merge BUSCO orthologous genes for phylogenetic analysis.

## Scripts
<i>QKbusco_merge.py</i> parses several runs of BUSCO, identify orthologous groups of genes, selects open reading frames, and prepares files for PRANK alignment.

<i>QKbusco_phylogeny.py</i> parses PRANK alignments and concatenates aligned sequences. Uses *QKphylogeny_alignment.py* for alignment assessment.

## Use
First, a white space delimited text file must be generated that contains the following:
  1) Abbreviation of the species (or accession)
  2) FASTA file of transcriptome (or ORFs) sequence
  3) Full table from BUSCO analysis
  4) Path to the HMMer output folder from BUSCO
  5) Longest ORFs (coding sequence) from TransDecoder

The script is current designed to only take TransDecoder output from version 4.1.0.

```
Acs	Achnatherum_splendens_trinity_assembly_v3b.fa run_Achnatherum_splendens_busco/full_table_Achnatherum_splendens_busco.tsv run_Achnatherum_splendens_busco/hmmer_output Achnatherum_splendens_trinity_assembly_v3b.fa.transdecoder_dir/longest_orfs.cds
Aet	Aegilops_tauschii_trinity_assembly_v3b.fa run_Aegilops_tauschii_busco/full_table_Aegilops_tauschii_busco.tsv run_Aegilops_tauschii_busco/hmmer_output Aegilops_tauschii_trinity_assembly_v3b.fa.transdecoder_dir/longest_orfs.cds
Agc	Agropyron_cristatum_trinity_assembly_v3b.fa run_Agropyron_cristatum_busco/full_table_Agropyron_cristatum_busco.tsv run_Agropyron_cristatum_busco/hmmer_output Agropyron_cristatum_trinity_assembly_v3b.fa.transdecoder_dir/longest_orfs.cds
Agd	Agropyron_desertorum_trinity_assembly_v3b.fa run_Agropyron_desertorum_busco/full_table_Agropyron_desertorum_busco.tsv run_Agropyron_desertorum_busco/hmmer_output Agropyron_desertorum_trinity_assembly_v3b.fa.transdecoder_dir/longest_orfs.cds
...	...	...	...	...
```

To run *QKbusco_merge.py* you use the following command:

```bash
python QKbusco_merge.py -b embryophyta_3193_OrthoDB9_orthogroup_info.txt -m poales_master_file.txt -s complete -o test -t 40 -p 64
```

The `-b` option is to specify the orthogroup file to be used, which in this case is for plants (embryophyta_3193_OrthoDB9_orthogroup_info.txt). The option `-m` specifies the master file described above. The option `-s` can be either `complete`, `allcomplete`, or `allfragmented`, which will include only genes that are complete and single copy, all complete genes (even if duplicated), anor all genes, including duplicated and fragmented versions, respectively. The option `-o` specifies the folder to send alignments. The option `-t` specifies the minimum number of accessions/species required for inclusion of a gene for alignments. The option `-p` will determine how many shell scripts to generate for PRANK alignment (a simple serial parallelization approach).

After running *QKbusco_merge.py*, the outputs from PRANK are merged using *QKbusco_phylogeny.py*.

```bash
python QKbusco_phylogeny.py -d 0.4 -s superalignmeny.phy *.phy
```

*QKbusco_phylogeny.py* will call *QKphylogeny_alignment.py* (from the [QKphylogeny](https://github.com/matthewmoscou/QKphylogeny) suite of scripts) to screen alignments for the required depth in coverage.

## Improvements
Additional improvements will be added to *QKbusco_merge.py* to provide parameters of the quality of alignments and overall statistics on coverage, species represents, etc. In addition, the use of updated TransDecoder versions.
