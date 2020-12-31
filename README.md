# Homologizer

> "It's too late to homologize..." -Not One Republic

It's not too late! Using the setHomeologPhase function in RevBayes, gene alignments with pairs of sequences can be adjusted so that labels across genes are homologized. Once adjusted, the data from many genes can be analyzed together, for example using a concatenated supermatrix or summary species tree analysis.

## Dependencies

RevBayes

R (version 4.0 or greater)

## Step 0: Input files

Gene alignments that include samples that have two homeolog sequences. The sequences within each gene should be phased. For target capture data, see the alleles_workflow and haplonerate methods in this same Phyloscripts repository for tips on generating phased haplotypes within a gene sequence. The gene alignments can be in FASTA or NEXUS format.

### genelist.txt
Create a file with the FULL path to each sequence file, one per line.

### genecopymap.txt
Create a comma-delimited (CSV) file that will be used to map your gene copies (sequence identifiers) to subgenomes and samples in homologizer. The first row of the file should have a header:

* Column 1: Sample
* Column 2: Subgenome
* Columns 3-N: GeneNames

For each sample, indicate a sample name, a subgenome name, and gene sequence names for each gene. If you wish to create a "dummy tip" for a sample, the sequence name must include the string `BLANK`. 

An example of the genecopymap can be found in the GitHub respository.


### revbayes_template.rev

This is a revbayes script that will be used for the joint inference of phylogeny. Parameters for the MCMC, including move proposals and substitution models can be set in this file, prior to running the first scripts. At minimum, the location of the genelist.txt and genecopymap.txt files should be edited in this file.

An example of the revbayes template file can be found in the GitHub repository.

## Step 1: Making RevBayes scripts

Given a set of gene alignments, we will need to prepare RevBayes files that have the appropriate alleles to switch. The paired labels are taken from the label swap file, and is added to a basic template for RevBayes.

### `revscript_maker.R`

This script will use the genecopymap to generate the RevBayes code to initialize the phase for every sample and the homeolog phasing proposals for the MCMC. Run the script from the command line with:

```
RScript --vanilla revscript_maker.R genecopymap.csv
```

The script will generate two output files in the current directory:

* `InitialPhase.Rev` contains commands to initialize phase for each sequence at each gene.
* `PhaseMoves.Rev` contains all of the homeolog phase move proposals needed for the MCMC.

Both of these files need to be in the directory where you execute RevBayes.

## Step 2: Running RevBayes

Edit the `revbayes_template.txt` file to have the location of your `genelist.txt` file, the two output files from the `revscript_maker.R` script, and indicate a prefix for the output files in the first lines of this file. 

From the command line you can then run:

```
rb revbayes_template.txt
```

## Step 3: Examining Phase Posteriors

After RevBayes has finished, evaluate the MCMC convergence with standard metrics-- for example, you load the `.log` files into Tracer to examine stationarity and Effective Sample Size (ESS).

If the run has converged, the `summarize_posteriors_new.R` script will evaluate the run for both joint and marginal posterior probabilities of phase inference across loci. The script requires the prefix used in revscript_maker.R (i.e., the prefix part of the phase log files) and the location of the `genecopymap.csv`:

```
Rscript --vanilla summarize_posteriors_new.R prefix genecopymap.csv
```

The script will output a joint phase probability file for each multi-tip sample, and a marginal phase probability file for each sample and gene.

[ **ToDo**: add some description of what these files actually mean ]

## Step 4: Visualization

## Step 5: Generating phased sequences

Fixing the names in sequence alignments and/or trees to match the MAP phase with some threshold, so they can be fed into other software (e.g. ASTRAL). Scripts not currently written.