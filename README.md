# **OrthoPhyl2 (v2.2.1)**: Bigger and Better Orthology-based Phylogenomics
# Now Freaturing **ReLeaf!**: Rapidly Add Samples to Previous OrthoPhyl Runs
<br /> <br />
## Table of Contents

1. [About](#About)
2. [Getting Started](#GettingStarted)
    - [Singularity](#Singularity)
    - [Manual Install](#ManualInstall)
    - [Test Install](#TestInstall)
3. [Running OrthoPhyl](#RunningOrthoPhyl)
    - [Run Examples](#RunExamples)
4. [More Notes on OrthoPhyl](#NotesOnOrthoPhyl)
5. [Known Errors, Issues, and What-have-yous](#KnownErrors)
6. [Future Capabilities](#FutureCapabilities)
7. [OrthoPhyl Citation](#OrthoPhylCitation)

<a name="About"></a>

## About
### For the version used in *OrthoPhyl – Streamlining large scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales* see the OrthoPhyl_1.0 branch.
### Developed at Los Alamos National Labs (LANL - C22064)
#### Written by Earl Middlbrook with input from Robab Katani at Penn State and Jeanne Fair at LANL.
#### The software is available through a GLPv3 open source licence. 
#### Purpose
This software is designed to generate phylogenetic trees from bacterial genome assemblies. While many methods use whole genome alignments to generate informative sites to base tress on, OrthoPhyl annotates bacterial genes, identifies orthologous sequences, aligns related proteins to inform transcript alignments, then builds species trees with two methods. The first is a conventional gene concatenation and ML tree estemation method. The second attempts to reconcile gene trees with a unified species tree using quartets (ASTRAL). Both methods allow filtering of gene lists on number of species represented, length, and gappiness in order to tune signal-to-noise ratio for tree estimation. 
The main advantages of this software pipeline are three fold: 1) It extends the evolutionary distance input species can represent (over whole genome alignment and k-mer methods) while maintaining phylogenetic resolution, 2) this software is designed to be very user friendly, requiring just a single to estimate trees from a directory of assemblies. Additionally a Singularity image is now available to avoid dependancy hell and 3) this pipeline is amenable to estimating trees for 1000s of bacterial genome assemblies. To handle large numbers of genomes, the pipline calculates a diversity-representing subset of genomes to run OrthoFinder on, then expands the found OrthoGroups to all assemblies with an iterative HMM search strategy.
![screenshot](/img/OP2.0_workflow.png)
Grey boxes indicate processes. Orange, tan, and purple boxes represent user input, intermediate files, and species tree outputs, respectively. Purple arrows show iterative approaches. The workflow is divided into four main tasks: a) annotate assemblies, clean-up files, and remove identical CDSs. If more than “N” assemblies are being analyzed, b1) identify a subset of diversity-spanning assemblies, b2) pass them through OrthoFinder to generate orthogroups, and b3) expand the OrthoFinder-identified orthogroups to the full dataset of assemblies through iterative HMM searches. c) Align full orthogroup protein sets, generate and trim matching codon alignments, then filter orthogroups by taxon representation. Finally, d) estimate species tree topologies with concatenated codon alignment supermatrices along with a gene tree to species tree consensus method.

<a name="GettingStarted"></a>
## Getting Started
### Compatability
#### This workflow has been tested on CentOS8 and debian(mint) machines, but should be pretty portable to other linux systems. 

<a name="Singularity"></a>
### Singularity
#### If you have Singularity you can run a prebuilt container. Grab the container which requires ~1.7gb space (at the moment):
```
singularity_images=~/singularity_images/
mkdir ${singularity_images}
cd ${singularity_images}
singularity pull library://earlyevol/default/orthophyl
```
Congradulations! You can skip down to testing the "install"!

<a name="ManualInstall"></a>
### Manual Install
#### Clone the Orthophylo repo
Cloning takes a minute because of the large test files
```
Path_to_gits=~/gits
mkdir ${Path_to_gits}
cd ${Path_to_gits}/
git clone https://github.com/eamiddlebrook/OrthoPhyl.git
cd OrthoPhyl
```

#### Conda/mamba installable dependencies
If you need to install conda or mamba we recomend mamba, however the commands are exactly the same (exept running the install)
```
cd ~/Downloads/ # or wherever you want to put the installer
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh
```
...or
```
cd ~/Downloads/ # or wherever you want to put the installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
#### Follow on screen instructions
Set up auto initialize. If you don't, it is up to you to figure out how to make the script happy with your decision, you can also run:
```
mamba init bash
```
On some Linux distros, using ```~/.bashrc``` only for interactive shells is enforced (Linux Mint). In this case, mamba will not be initialized correctly within OrthoPhyl. To fix, just copy the conda/mamba initialization block from the bottom of your ```~/.bashrc``` file (below, replace "..." with what is in your file), and place them at the bottom of ```~/.bash_profile```. 
```
# >>> conda initialize >>>
...
# <<< conda initialize <<<
```


#### Create conda environment and install dependencies
It is highly recomended that you create an separate environment for this install. Might take some time....
```
mamba create -n orthophyl --file $Path_to_gits/OrthoPhyl/orthophyl_env.XXX.txt -c bioconda -c conda-forge
```
...or you can install different versions if neccessary. Compatability should be pretty good, but no promises.
```
mamba create -n orthophyl -c bioconda -c conda-forge \
git prodigal=2.6.3 orthofinder=2.5.4 \
bbmap=39.01 fasttree=2.1.11 hmmer=3.3.2 pal2nal=14.1 \
ete3 raxml=8.2.12 trimal=1.4.1 \
parallel=20160622 r-essentials=4.1
```
#### Sometimes R (r-essentials) can be a pain. 
If it is causing problems with the conda/mamba install, remove it and install manually:
Go to https://cran.r-project.org/mirrors.html and pick an appropriate mirror. Then chose the linux distro you are using, download package files, and follow install instructions.

#### To allow conda changes to take effect:
```
source ~/.bashrc
```
Alternatively, close and reopen terminal. 
#### To remove parallel's citation reminder (BUT DONT FORGET TO CITE!) run:
```
mamba activate orthophylo
parallel --citation
```
#### Test R install
```
Rscript --help
```

#### Other dependencies
This reflects how I like to organize my machine, pick what works for you. The control_file.paths reflects this setup. If you choose to install the below packages in different locations, just change control_file.paths to reflect this.
```
cd ${Path_to_gits}/
# Install ASTRAL
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL
unzip Astral.5.7.8.zip #change to curren version if needed

cd ${Path_to_gits}/
git clone https://github.com/nylander/catfasta2phyml.git
git clone https://github.com/dportik/Alignment_Assessment.git
cd Alignment_Assessment/
# convert script to python3
2to3 -w Alignment_Assessment_v2.py 

# install fastANI
cd ~/Downloads
wget https://github.com/ParBLiSS/FastANI/releases/download/v1.33/fastANI-Linux64-v1.33.zip
unzip fastANI-Linux64-v1.33.zip
# For macOS !!!!! NOT SOLVED !!!!!!
# git clone https://github.com/ParBLiSS/FastANI.git
# and follow install instructions. This has been a huge pain on my machine...not solved
mkdir ~/apps/
mv fastANI ~/apps/ # or anywhere else you would like to put it. Change control_file.required to reflect path
```

#### Edit ```${Path_to_gits}/OrthoPhyl/control_file.paths``` to reflect system specific locations and conda environment name if binaries and conda env are in different locations than discribed above

<a name="TestInstall"></a>
## Test Install 
### Test Singularity container
You must specify an output directory (-s) if running through singularity, because the script will try to write directly to the container (wont work)
```
singularity run ${singularity_images}/OrthoPhyl.0.9.3.sif -T TESTER_chloroplast -s ./tester_chloroplast_output -t 4
```
#### Test Manual install
To test 'conda activate' within OrthoPhyl, make sure the conda environment is not activated
```
mamba deactivate
```
#### Run super "fast" test locally
Tested on RHEL 8.5 machine with Intel Core i7-8700 CPU and 16gb ram (~8 minute runtime using 3 cores)
```
bash ${Path_to_gits}/OrthoPhyl/OrthoPhyl.sh -T TESTER_fasttest -t 3
```
There should be a directory created at OrthoPhyl/TESTER/Workflow_test.fasttest$(date +%m-%d-%Y)
If the test was successful, there should be 4 species trees found in the FINAL_SPECIES_TREES

#### run Chloroplast test locally
Tested on RHEL 8.5 machine with Intel Core i7-8700 CPU and 16gb ram (~8 minute runtime using 3 cores)
```
bash ${Path_to_gits}/OrthoPhyl/OrthoPhyl.sh -T TESTER_chloroplast -t 3
```
There should be a directory created in OrthoPhyl/TESTER/Workflow_test.chloroplast$(date +%m-%d-%Y)
If the test was successful, there should be 4 species trees found in the FINAL_SPECIES_TREES

#### run bigger test locally
This script takes about 20 hr to complete with 20 cores. Most of this is ML tree building. Will make an artificial set of truncated genomes later
```
bash OrthoPhyl.sh -T TESTER -t 3
```

<a name="RunningOrthoPhyl"></a>

## Running OrthoPhyl 
```
USAGE: OrthoPhyl.sh -g Path_to_directory_of_assemblies -s directory_to_store_output
# ALL arguments are optional if set with \"-c control_file.your_args\"
#   Many default parameters are set in control_file.defaults
#   I will work to expose the more useful ones in later versions of OP
Required:
-g  path to genomes directiory
or
-a  paths to protien and transcript directories.
       They should be delared as \"-a path_to_transcript_dir,path_to_prot_dir\"
-s  path to the main directory for output
Optional:
-t  threads to use [4]
-R|--rigor	Set the overall analysis rigor. This overrides most conflicting parameters set with other arguments except "-o\|--omics" and "-m\|--min_frac_orthos" for a fast run. (fast, medium, full) 
	full: Run Iqtree2 on CDS and PROT concatenated alignments with partion merging and GTR + freerate family model testing for CDS and standard PROT models. 
		Additionally, run ASTRAL gene to species tree consensus on CDS and PROT sequences
		Analyses will be performed on strict and relaxed single copy orthologs (found in >=30% taxa)
	medium: Run Iqtree2 with GTR and FreeRate models for CDS and standard models for Prots
	fast: Run FastTree on concatenated PROT alignments using only strict single copy orthologs 
-p  phylogenetic tree software to use astral, fasttree, raxml, and/or iqtree [\"fasttree iqtree astral\"]
	i.e. -p \"fasttree iqtree astral\"
-o  "omics" data to use for tree building ([CDS], PROT, BOTH)
	for divergent sequences, it is good to compare protein trees to 
	nucleotide trees to identify artifacts of saturation (long branch attraction)
-c	path to a control file with required variables and any optional ones to override defaults.
	Will override values set on command line! [NULL]
-x  trimal paramerter string (in double \"quotes\")
-r  flag to rerun orthofinder on the ANI_shorlist (true/[false])
-n  Max number of proteomes to run through OrthoFinder.
	If more than this many assemblies are provided, a subset of proteomes (based on genomes/transcripts ANI) will be chosen for OrthoFinder to chew on [20]
-m  Minimum fraction of total taxa per orthogroup to consider it for the relaxed SCO dataset.
        Expects a float from 0-1
        A value of 0 or 1 will lead to only estimating trees for the SCO_stict dataset.
        [0.30]
-d  Force ANI subsetting to run on transcript or genome Datasets. ([genome],transcript)
	Using \"-a\" implies \"-d transcript\".
	If \"-a\" is declared but you want to use the assemblies you have also provided, set \"-d genome\"
	If \"-a\" not used but you want to use transcripts (annotated within OrthoPhyl by Prodigal) for ANI subsetting, set \"-d transcript\"
-T  run test dataset, incompatable with -g|s|a (TESTER,TESTER_chloroplast)
-h  display a description and a super useful usage message
###############################################################\n
To run test datasets:
# Test of full bacterial genomes
bash OrthoPhyl.sh -T TESTER -t #threads
# Big test with ~100 orchid chloroplasts
bash OrthoPhyl.sh -T TESTER_chloroplast -t #threads
# reduced chloroplast dataset
bash OrthoPhyl.sh -T TESTER_fasttest -t #threads
# When running through Singularity an output directory is required:
singularity run \${singularity_images}/OrthoPhyl.XXX.sif -T TESTER -s output_dir -t #threads
```
<a name="RunExamples"></a>
### Run Examples 
#### Example1: Run OrthoPhyl on assemblies in ~/Projects/ASMS/ecoli/ using 12 cores within singularity and place all results and intermediate files in ~/Projects/phylogenetics/ecoli/
```
singularity run ${singularity_images}/OrthoPhyl.X.X.X.sif -g ~/Projects/ASMS/ecoli/ -s ~/Projects/phylogenetics/ecoli/ -t 12

```
#### Example2: Run the same OrthoPhyl command from the manual install
```
./OrthoPhyl.sh -g ~/Projects/ASMS/ecoli/ -s ~/Projects/phylogenetics/ecoli/ -t 12

```
#### Example3: Run the same OrthoPhyl command from the manual install, also incorporate preannoated samples with protein seqs in ```~/Projects/annots/protiens/``` and transripts in ```~/Projects/annots/transcripts/``` 
```
./OrthoPhyl.sh -g ~/Projects/ASMS/ecoli/ -a ~/Projects/annots/protiens/,~/Projects/annots/transcripts/ -s ~/Projects/phylogenetics/ecoli/ -t 12 

```
#### Example4: Run OrthoPhyl to build a species tree with 
+ ASTRAL (-p astral) 
+ from CDS (-o CDS) 
+ codon alignment gene trees build by iqtree used by ASTRAL (control_file.user contains just the line "export gene_tree_methods=("iqtree")") 
+ using only single copy orthologs found in 1/1 assemblies (-m 1)
+ using only 15 asssemblies to generate HMMs with OrthoFinder (-n 15)  
```
./OrthoPhyl.sh -g TESTER/genomes_chloroplast/ -s TESTER/TESTER_custom.chloro.B -m 1 -n 15 -c control_file.user -p astral -o CDS
```
#### Example5: Run OrthoPhyl to build a species tree with 
+ iqtree and ASTRAL (-p "iqtree astral") 
+ from protein alignments (-o PROT) 
+ protein alignment gene trees build by fasttree used by ASTRAL (default) 
+ using only single copy orthologs found in .3/1 assemblies (-m .3)
+ using 30 asssemblies to generate HMMs with OrthoFinder (-n 30)  
```
./OrthoPhyl.sh -g TESTER/genomes -s TESTER/TESTER_full_genomes -m .3 -n 30 -p "iqtree astral" -o PROT
```
<a name="NotesOnOrthoPhyl"></a>

## More Notes on OrthoPhyl: 
This software wraps many open source bioinformatic tools together with several custom programs. Genomes are annotated with Prodigal.  If a large number of genomes are to be analyzed, fastANI is used to estimate the pairwise genomic Average Nucleotide Identity (ANI) and then an OrthoPhyl function subsets the genomes to a minimal number which represent the diversity of the full set. Protein sequences from this subset (or full set for small numbers of genomes) are used as input for OrthoFinder to identify orthologous gene families. For the subset method, a HMM search, HMMER, is used to generalize the resulting orthogroup to the full data set. From there, orthogroup proteins are realigned with MAFFT, and these aignments are used to "codon" align the transcript sequences (generated by earlier annotation). These orthogroup transcript alignments are then trimmed (with TRIMAL) and filtered for strict single copy orthologs (SCO_strict) or SCOs found in at least X% of the input genomes (with X being tunable). Next, transcript alignmants are concatenated to generate super-maticies and used as input for species tree generation with either RAxML or fastTREE. Aditionally, gene trees are generated from the individual transcript alignments, which are used in gene tree to species tree estimation with ASTRAL.  

### Cleaning input genomes
The script OrthoPhyl/utils/gather_filter_asms.XX.sh was writen to streamiline aquiring all available assemblies for a specific taxon. It takes NCBI's taxID as and input (BrucellaTaxID=234). After the genomes are DL'd there are several simple genome filtering steps: length,N50,GC, etc. There are also a few more soficticated filtering methods: CheckM is used to assess completeness and contamination, while dedupe from BBmap is used to reduce duplicated contigs. 
It is absolutely imparative to clean up assemblies to get the best results from OrthoPhyl. For instance, when looking at the output for all Brucella accessions, checkM output shows many assemblies have duplicated "marker genes", if these are from falsely duplicated contigs in the asm, they will lead to removal of the ortholog group from both the strict SCO and relaxed SCO gene sets within the OrthoPhylo workflow. Furthermore, dedupe.sh (from bbmap) identifies many Duplicated or Contained contigs, some of them >100kb long. Again, if any of the duplicated contigs contain what would otherwise be SCOs, the SCOs will be removed. In essence, having duplicates poisons the analysis by severely reducing the number of Orthologs for downstream analysis. This problem is exacerbated by including 1) less curated assemblies (GenBank) and 2) the total number of assemblies fed to the analysis pipeline.

### What this software does not do: 
Generate trees that are ready for publication without parameter tuning or manual inspection. Reconstructing trees from whole genomes (gene tree reconciliation) requires many many steps, all of which have parameters that will differ based on the input sequences. Some importent outputs too look at: input genomes quality (checkM output), assembly subset used for Ortholog model generation (assembly shortlist), number of strict/relaxed single copy orthologs (drops quickly with additional assemblies), phylogenetic signal for transcript alignments, missing data in alignments (per gene and concatenated alignments)...to name a few. #### This pipeline also does not robunstly compute SCOs (like BUSCOs). It grabs all the SCOs it can from a dataset, but does not do any modeling to ensure species tree-like behavior of the individual genes and also does not deal with paralogs at all.

<a name="KnownErrors"></a>
## Known Errors:
#### If you have special characters in you fasta file names:
Filenames with special characters will likely make this workflow fail. As I encountered these characters, i will add fixes accordingly in Orthophyl.XX.sh, function SET_UP_DIR_STRUCTURE.
#### If the combination of fasta file name and contigs within are very long: 
Mafft will truncate the sequence names generated by prodigal. This will make the PAL2NAL module fail on many orthogroups. It manifests as a lot of cat commands failing during trimming because they cant find the transcript file PAL2NAL is suposed to spit out. This can be seen in the log file (in slurm_out). It's normally caused by very redundant sequence identifiers in the contig names. Can be fixed with something like this:
```
sed -i 's/REDUNDANT_STUFF/_/g' < SEQUENCE.fasta
```
#### If genomes are less than ~75% identity... 
FastANI doesnt calculate an ANI value for very divergent genomes. In this case my script ANI_picking.py assigns an ANI of 50 (hard coded) for these comparisons. This should not pose a problem if the number of clades with ~75% divergence internally are ~= ANI_shortlist number. The goal of the ANI stuff is to find representetives of the evolutionaly breadth of the dataset, so missing representatives from some clades will likelly not affect the final result (there is a lot of filtering later on in the workflow).
I will implement a more robust ANI estimator soooon... 
#### If you get an error like: 
``` 
/cm/local/apps/slurm/var/spool/job49381/slurm_script: line 511: J % percent: division by 0 (error token is "percent") 
```
It means that an upstream process failed and the directory used to enumerate a loop is empty

I will attempt to make errors easier to track...

#### I have gotten an error with gather_filter_genomes when pulling extra genomes via wget:
```
dyld: Symbol not found: _libiconv_open
```
This appeared to be a problem with wget, which was resolved by updating through conda/mamba.
#### Segmentation fault during Trimming/backtranslating
If trimming during backtranslating removes all sites you will get a segmentation fault:
```
gits/OrthoPhyl/script_lib/functions.sh: line 973: 3034614 Segmentation fault      (core dumped) trimal -in $i $trimal_parameter -ignorestopcodon -backtrans $SequencesCDS/${base}.CDS.fa -out $CDS_codon_alignment/${base}.codon_aln.fa >> $wd/logs/trimal.TRIMAL_backtrans 2>&1
```
This just results in there being a missing alignment. Since later steps are enumerated on what alignments are present, this doesnt cause a further error. Since the missing alignment is garbage, it will not be missed.
#### SyntaxWarning during OrthoFinder run
Due to a python version update, ETE3 throws and error about escape characters during the OrthoFinder run:
```
/panfs/biopan03/home/earlm/mambaforge/envs/orthophyl5/bin/scripts_of/tree.py:1422: SyntaxWarning: invalid escape sequence '\-'
  """
/panfs/biopan03/home/earlm/mambaforge/envs/orthophyl5/bin/scripts_of/newick.py:54: SyntaxWarning: invalid escape sequence '\['
  _ILEGAL_NEWICK_CHARS = ":;(),\[\]\t\n\r="
```
etc...
This DOES NOT affect the results. In future python versions this warning will upgrade to an error and kill ete3.  I will see if there is an update that 

<a name="FutureCapabilities"></a>
## Future Capabilities:
+ DONE: add assembly filtering, now very incomplete assemblies will be used, so the number of strict single copy orthologs could be drastically reduced (maybe to zero). Filtering done with gather_filter_asm
+ DONE: Quantify phylogenetic info per gene (https://github.com/dportik/Alignment_Assessment.git)
+ Check genome assembly file not empty at start of Main pipeline....will crash pipe at Orthofinder run. This happened because of a file transfer mishap for me :( 
+ DONE: Allow users to submit pre-annotated assemblies. This will open up OrthoPhyl to all taxa (after making sure all software is compatable with alternative codon tables)
+ Test evolutionary models of genes with ete3. Then test GO enrichment, AMR, Virulence genes
+ Look at tree wide paralog numbers. Do GO analysis...
+ Identify HGT from Transcript alignments. HGTs specific to any group?
+ allow "protected assemblies" when filtering. Perhaps your favorite assembly is crap, but you really want it in the tree. A couple of bad assemblies shouldnt reduce the number of SCOs that much
+ allow usage of precomputed Orthogroup HMMs, build new models from hits, then grab novel orthogroups from the remaining protiens
+ DONE!!!: Related to above: allow adding assemblies to pre-run pipeline. i.e. use precomputed hmms to identify orthologs, and add them to alignemnts and regenerate trees
+ DONE: allow users to pick which dataset to use for tree building (currently strict and relaxed are used). This would greatly reduce pipeline run time.
+ change how the list of OGs for HMM  searching is inumerated. Currently tries with OG that have be filtered out because of paralogs, then gets mad because there is no multifasta to align (doesnt change outcome, just messy)
+ Depending on the step, keyboard interupt just hangs. Need to figure this out. Should be simple in the loops. Why it hangs in IQTREE, I don't know. 


<a name="OrthoPhylCitation"></a>
## OrthoPhyl Citation:
Earl A Middlebrook, Robab Katani, Jeanne M Fair, OrthoPhyl—streamlining large-scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales, G3 Genes|Genomes|Genetics, Volume 14, Issue 8, August 2024, jkae119, https://doi.org/10.1093/g3journal/jkae119
## Citations for dependencies
### If there are nested dependencies that I am missing citations for, please let me know
### [Coming soon]
+ prodigal 
+ orthofinder 
+ bbmap
+ fastTree*
+ hmmer
+ mafft
+ prodigal 
+ ete3
+ IQTree*
+ raxml*
+ trimal
+ parallel
+ catfasta2phyml
+ ASTRAL
+ fastANI
+ R
+ Alignment_Assessment
+ blast
+ scipy
+ numpy
+ muscle
+ mmseq2
+ mcl
+ fastme
+ diamond

