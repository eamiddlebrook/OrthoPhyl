#!/bin/bash

USAGE () {
echo "
USAGE: OrthoPhyl.sh -g Path_to_directory_of_assemblies -s directory_to_store_output
# ALL arguments are optional if set with \"-c control_file.your_args\"
#   Many default parameters are set in control_file.defaults
#   I will work to expose the more useful ones in later versions of OP
Required:
-g|--genome_dir  path to genomes directiory
or
-a|--annotations  paths to protien and transcript directories.
       They should be delared as \"-a path_to_transcript_dir,path_to_prot_dir\"
-s|--storage_dir  path to the main directory for output
Optional:
-t|--threads  threads to use [4]
-p|--phylo_tool  phylogenetic tree software to use astral, fasttree, raxml, and/or iqtree [\"fasttree iqtree astral\"]
	i.e. -p \"fasttree iqtree astral\"
-o|--omics  "omics" data to use for tree building ([CDS], PROT, BOTH)
	for divergent sequences, it is good to compare protein trees to 
	nucleotide trees to identify artifacts of saturation (long branch attraction)
-c|--control_file	path to a control file with required variables and any optional ones to override defaults.
	Will override values set on command line! [NULL]
-x|-trimal_param  trimal paramerter string (in double \"quotes\")
-r|--rerun  flag to rerun orthofinder on the ANI_shorlist (true/[false])
-n|--num_OF_annots  Max number of proteomes to run through OrthoFinder.
	If more than this many assemblies are provided, a subset of proteomes (based on genomes/transcripts ANI) will be chosen for OrthoFinder to chew on [20]
-m|--min_num_orthos  Minimum fraction of total taxa per orthogroup to consider it for the relaxed SCO dataset.
        Expects a float from 0-1
        A value of 0 or 1 will lead to only estimating trees for the SCO_stict dataset.
        [0.30]
-d|--ani_Data  Force ANI subsetting to run on transcript or genome Datasets. ([genome],CDS)
	Using \"-a\" implies \"-d CDS\".
	If \"-a\" is declared but you want to use the assemblies you have also provided, set \"-d genome\"
	If \"-a\" not used but you want to use CDS (annotated within OrthoPhyl by Prodigal) for ANI subsetting, set \"-d CDS\"
-T|--test  run test dataset, incompatable with -g|s|a (TESTER,TESTER_chloroplast,TESTER_fasttest)
-h|--help  display a description and a super useful usage message
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
"
}


######################
# set up environment #
######################
echo "Loading required files..."

# grab the Orthophyl directory no matter where script was run from
#used for loading stuff from git repo
export script_home=$(dirname "$(readlink -f "$0")")

#############################################
######## Import funtions and stuff ##########
#############################################
# import main_script functions
source $script_home/script_lib/functions.sh

#load custom aliases...might get rid of
shopt -s expand_aliases
source $script_home/script_lib/bash_utils_and_aliases.sh


############################################
####  import control file variables  #######
############################################

# if not in a singularity or docker container, grab the "conda init" section of .bash_profile or .bashrc 
#    and import paths to external programs and the conda environment
if [[ -z ${SINGULARITY_CONTAINER+x} ]] && [[ -z ${DOCKER+x} ]]
then
    source $HOME/.bash_profile
    source $HOME/.bashrc
    source $script_home/control_file.paths || { echo "cannot open control_file.paths. This file is required for telling OrthoPhyl where to look for manually installed programs"; exit 1 ;}
fi


# Used to Catch all set variables later
tmpfile=$(mktemp)
declare -p > "$tmpfile"


# import pipeline defaults (will be overridden by command line args or -c control_file)
source $script_home/control_file.defaults || { echo "cannot open control_file.refaults. This file is required for telling OrthoPhyl what parameters to use if not provided on the CMD line"; exit 1 ;}


##############################################
##### if -T "TESTER" is set, run a small #####
#######  pipeline to test...pipeline  ########
##############################################
tester () {
if [[ $TESTER == "TESTER" ]]
then
	echo "
	################################################
	#####  Testing Workflow with Control Files  ####
	####  and genomes from the TESTER directory  ###
	################################################
	"
	source $script_home/TESTER/control_file.user
	if [ -d $store ]
	then
		rm -r $store
	fi
elif [[ $TESTER = "TESTER_chloroplast" ]]
then
	echo "
	################################################
	#####  Testing Workflow with Control Files  ####
	####  and genomes from the TESTER directory  ###
	##########   Chloroplast Edition!!!!   #########
	################################################
    "
	# NEED TO ADD COMPRESSED GENOME FILES TO TEST NEW FUNC
	source $script_home/TESTER/control_file.user_chloroplast
	if [ -d $store ]
    then
    	rm -r $store
    fi
elif [[ $TESTER = "TESTER_fasttest" ]]
then
	echo "
	################################################
	#####  Testing Workflow with Control Files  ####
	####  and genomes from the TESTER directory  ###
	##########    fasttest Edition!!!!     #########
	##########    (well not that fast)     #########
	################################################
	"
	# NEED TO ADD COMPRESSED GENOME FILES TO TEST NEW FUNC
	source $script_home/TESTER/control_file.user_fasttest
	if [ -d $store ]
    then
        rm -r $store
    fi
else 
	BAIL "PANIC: TESTER was set, but didnt equal 'TESTER', 'TESTER_chloroplast' or TESTER_fasttest"
fi
}


###################################################
#######  Handle command line Arguments ############
###################################################
#Setting Variables manually to override ones in control file
ARGS_SET=""

N=1 
L=${#1} 
while [[ $# -gt 0 ]]; do
	case ${1} in 
		-h|--help) if [[ ! -n ${2} ]] ; then 
				echo -e "\nHere is your help message!\n" 
				USAGE
				exit 0
			fi 
			# run function to print usage and exit
			echo -e "\nLooks like you used the '-h' with an argument...thats not right"
			echo "It's ok, here is your help message anyway" 
			USAGE
			exit 1
			;; 

		-g|--genome_dir) if [[ ! -n ${2} ]]; then 
				echo "looks like the argument for -g is messed up" 
				echo $USAGE
				exit 1
			fi
			input_genomes="$( cd "$(relative_absolute ${2})" && pwd )"
			genomes_provided=true
			ANI_genome=true
			echo $input_genomes
			if [[ ! -d "${input_genomes}" ]]; then
				echo "WARNING: -g declaring an input genome directory that does not exist...maybe check on that"
				exit 1
			fi
			ARGS_SET+=g
			shift
			shift ;;

		-s|--storage_dir) if [[ ! -n ${2} ]] ; then 
				USAGE 
				exit 1 
			fi
			store=$(relative_absolute ${2})
			echo $store
			ARGS_SET+=s
			shift
			shift ;;
		
		-o|--omics) if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			if [[ ${2} == "CDS" || ${2} == "PROT" || ${2} == "BOTH" ]] ; then
				TREE_DATA=${2}
				echo "Will build trees with $TREE_DATA (CDS,PROT, or BOTH)"
				ARGS_SET+=o
			else
				USAGE
				echo "WARNING: -o flag was used with a value other than CDS,PROT, or BOTH. You put '-o ${2}', should be e.g. '-o BOTH'"
				exit 1
			fi
			shift
			shift ;;
		
		-p|--phylo_tool) if [[ ! -n ${2} ]] ; then 
				USAGE 
				exit 1 
			fi
			tree_method=(${2})
			echo "#####"
			echo ${2}
			echo "#####"
			ARGS_SET+=p
			shift
			shift ;;

		-a|--annotations) if [[  ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			annots_provided=true
			ANI_trans=true
			trans_tmp=$(awk  -F ',' '{print $1}' <<< ${2})
			prots_tmp=$(awk  -F ',' '{print $2}' <<< ${2})
			# deal with relative or absolute paths correctly
			eval prots_tmp=$prots_tmp
			eval trans_tmp=$trans_tmp
			input_prots=$(relative_absolute ${prots_tmp})
			input_trans=$(relative_absolute ${trans_tmp})
			echo "Building codon based tree for protein in "$input_prots" and CDSs in "$input_trans
			# running some simple checks...
			# test if the CDS file is nucleotide and Prot is Prot
			CDS_chars=$(cat $input_trans/* | \
				head -n 20 | \
				grep -v ">" | \
				sed 's/\(.\)/\1\n/g' | \
				sort | \
				uniq | \
				wc -l ) 
			PROT_chars=$(cat $input_prots/* | \
				head -n 20 | \
				grep -v ">" | \
				sed 's/\(.\)/\1\n/g' | \
				sort | \
				uniq | \
				wc -l ) 
			if [ $CDS_chars -eq 22 ] ; then
				echo -e "Check your input CDS annotation file.\n Is it actually the CDS file and not the protein file?\n They should be delared as \"-a path_to_transcript_dir,path_to_protien_dir\""
				exit 1
			elif [ $CDS_chars -lt 4 ] ; then
				echo -e "Check your input CDS annotation file.\n Is the format a correct fasta file?"
				exit 1
			fi 
			if [ $PROT_chars -eq 4 ] ; then
				echo -e "Check your input PROT annotation file.\n Is it actually the PROT file and not the CDS file?\n They should be delared as \"-a path_to_transcript_dir,path_to_protien_dir\""

				exit 1
			elif [ $PROT_chars -gt 23 ] ; then
				echo -e "Check your input PROT annotation file.\n Is the format a correct fasta file?"
				exit 1
			fi 

			#test if number of CDS and PROT sequences are the same
			prot_num=$(ls $input_prots | wc -l)
			nucl_num=$(ls $input_trans | wc -l)
			if [[ ! $prot_num == $nucl_num ]]
			then
				echo "Number of provided Protein files and nucleotide (CDS) files do not match"
				echo "Perhaps there are extra files in the directories?"
				exit 1
			elif [[ ! -d "${input_prots}" ]] || [[ ! -d "${input_trans}" ]] ; then
				echo "WARNING: -a declared input annoation directories that do not exist...maybe check on that"
				echo " They should be delared as \"-a path_to_transcript_dir,path_to_protien_dir\""
				exit 1
			fi
			# add a test for nucleotide_file,AA_file 
			ARGS_SET+=a
			shift
			shift ;;
		-m|--min_num_orthos) if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			if (( $(echo ${2} | awk '{if ($1 > 1 || $1 < 0 ) print 1;}') ))
			then
				echo "-m set to $2, but should be between 0-1"
				USAGE
				exit 1
			fi
			if (( $(echo ${2} | awk '{if ($1 == 1 || $1 == 0 ) print 1;}') ))
			then
				relaxed=false
			else
				min_num_orthos=${2}
				ARGS_SET+=m
			fi
			shift
			shift ;;
		
		-t|--threads) if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			threads="${2}"
			ARGS_SET+=t
			shift
			shift ;;

		-c|--control_file) if [[ ! -n ${2} ]] ; then 
				USAGE
				exit 1
			fi
			# tests if control file given is relative for an absolute path
			control_file=$(relative_absolute ${2})
			echo $control_file
			ARGS_SET+=c
			# test if control_file is a file
			if [[ ! -f $control_file ]]; then
				echo "given control file (-c ${2}) does not exist"
				exit 1
			fi
			shift
			shift ;;

		-x|--trimal_params) if [[ ! -n ${2} ]] ; then 
				USAGE
				exit 1
			fi
			trimal_parameter=${2} 
			ARGS_SET+=x
			shift
			shift ;;

		-r|--rerun)	if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			orthofinder_rerun=${2}
			ARGS_SET+=r
			shift
			shift ;;

		-d|--ani_Data) 	if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			fi
			if [ ${2} == "genome" ] ; then
				ANI_genome=true
			elif [ ${2} == "CDS" ] ; then
				ANI_trans=true
			else
				echo "Looks like your arg for -d is invalid, check usage\n"
				USAGE
				exit 1
			fi
			ARGS_SET+=d
			shift
			shift ;;

		-n|--num_OF_annots) if [[ ! -n ${2} ]] ; then
				USAGE
				exit 1
			elif [[ ${2} -lt 4 ]] ; then
				echo "User asked for less than 4 input proteomes to be put through OrthoFinder." 
				echo "Please set \"-n 4\" or greater"
			fi
			ANI_shortlist=${2} 
			ARGS_SET+=n
			shift
			shift ;;

		-T|--test) if [[ ! -n ${2} ]] ; then
			USAGE
			exit 1
			fi
			export TESTER=${2}
			ARGS_SET+=T
			# set variables specific for running OrthoPhyl on a test data set.
			tester
			shift 
			shift;;

		\?) # Invalid option
			echo "Error: Invalid option"
			USAGE
			exit;;
		*) USAGE
		exit 1 ;;
	esac
done
if [[ -n ${1} ]] ; then
	echo "!!!!: you have a hanging argument without a flag!"
	USAGE
	exit 1
fi

####################################################################################
# test for incompatable args (if) and makes sure required args are present (elif) ##
####################################################################################
if [[ "$ARGS_SET" == *@(g|c|a)*@(T)* ]] || [[ "$ARGS_SET" == *@(T)*@(g|c|a)* ]]
then
	echo "!!!!: -T was set along with -g/c/a, which are incompatable args"
	USAGE
	exit 1
elif [[ ! "$ARGS_SET" == *c* ]] && [[ ! "$ARGS_SET" == *@(g|s)*@(g|s)* ]] && [[ ! "$ARGS_SET" == *T* ]] && [[ ! "$ARGS_SET" == *@(a|s)*@(a|s)* ]]
then
	echo -e "WARNING: Required arguments were not given, you need to provied either \n\ta control file with \"-c control_file\"\n\t-g genome_directory -s storage_directory \n\t-n prot_dir,trans_dir -s storage_directory \n\t-g genome_directory -a prot_dir,trans_dir -s storage_directory  \n\tor \"-T TESTER\" to start a test run"
	echo "Arguments that you set are -"${ARGS_SET}
	USAGE
	exit 1
elif [[ "$ARGS_SET" == *c* ]]
then
	echo "Args comming in through file set with -c, other args seems good"
fi

#############################################
# test if "-p tree_method" is set correctly #
#############################################
if [[ " ${tree_method[*]} " =~ " raxml " ]] || [[ " ${tree_method[*]} " =~ " fasttree " ]] || [[ " ${tree_method[*]} " =~ " iqtree " ]] || [[ " ${tree_method[*]} " =~ " astral " ]]
then
	echo Running ${tree_method[@]} to generate Species trees
else
	echo "Species tree estemation from $input_alignment (concatenated genes) not done; tree_method not set to either fasttree, raxml, and/or iqtree"
	USAGE
	exit 1
fi

############################################################
# print all variables set up to this point                #
# for debugging arg parsing                               #
# compare the output of declare -p at the begining and now#
# $tmpfile was created at the top of script               #
############################################################
#declare -p | diff "$tmpfile" - | grep "declare" | cut -d " " -f 4- > varables_set.txt
rm -f "$tmpfile"


#################################
# import control_file arguments #
#################################
if [[ -z ${control_file+x} ]]
then
	#nothing
	echo "All required args taken from command line"
else
	echo ""
	echo "##########################################"
	echo "  setting variable found in $control_file"
	echo "   This will overwrite any args set on the command line"
	echo "##########################################"
	echo ""
	source $control_file || exit 1
fi

#######################
## set up some stuff ##
#######################

# outputs are held in:
if [ ! -d $store ]; then
	mkdir -p $store || exit 1
fi
cd $store || exit 1

export genome_dir=$store/genomes
export wd=$store/phylo_current
export run_notes=$store/run_notes.txt
export trans=$store/annots_nucls
export prots=$store/annots_prots
export annots=$store/annots_details

# internal repo programs
export ANI_genome_picking=$script_home/python_scripts/ANI_genome_picking.py
export OG_sco_filter=$script_home/python_scripts/OG_sco_filter.py

#load custom aliases...might get rid of
shopt -s expand_aliases
source $script_home/script_lib/bash_utils_and_aliases.sh

# if present, clean-up a file detailing the runtimes for each step 
if [ -f $store/timing ]
then
	rm $store/timing
fi

# need to be here to deal with changes in ANI_shortlist in control_file.user
#  will fix later with an ANI_shortlist_frac variable (in place of "3")
export ANI_shortlist_min_OGs=`expr $ANI_shortlist / 3`

echo "
##################################################
########### Declare MAIN_PIPE to run  ############
##################################################
"
MAIN_PIPE () {
	func_timing_start
	export ANI=false
	SET_UP_DIR_STRUCTURE
	if [[ ${genomes_provided+x} ]]
	then
		PRODIGAL_PREDICT $genome_dir
	fi
	# handle user inputting their own annotations
	if [[ ${annots_provided+x} ]]
	then
		echo "trying to move gene seqs"
		#move provided prots and trans to their respecive folders
		for I in $(ls $input_prots/)
		do
			echo "moving $I"
 			cp $input_prots/${I}   $prots/${I%.*}.faa || exit
 		done
 		for I in $(ls $input_trans/)
		do
			echo "moving $I"
 			cp $input_trans/${I}   $trans/${I%.*}.fna || exit
 		done
 	fi
	DEDUP_annot_trans $store $trans $prots 
	FIX_TRANS_NAMES $trans
	FIX_PROTS_NAMES $prots

	# test if user input preannotated transcripts 
        #   or wants to use transcripts for the ANI subsetting
	if [ "$ANI_trans" = true ]
	then
		ANI_dataset=$trans
	elif [ "$ANI_genome" = true ]
	then
		ANI_dataset=$genome_dir
	else
		echo "PANIC: no dataset declared for ANI_species_shortlist to chew on"
		exit
	fi
    
	# find subset of genomes or transcripts that represents diversity of full set
	#   if number of sequence files is greater than the max number to send through orthofinder
	if [[ $(ls $ANI_dataset | wc -l)  -gt $ANI_shortlist ]]
	then
		export ANI=true
		#ANI script makes a prot directory from shortlist for orthofinder ($prots.shortlist)
		ANI_species_shortlist $ANI_dataset $ANI_shortlist
		prots4ortho=${prots}.shortlist
		ORTHO_RUN $prots4ortho
		# find genes from full set for each OG (make HMM profile and search against all prots)
		ANI_ORTHOFINDER_TO_ALL_SEQS
	else
		prots4ortho=${prots}.fixed
		ORTHO_RUN ${prots4ortho}
		REALIGN_ORTHOGROUP_PROTS
	fi

	# trim protein sequences
	TRIM $wd/AlignmentsProts "PROT"

	# identify single copy orthologs in protein alignments
	#  create files with OG names 
	if [[ "$relaxed" != false ]]
		then
			SCO_MIN_ALIGN $wd/AlignmentsProts.trm $min_num_orthos
			#ALIGNMENT_STATS $wd/SCO_${min_num_orthos}.align
	fi
	SCO_MIN_ALIGN $wd/AlignmentsProts.trm $(cat $store/all_input_list | wc -l)

	# Run Prot workflow
	if [[ "$TREE_DATA" = "PROT" || "$TREE_DATA" = "BOTH" ]]
	then
		echo "
		##########################################
		### Running Protein specific  workflow ###
		##########################################
		"
		prot_rename $wd/AlignmentsProts.trm
		if [[ "$relaxed" != false ]]
		then
			echo "$wd/SCO_$min_num_orthos"
			cat_alignments \
				$wd/SCO_$min_num_orthos \
				$wd/AlignmentsProts.trm.nm \
				$wd/AlignmentsConcatenated \
				PROT
			ALIGNMENT_STATS $wd/SCO_${min_num_orthos}.PROT.align
		fi
		ALIGNMENT_STATS $wd/AlignmentsProts.trm.nm
		cat_alignments \
			$wd/SCO_strict \
			$wd/AlignmentsProts.trm.nm \
			$wd/AlignmentsConcatenated \
			PROT
		ALIGNMENT_STATS $wd/SCO_strict.PROT.align
		for I in $wd/AlignmentsConcatenated/*.PROT.*.phy
			do
			TREE_BUILD $wd/SpeciesTree/ $I $threads "PROT"
		done
		# for later...
		#build_gene_trees
		#build_ASTRAL_trees
	fi

	# Run CDS workflow
	if [[ "$TREE_DATA" = "CDS" || "$TREE_DATA" = "BOTH" ]]
	then
		# make Codon alignments for all OGs (could just do SCO_relaxed and that would grab all SCO_strict)
		#  Its pretty fast, who cares, and might want to make SCO_very_relaxed trees
		TRIMAL_backtrans $wd $wd/AlignmentsProts.trm $wd/all_trans.nm.fa $wd/AlignmentsTrans.trm $wd/OG_names $wd/SequencesCDS
		
		# concatenate SCO alignments (only do SCO_$min_num_orthos if relaxed != false)
		if [[ "$relaxed" != false ]]
		then
			cat_alignments \
				$wd/SCO_$min_num_orthos \
				$wd/AlignmentsTrans.trm.nm \
				$wd/AlignmentsConcatenated \
				CDS
			ALIGNMENT_STATS $wd/SCO_${min_num_orthos}.CDS.align
		fi
		cat_alignments \
			$wd/SCO_strict \
			$wd/AlignmentsTrans.trm.nm \
			$wd/AlignmentsConcatenated \
			CDS
		ALIGNMENT_STATS $wd/SCO_strict.CDS.align

		# Build ML Species Treeeeees
		for I in $wd/AlignmentsConcatenated/*.CDS.*.phy
		do
			#ALIGNMENT_STATS $wd/AlignmentsTrans.trm.nm/
			TREE_BUILD $wd/SpeciesTree/ $I $threads "CDS"
		done

		# Not really useful, will probably remove (sort of redundant)
		if [ ! $ANI = true ]
		then
					# at the moment it doesnt make sense to make trees from OF if using ANI
					# (only a shortlist genomes will be in the tree)
					orthofinderGENE2SPECIES_TREE
		fi
	fi	
	
	# runs Astral workflow on CDS, PROT, or BOTH
	if [[ " ${tree_method[*]} " =~ " astral " ]]
	then
		if [[ "$TREE_DATA" = "CDS" || "$TREE_DATA" = "BOTH" ]]
		then
			tree_data=CDS
			# build gene trees for all CDS alignments
			allGENE_TREEs $tree_data $wd/AlignmentsTrans.trm.nm $wd/${tree_data}_gene_trees

			# astral_allTransGENE2SPECIES_TREE #Still not written (needs a different ASTRAL )
			
			# build ASTRAL tree from SCO_relaxed gene trees
			if [[ "$relaxed" != false ]]
				then
				astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_$min_num_orthos
			fi

			# Build ASTRAL tree from SCO_strict
			astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_strict
		fi
		if [[ "$TREE_DATA" = "PROT" || "$TREE_DATA" = "BOTH" ]]
		then
			tree_data=PROT
			# build gene trees for all CDS alignments
			allGENE_TREEs $tree_data $wd/AlignmentsProts.trm.nm $wd/${tree_data}_gene_trees 

			# astral_allTransGENE2SPECIES_TREE #Still not written (needs a different ASTRAL )
			
			# build ASTRAL tree from SCO_relaxed gene trees
			if [[ "$relaxed" != false ]]
			then
				astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_$min_num_orthos
			fi

			# Build ASTRAL tree from SCO_strict
			astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_strict
		fi
	fi	
	# clean up some files
	# need to make specific for PROT or CDS run
	WRAP_UP

	# if running tester script, compare trees to reference to make sure stuff worked
	if [[ $TESTER == "TESTER" ]]
	then
		TESTER_compare
	fi

}
##################################################
##################################################




#call pipe modules from script_lib/functions.sh via the MAIN_PIPE function declared above
MAIN_PIPE

echo "
######################################################
Weelll it finished.  Check $wd/FINAL_SPECIES_TREES/ for output
######################################################
Also check Alignment figures in the AlignmentsTrans.trm.nm.vis directory
Good luck! 
If you have an issues, please go to https://github.com/eamiddlebrook/OrthoPhyl and open an issue.

If used in a publication, please cite:
Earl A Middlebrook, Robab Katani, Jeanne M Fair 
OrthoPhyl - Streamlining large scale, orthology-based phylogenomic studies of bacteria at broad evolutionary scales
G3 Genes|Genomes|Genetics, 2024;, jkae119, https://doi.org/10.1093/g3journal/jkae119
"
