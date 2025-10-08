#!/bin/bash
echo "
Running...

   ____       __  __                __          __   
  / __ \_____/ /_/ /_  ____  ____  / /_  __  __/ /   
 / / / / ___/ __/ __ \/ __ \/ __ \/ __ \/ / / / /    
/ /_/ / /  / /_/ / / / /_/ / /_/ / / / / /_/ / /     
\____/_/   \__/_/ /_/\____/ ,___/_/ /_/\__, /_/      
                         /_/          /____/     
				    ____       __               ____
				   / __ \___  / /   ___  ____  / __/
				  / /_/ / _ \/ /   / _ \/ __ '/ /_  
				 / _, _/  __/ /___/  __/ /_/ / __/  
				/_/ |_|\___/_____/\___/\__,_/_/     
                                    
"
######################
# set up environment #
######################
echo "Loading required files..."

# grab the Orthophyl directory no matter where script was run from
#used for loading stuff from git repo
export script_home=$(dirname "$(readlink -f "$0")")
source $HOME/.bash_profile
source $HOME/.bashrc
source $script_home/control_file.defaults || { echo "cannot open control_file.defaults. This file is required for telling OrthoPhyl what parameters to use if not provided on the CMD line"; exit 1 ;}
source $script_home/control_file.paths || { echo "cannot open control_file.paths. This file is required for telling OrthoPhyl where some exteranl programs are."; exit 1 ;}


#############################################
######## Import funtions and stuff ##########
#############################################
# import main_script functions
source $script_home/script_lib/functions.sh
source $script_home/script_lib/run_setup.sh


source $script_home/script_lib/functions_addem.sh
source $script_home/script_lib/arg_parse_addem.sh

#load custom aliases...might get rid of
shopt -s expand_aliases
source $script_home/script_lib/bash_utils_and_aliases.sh

if [[ ${##[@]} -lt 2 ]]
then 
	USAGE
	echo ""
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "No arguments provided to ReLeaf"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	echo "See GitHub page for useful tips"
	echo "github.com/eamiddlebrook/OrthoPhyl"
	exit 1
fi

###################################################
### Parameters (replace with arg parce section) ###
###################################################

ARG_PARSE_addem "$@"
#test_args_addem


# instead of $store, you can provide the files/directories directly
# tree=
# tree_dir=
# hmm_dir=
# OG_alignment_dir=
# iqtree_best_scheme



# maybe stuff needed for adding to fasttree
fasttree_CDS_model=GTR_gamma
fasttree_PROT_model=???
total_tree_rerun=FALSE

#########################################
### find pre-computed stuff in $store ###
#########################################
#Could speed this up by only searching for OGs in the SCO_sets

if [ -d "$store/OG_alignmentsToHMM/hmm_round2/hmms/" ]
then
	hmm_dir="$store/OG_alignmentsToHMM/hmm_round2/hmms/"
	ls $hmm_dir/*.hmm > /dev/null || echo "NO hmm files found in "$hmm_dir
	echo "Did the original OrthoPhyl run complete?"
elif [ -d "$store/annots_prots.fixed/OrthoFinder/Results_ortho/MultipleSequenceAlignments/" ]
then
	# to use later
	hmm_dir="$store/annots_prots.fixed/OrthoFinder/Results_ortho/MultipleSequenceAlignments/"
	build_hmm="True"
else
	echo -e "NO OLD HMM DIRECTORY FOUND within the provided OrthoPhyl run directory:\n$store "
	echo "Cannot find $store/OG_alignmentsToHMM/hmm_round2/hmms/"
	exit 1
fi
#partition_file="$store/add_seqs/partition_file.blah.blah"
OLD_prot_alignments="$store/phylo_current/AlignmentsProts"
OLD_prot_alignments_trm="$store/phylo_current/AlignmentsProts.trm"
OLD_CDS_alignments="$store/phylo_current/AlignmentsCDS"
OLD_CDS_alignments_trm="$store/phylo_current/AlignmentsCDS.trm"
OLD_trim_cols="$store/phylo_current/trimmed_columns"
OLD_iqtree_output="$store/phylo_current/SpeciesTree/iqtree/"
OLD_TREES=$store/FINAL_SPECIES_TREES/

##########################
#### output locations ####
##########################
addasm_dir=$store/ReLeaf_dir
wd=$addasm_dir
all_prots=$addasm_dir/all_prots.nm.fa
all_CDS=$addasm_dir/all_CDS.plus_old.nm.fa
genome_dir=$addasm_dir/genomes
prots=$addasm_dir/annots_prots
trans=$addasm_dir/annots_CDS
annots=$addasm_dir/annots
new_hmm_output=$addasm_dir/hmm_out
new_OG_prots=$addasm_dir/new_OG_prots
new_OG_CDS=$addasm_dir/new_OG_CDS
new_prot_alignments=$addasm_dir/new_prot_alignments
new_CDS_alignments=$addasm_dir/new_CDS_alignments
new_trees=$addasm_dir/new_trees
OG_names=$addasm_dir/OG_names
OLD_TREE_INFO=$addasm_dir/OLD_TREE_INFO
OLD_TREE_LIST=$addasm_dir/OLD_TREE_LIST
SpeciesTree=$addasm_dir/SpeciesTree
AlignmentsConcatenated=$addasm_dir/AlignmentsConcatenated
logs=$addasm_dir/logs

################################
### Set variable ###############
################################
# set some user specified stuff variable 
# !!!!! this is ignoreed until ASTRAl block...need change it to go off of available trees.
if [[ relaxed=="False" ]]
then
	relaxed="False"
elif [[ relaxed=="True" && -f $wd/OG_SCO_$min_frac_orthos ]]
then
	relaxed=="True"
elif [ -f $wd/OG_SCO_$min_frac_orthos ]
then
	relaxed="True"
fi


MAIN_PIPE () {
	func_timing_start
    ### 
    #SET_UP_DIR_STRUCTURE
	# use pre made OrthoPhyl $wd ($store/phylo_current)
	GET_OLD_TREE_INFO
	SET_UP_ADDASM_DIRS
	CLEAN_N_COPY_GENOMES \
		$genome_dir \
		$input_genomes
	# get/organize annotations
	ANNOTATIONS
	# search for orthologs in new annots based on old HMMs
	HMM_search \
		$new_hmm_output \
		$hmm_dir \
		$all_prots \
		$new_OG_prots
	
	# if running CDS ML or astral workflow, get the backtranslation
	if [[ " ${TREE_DATA_list[*]} " =~ " CDS " ]]
	then
		GET_OG_NAMES $wd/OG_names $new_OG_prots 
		filterbyname_subfunc $wd/OG_names $addasm_dir/all_trans.nm.fa $new_OG_CDS
		ADD_2_ALIGNMENTS \
			$OLD_CDS_alignments \
			$new_OG_CDS \
			$new_CDS_alignments \
			$threads
		
		# done in original OP run
		#GET_TRIMMED_COLS \
		#	$OLD_trim_cols
		#	$new_trim_cols
		#cat $addasm_dir/all_trans.nm.fa $store/phylo_current/all_trans.nm.fa \
		#	> $all_CDS
		#TRIMAL_backtrans \
		#	$addasm_dir \
		#	$new_prot_alignments \
		#	$all_CDS \
		#	$new_CDS_alignments \
		#	$OG_names \
		#	$new_OG_CDS \
		#	true		
		TRIM_selectcols $new_CDS_alignments "CDS" $OLD_trim_cols $wd/TRIM_selectcols_CDS.complete \
			&& touch $wd/TRIM_selectcols_CDS.complete
		CDS_rename $new_CDS_alignments.trm

	fi

	if [[ " ${TREE_DATA_list[*]} " =~ " PROT " ]]
	then 
		# add these OGs to old alignments
		ADD_2_ALIGNMENTS \
			$OLD_prot_alignments \
			$new_OG_prots \
			$new_prot_alignments \
			$threads 
		TRIM_selectcols $new_prot_alignments "prot" $OLD_trim_cols complete && touch $wd/TRIM_selectcols_prot.complete
		prot_rename ${new_prot_alignments}.trm
	fi

	
	
	# run standard ML workflow on PROT and/or CDS supermatrix
	if [[ " ${tree_method[*]} " =~ " raxml " || " ${tree_method[*]} " =~ " iqtree " || " ${tree_method[*]} " =~ " fasttree " ]]
	then
		RUN_ML_SUPERMATRIX_WORKFLOW 
		# takes global variables...ugh
		# $new_prot_alignments_trm and $new_CDS_alignments among others
	fi

	echo -e "#########\ngot to end\n##########" && exit

	# run ASTRAL workflow on PROT and/or CDS trees
	if [[ " ${tree_method[*]} " =~ " astral " ]]
	then
		ASTRAL_WORKFLOW
	fi
	# clean up some files
	# need to make specific for PROT or CDS run
	WRAP_UP

}
MAIN_PIPE