#!/bin/bash
####################################################################
# USAGE:
# $ bash orthophylo.XXX.sh
# to test install run:
# $ bash orthophylo.XXX.sh TESTER
# to run through slurm scheduler:
# $ sbatch orthophylo.XXX.slurm TESTER
# all variable can be found in the files control_file.*
# control_file.required contains variable that must be set for a custom run
#   (TESTER requires no editing of the control files unless you have a different conda env)
# control_file.user allows the user to declare variable to override control_file.defaults
#   these can be copy/pasted directly from the defaults file and edited as desired
###################################################################
# This script is the main wrapper for the orthophylo repo
# the goal of this script is to take in a directory of bacterial genomes (in the hundreds)
# and generate species trees.
# Trees generated from both concatenated (protien informed) transcript alignments
# and gene tree to species tree methods (ASTRAL)
##
# future versions will include evolutionary model testing of indevidual orthogroups with ete3


######################
# set up environment #
######################

# sometimes reqired by HPC submission software
# also on a mac ~/.bash*s are not sourced automatically
# (just covering all the bases)
source $HOME/.bash_profile
source $HOME/.bashrc

#used for loading stuff from git repo
export script_home=$(pwd)
TESTER=""
if [[ ! -z ${1+x} ]]
then
	export TESTER="$1"
	echo $TESTER
fi

############################################
####  import control file variables  #######
############################################

# Will be overwritten if provided on the command line

# import required user defined variables
source $script_home/control_file.required

# import pipeline defaults
source $script_home/control_file.defaults

# import used defined parameters to override defaults
source $script_home/control_file.user

# import main_script functions
source $script_home/script_lib/functions.sh


# if first arg is "TESTER" run a small pipeline to test...pipeline
if [[ $TESTER == "TESTER" ]]
then
	echo "
	################################################
	#####  Testing Workflow with Control Files  ####
	####  and genomes from the TESTER directory  ###
	################################################
	"
	source $script_home/TESTER/control_file.required
	source $script_home/control_file.defaults
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
       	##########   Chloroplast Edition!!!!   ##########
	################################################
        "
	# NEED TO ADD COMPRESSED GENOME FILES TO TEST NEW FUNC
	source $script_home/TESTER/control_file.required_chloroplast
        source $script_home/control_file.defaults
        source $script_home/TESTER/control_file.user
        if [ -d $store ]
        then
                rm -r $store
        fi
else


	###################################################
	#######  Handle command line Arguments ############
	###################################################
	DESCRIPTION="OrthoPhylo makes trees and stuff...link this to an external file"
	USAGE="orthophyl.sh -hgstxrA
	# all arguments are optional if set in control_file.required
	-h 	display a description and a super useful usage message
	-g	full path to genomes directiory
	-s	full path to the main directory for output
	-t	threads to use
	-x	trimal paramerter string (in double "quotes")
	-r	flag to rerun orthofinder on the ANI_shorlist (true/false)
	-a	max number of genomes to run through OrthoFinder. If more than this many assemblies are profided, a subset of genomes will be chosen for OrthoFinder to chew on

	"

	#Setting Variables manually to override ones in control file
	while [[ ${1:0:1} = '-' ]] ; do 
	N=1 
	L=${#1} 
	while [[ $N -lt $L ]] ; do 
	  case ${1:$N:1} in 
	     'h') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 'g' 
	            exit 1 
	          fi 
	          echo $DESCRIPTION
	          echo $USAGE
	          exit $1
	          shift ;; 

	     'g') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo "looks like the argument for -g is messed up" 
	            echo $USAGE
	            exit 1 
	          elif [[ ! -d {2} ]]; then
	            echo "WORNING: -g declaring an input genome directory that does not exist...maybe check on that"
	            exit 1
	          fi 
	          input_genomes=${2} 
	          shift ;; 

	     's') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 's' 
	            exit 1 
	          fi 
	          store=${2} 
	          shift ;; 
	    
	     't') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 't' 
	            exit 1 
	          fi 
	          threads=${2} 
	          shift ;; 
	     
	     's') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 's'
	            exit 1 
	          fi 
	          store=${2} 
	          shift ;; 
	  
	     'x') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 'x' 
	            exit 1 
	          fi 
	          trimal_parameter=${2} 
	          shift ;; 
	     'r') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 'r' 
	            exit 1 
	          fi 
	          orthofinder_rerun=${2} 
	          shift ;; 
	     'a') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
	            echo $USAGE 'A' 
	            exit 1 
	          fi 
	          ANI_shortlist=${2} 
	          shift ;; 
	     \?) # Invalid option
	         echo "Error: Invalid option"
	         echo $USAGE
	         exit;;
	     *) echo $USAGE 't2' 
	        exit 1 ;; 
	  esac 
	  N=$(($N+1)) 
	done 
	shift 
	done 
	if [[ ! -n ${1} ]] ; then 
	echo $USAGE 
	exit 1 
	fi 
fi

# outputs are held in:
mkdir $store
export genome_dir=$store/genomes
export wd=$store/phylo_current
export run_notes=$store/run_notes.txt
export genome_list=($(ls $genome_dir))
cd $store || exit

export trans=$store/prodigal_nucls/
export prots=$store/prodigal_prots
export annots=$store/prodigal_annots

# internal repo programs
export ANI_genome_picking=$script_home/python_scripts/ANI_genome_picking.py
export OG_sco_filter=$script_home/python_scripts/OG_sco_filter.py

#load custom aliases...might get rid of
shopt -s expand_aliases
source $script_home/script_lib/bash_utils_and_aliases.sh

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
	PRODIGAL_PREDICT $genome_dir
	DEDUP_annot_prots
	FIX_TRANS_NAMES $trans
	FIX_PROTS_NAMES $prots
	# find subset of genomes ($ANI_shortlist) that represents diversity of full set
	if [[ $(cat $store/genome_list | wc -l)  -gt $ANI_shortlist ]]
	then
		export ANI=true
		#ANI script makes a prot directory from shortlist for orthofinder ($prots.shortlist)
		ANI_species_shortlist $genome_dir $ANI_shortlist
		prots4ortho=${prots}.shortlist
		ORTHO_RUN $prots4ortho # do not comment out, robust to reruns
		# find genes from full set for each OG (make HMM profile and search against all prots)
		ANI_ORTHOFINDER_TO_ALL_SEQS
	else
		prots4ortho=${prots}.fixed
		ORTHO_RUN ${prots4ortho}
		REALIGN_ORTHOGROUP_PROTS
	fi
	PAL2NAL
	TRIM_TRANS
	ALIGNMENT_STATS $wd/AlignmentsTrans.trm.nm/
	SCO_MIN_ALIGN $min_num_orthos
	ALIGNMENT_STATS $wd/OG_SCO_${min_num_orthos}.align
	SCO_strict
	ALIGNMENT_STATS $wd/OG_SCO_strict.align
	for I in $wd/SpeciesTree/*.trm.sco.nm.phy
        do
		TREE_BUILD $wd/SpeciesTree/ $I $threads
	done
	if [ ! $ANI = true ]
        then
                # at the moment it doesnt make sense to make trees from OF if using ANI
                # (only a shortlist genomes will be in the tree)
                orthofinderGENE2SPECIES_TREE
        fi
	allTransGENE_TREEs
	# astral_allTransGENE2SPECIES_TREE #Still not written (needs a different ASTRAL )
	astral_TransGENE2SPECIES_TREE $wd/OG_SCO_$min_num_orthos
	astral_TransGENE2SPECIES_TREE $wd/OG_SCO_strict
	WRAP_UP
	if [[ $TESTER == "TESTER" ]]
	then
		TESTER_compare
	fi

}
##################################################
##################################################




#call pipe modules via the MAIN_PIPE function declared at the top of script (for convieniance)
MAIN_PIPE

echo '
######################################
Weelll it finished.  I doubt it worked
######################################
'

