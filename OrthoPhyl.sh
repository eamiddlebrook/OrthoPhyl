#!/bin/bash

USAGE () {
echo "
USAGE: OrthoPhyl.sh -g Path_to_directory_of_assemblies -s directory_to_store_output
# ALL arguments are optional if set in control_file.required
#   Many default parameters are set in control_file.defaults
Required:
-g	path to genomes directiory
or
-a	paths to protien and transcript directories. 
       They should be delared as \"-a path_to_protien_dir,path_to_transcript_dir\"
-s	path to the main directory for output
Optional:
-t	threads to use [4]
-p     phylogenetic tree software to use [fasttree, raxml, and/or iqtree] 
	i.e. -p \"fasttree iqtree\"
-c	path to a control file with required variables and any optional ones to override defaults. 
	Will override values set on command line! [NULL]
-x	trimal paramerter string (in double \"quotes\")
-r	flag to rerun orthofinder on the ANI_shorlist (true/[false])
-n	max number of genomes to run through OrthoFinder. 
	If more than this many assemblies are provided, a subset of genomes will be chosen for OrthoFinder to chew on [20]
-T	run test dataset, incompatable with -g|s (TESTER,TESTER_chloroplast)
-h	display a description and a super useful usage message
###############################################################\n
To run test datasets:
# Test of full bacterial genomes
bash OrthoPhyl.sh -T TESTER -t #threads
# Big test with ~100 orchid chloroplasts
bash OrthoPhyl.sh -T TESTER_chloroplast -t #threads
# When running through Singularity an output directory is required:
bash OrthoPhyl.sh -T TESTER -s output_dir -t #threads
"
}


######################
# set up environment #
######################


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
    source $script_home/control_file.paths
fi


# Used to Catch all set variables later
tmpfile=$(mktemp)
declare -p > "$tmpfile"


# import pipeline defaults (will be overridden by command line args or -c control_file)
source $script_home/control_file.defaults


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
        ##########   fasttest Edition!!!!   #########
	################################################
        "
	# NEED TO ADD COMPRESSED GENOME FILES TO TEST NEW FUNC
	source $script_home/TESTER/control_file.user_fasttest
	if [ -d $store ]
       then
              rm -r $store
       fi
else
	echo "PANIC: TESTER was set, but didnt equal 'TESTER' or 'TESTER_chloroplast'" 
	exit
fi
}


###################################################
#######  Handle command line Arguments ############
###################################################
#Setting Variables manually to override ones in control file
ARGS_SET=""
while [[ ${1:0:1} = '-' ]] ; do
N=1 
L=${#1} 
while [[ $N -lt $L ]] ; do 
  case ${1:$N:1} in 
     'h') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE
            exit 0
          fi 
          # run function to print usage and exit
	  echo "Looks like you used the '-h' with an argument...thats not right" 
	  USAGE
          exit 1
          shift ;; 

     'g') if [[ $N -ne $(($L-1)) || ! -n ${2} ]]; then 
            echo "looks like the argument for -g is messed up" 
            echo $USAGE
            exit 1
          fi
          input_genomes="$( cd "$(relative_absolute ${2})" && pwd )"
          genomes_provided="TRUE"
          if [[ ! -d "${input_genomes}" ]]; then
            echo "WARNING: -g declaring an input genome directory that does not exist...maybe check on that"
            exit 1
          fi
	  ARGS_SET+=g
	  shift ;;

     's') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE 
            exit 1 
          fi
          store=$(relative_absolute ${2})
          ARGS_SET+=s
          shift ;;

     'p') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE 
            exit 1 
          fi
          tree_method=(${2})
          echo "#####"
          echo ${2}
          echo "#####"
          ARGS_SET+=p
          shift ;;

     'a') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then
            USAGE
            exit 1
          fi
          annots_provided="TRUE"
          prots=$(awk  -F ',' '{print $1}' <<< ${2})
          trans=$(awk  -F ',' '{print $2}' <<< ${2})
	   # deal with relative or absolute paths correctly
	   input_prots="$( cd "$(relative_absolute ${prots})" && pwd )"
	   ls $input_prots
          input_trans="$( cd "$(relative_absolute ${trans})" && pwd )"
          if [[ ! -d "${input_prots}" ]] || [[ ! -d "${input_trans}" ]] ; then
            echo "WARNING: -a declared input annoation directories that do not exist...maybe check on that"
            echo " They should be delared as \"-a path_to_protien_dir,path_to_transcript_dir\""
            exit 1
          fi
          ARGS_SET+=a
          shift ;;

     't') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then
            USAGE
            exit 1
          fi
          threads="${2}"
          ARGS_SET+=t
          shift ;;

     'c') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
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
	  shift ;;

     'x') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE
            exit 1
          fi
          trimal_parameter=${2} 
          ARGS_SET+=x
          shift ;;

     'r') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE
            exit 1 
          fi 
          orthofinder_rerun=${2} 
          ARGS_SET+=r
          shift ;;

     'n') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then 
            USAGE
            exit 1
          fi
          ANI_shortlist=${2} 
          ARGS_SET+=a
          shift ;;

     'T') if [[ $N -ne $(($L-1)) || ! -n ${2} ]] ; then
	    USAGE
	    exit 1
	  fi
	export TESTER=${2}
        ARGS_SET+=T
	# set variables specific for running OrthoPhyl on a test data set.
	tester
	shift ;;

     \?) # Invalid option
         echo "Error: Invalid option"
         USAGE
         exit;;
     *) USAGE 
        exit 1 ;; 
  esac 
  N=$(($N+1)) 
done 
shift 
done 
if [[ -n ${1} ]] ; then 
	echo "!!!!: you have a hanging argument without a flag!"
	USAGE 
	exit 1 
fi 
# test for incompatable args (if) and makes sure required args are present (elif)
if [[ "$ARGS_SET" == *@(g|c|a)*@(T)* ]] || [[ "$ARGS_SET" == *@(T)*@(g|c|a)* ]]
then
	echo "!!!!: -T was set along with -g/c/a, which are incompatable args"
	USAGE
	exit 1
elif [[ ! "$ARGS_SET" == *c* ]] && [[ ! "$ARGS_SET" == *@(g|s)*@(g|s)* ]] && [[ ! "$ARGS_SET" == *T* ]] && [[ ! "$ARGS_SET" == *@(a|s)*@(a|s)* ]]
then
	echo -e "WARNING: Required arguments were not given, you need to provied either \n\ta control file with \"-c control_file\"\n\t-g genome_directory -s storage_directory \n\t-n prot_dir,trans_dir -s storage_directory \n\t-g genome_directory -n prot_dir,trans_dir -s storage_directory  \n\tor \"-T TESTER\" to start a test run"
	echo "Arguments that you set are -"${ARGS_SET}
	USAGE
	exit 1

fi

# test if "-p tree_method" is set correctly
if [[ " ${tree_method[*]} " =~ " raxml " ]] || [[ " ${tree_method[*]} " =~ " fasttree " ]] || [[ " ${tree_method[*]} " =~ " iqtree " ]]
then
	echo Running $tree_method to generate Species trees
else
	echo "Species tree estemation from $input_alignment (concatenated genes) not done; tree_method not set to either fasttree, raxml, and/or iqtree"
	USAGE
	exit 1
fi

# print all variables set up to this point
# for debugging arg parsing
# compare the output of declare -p at the begining and now
# $tmpfile was created at the top of script
#declare -p | diff "$tmpfile" - | grep "declare" | cut -d " " -f 4- > varables_set.txt
rm -f "$tmpfile"



#import control_file arguments
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


# outputs are held in:
mkdir -p $store || exit 1
cd $store || exit 1

export genome_dir=$store/genomes
export wd=$store/phylo_current
export run_notes=$store/run_notes.txt
export trans=$store/prodigal_nucls
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
	if [[ ${genomes_provided+x} ]] 
	then
		PRODIGAL_PREDICT $genome_dir
	fi
	# handle user inputting their own annotations
	if [[ ${annots_provided+x} ]]
	then	
		echo "trying to move genes seqs"
		#move provided prots and trans to thier respecive folders
		for I in $(ls $input_prots/)
		do
 			cp $input_prots/${I}   $prots/${I%.*}.faa || exit
 		done
 		for I in $(ls $input_trans/)
		do
 			cp $input_trans/${I}   $trans/${I%.*}.fna || exit
 		done
 	fi
	DEDUP_annot_trans
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

echo "
######################################
Weelll it finished.  I doubt it worked
######################################
"
