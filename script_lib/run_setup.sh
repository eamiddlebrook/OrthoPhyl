#!/bin/bash


#######################
##### arg parser ######
#######################
ARG_PARSE () {
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
}

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

####################################################################################
# test for incompatable args (if) and makes sure required args are present (elif) ##
####################################################################################
test_args () {
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

}