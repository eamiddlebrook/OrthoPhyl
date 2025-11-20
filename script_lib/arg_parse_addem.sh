#!/bin/bash

USAGE () {
echo "
USAGE: ReLeaf.sh -g Path_to_directory_of_assemblies -a path_to_CDS_dir,path_to_PROTS_dir -s Previous_OP_output
***Check out github.com/eamiddlebrook/OrthoPhyl for lots of details***
# ALL arguments are optional if set with \"-c control_file.your_args\"
#   Many default parameters are set in control_file.defaults
Required:
-g|--genome_dir  path to genomes directiory
or
-a|--annotations  paths to protien and transcript directories.
       They should be delared as \"-a path_to_CDS_dir,path_to_prot_dir\"
-s|--storage_dir  path to the main directory output from original OrthoPhyl run
Optional:
-t|--threads  threads to use [4]
-p|--phylo_tool  phylogenetic tree software to use astral, fasttree, and/or iqtree [\"fasttree iqtree astral\"]
	i.e. -p \"fasttree iqtree astral\"
	Default: Will be taked from trees available in storage_dir
-o|--omics  "omics" data to use for tree building ([CDS], PROT, BOTH)
	for divergent sequences, it is good to compare protein trees to 
	nucleotide trees to identify artifacts of saturation (long branch attraction)
	Default: Will be taked from trees available in storage_dir
-h|--help  display a description and a super useful usage message
"
}

#######################
##### arg parser ######
#######################
ARG_PARSE_addem () {
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
			-m|--min_frac_orthos) if [[ ! -n ${2} ]] ; then
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
					min_frac_orthos=${2}
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

			-x|--trimal_params) if [[ ! -n ${2} ]] ; then 
					USAGE
					exit 1
				fi
				trimal_parameter=${2} 
				ARGS_SET+=x
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
echo $ARGS_SET
echo $tree_method
}
