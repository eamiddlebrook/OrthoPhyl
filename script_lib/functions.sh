#!/bin/bash

# functions used in orthophylo.sh
# i.e.
# source $script_home/scipt_lib/functions.sh

BAIL () {
	# print function that had error in it and any associated message ${1}
	echo -e "WARNING: Function \"${FUNCNAME[2]}\" failed"
	echo ${1}
	exit 1
}

SET_UP_DIR_STRUCTURE () {
	echo '
       	###################################################
       	######### Setting up Directory structure  #########
       	########## and Declaring some variables ###########
        ###################################################
	'
	date
	func_timing_start

	if [ -f $store/setup.complete ]
	then
		echo "Output folder setup seems to be already completed"
		echo "   Resetting min_frac_orthos in case OrthoPhyl crashed between Setup and SCO_min_frac_orthos"
		export min_frac_orthos=$(printf %.0f $(echo "$(cat $store/all_input_list | wc -l)*$min_frac_orthos" | bc))
		echo " If this is not correct please rm $store/setup.complete and restart OrthoPhyl"
	else
		# set up annotation output directories
		if [ ! -d $store ] ; then
			mkdir $store
		fi
		if [ ! -d $trans ] ; then
			mkdir $trans 
		fi
		if [ ! -d $prots ] ; then
			mkdir $prots
		fi 
		if [ ! -d $annots ] ; then
			mkdir $annots
		fi

		CLEAN_N_COPY_GENOMES $genome_dir $input_genomes

		#make a file with the genome and protein file names (if input) for later use
		if compgen -G "$genome_dir/*.*a" > /dev/null
                then
                    	ls > $store/genome_list
			genome_list=$store/genome_list
                fi
		if [[ ${annots_provided+x} ]]
		then
			ls $input_prots > $store/pre_annotated_list
			preannotated_list=$store/pre_annotated_list
		fi
		# make final input list
		cat $genome_list $preannotated_list > $store/all_input_list

		# calculate and report how many input assemblies there are
		echo "Building phylogeny for $(cat $store/all_input_list | wc -l) assemblies" | tee -a $run_notes
		export min_frac_orthos=$(printf %.0f $(echo "$(cat $store/all_input_list | wc -l)*$min_frac_orthos" | bc))
		echo "All Single copy orthologs and SCOs represented in $min_frac_orthos will be used to create species trees with ${tree_method[@]}" | tee -a $run_notes
		###########################################
		# just setting up the directory structure #
		###########################################
		mkdir $wd
		cd $wd || exit
		mkdir AlignmentsProts/
		mkdir AlignmentsCDS/
		mkdir AlignmentsCDS.trm/
		mkdir AlignmentsCDS.trm.nm/
		mkdir SequencesCDS/
		mkdir SequencesProts/
		mkdir OG_names/
		mkdir RAxMLtrees/
		mkdir SpeciesTree
		mkdir logs
		mkdir trimmed_columns
		mkdir OG_names

		touch $store/setup.complete
	fi
}

# functionalized for use in addasm
CLEAN_N_COPY_GENOMES () {
	echo '
	#######################################################
	############# Clean up and Copy Assemblies ############
	################ To Working Directory  ################
	#######################################################
	'
	genome_dir=$1
	input_genomes=$2
	#######################################################################
	#### Make a directory for links to genomes to place in a phylogeny
	#### fix names that contain \( or \)
	#### Can add other chars to fix as needed
	if [ -d $genome_dir ]
	then
			rm -r $genome_dir
	fi
	mkdir $genome_dir
	# If genomes are gzipped, uncompress into working dir
	if compgen -G "$input_genomes/*.gz" > /dev/null
	then
		for I in $(ls "${input_genomes}"/*.gz)
		do
			base=$(basename ${I%.*})
			gunzip -c $I > $genome_dir/$base
		done
	fi
	# make symbolic links for input genomes in $genome_dir
	if compgen -G "$input_genomes/*.fna" > /dev/null
	then
			ln -s $input_genomes/*.fna $genome_dir/ # make links to input genome dir (should switch to cp'n all the files of m$
	fi
	if compgen -G "$input_genomes/*.fa" > /dev/null
	then
			ln -s $input_genomes/*.fa $genome_dir/
	fi
	if compgen -G "$input_genomes/*.fasta" > /dev/null
		then
				ln -s $input_genomes/*.fasta $genome_dir/
		fi
	cd $genome_dir || exit
	#remove parentheses from genomes names
	#could add additional lines subbing out the "(" for the character to replace
	# get rid of [(,),:,=] in genome file Names
	for I in $(ls ./ | grep \)); do   mv $I ${I//\)/}; done
	for I in $(ls ./ | grep \(); do   mv $I ${I//\(/}; done
	for I in $(ls ./ | grep \:); do   mv $I ${I//\:/}; done
	for I in $(ls ./ | grep "_\=_" ); do mv $I ${I//_\=_} ; done
	# fix contig names with a "|" or ":" that break everything post annotion
	for I in $(grep -lm 1 -d recurse ./ -e "|");do sed -i 's/|/_/g' $I; done
	for I in $(grep -Rlm 1 -e ":");do sed -i 's/:/_/g' $I; done
}

PRODIGAL_PREDICT () {
	echo '
	###################################################
	############# Annotate Genomes For ################
	################ Phylogenetics ####################
	###################################################
	'
	date
	func_timing_start
	# pull genomes from $genome_dir and annotate. Prots and trans are left in their respective directories.
	# Added (-m) to stop gene models from crossing gaps in assemblies (and creating prots with poly Xs)
	# could expose this variable....
	local_genome_dir=$1
	cd $local_genome_dir || exit
	# need to switch this to use genome_list...it creates a bunch of *.suffix files if either .fna or .fa files are absent
	# enumerate list of genomes to annotate, wait if $threads are already running
	for genome in $(ls $local_genome_dir)
	do
		if test "$(jobs | wc -l)" -ge $threads
		then
		        wait -n
				trap control_c INT
		fi
		if [ -f $annots/"${genome%.*}".complete ]
		then
			echo "${genome%.*} already annotated"\
			> $annots/notes
		else
			{
			prodigal $prodigal_options \
			 -i ${genome} \
			 -o $annots/${genome%.*}.genes \
			 -a $prots/${genome%.*}.faa \
			 -d $trans/${genome%.*}.fna \
			 && touch $annots/${genome%.*}.complete
			} &
		fi
	done 
	wait
}
: << 'COMMENT'
CHECK_ANNOTS () {
	echo "###################################"
	echo "#### Check provided annotaions ####"
	echo "##### For compatability with ######"
	echo "####### downstream methods ########"
	echo "###################################"
	
	bad_chars=$script_home/gen_files/bad_chars
	annot_char_check=$store/annot_char_check.complete
	annot_name_check=$store/annot_name_check.complete
	wd=$store/annot_check.tmp
	char_check_failed=$wd/char_check_failed
	name_check_failed=$wd/name_check_failed
	uniq_check_failed=$wd/uniq_check_failed
	mkdir $wd
	cd $wd || exit 1
	for I in $(cat $store/pre_annotated_list | sed 's/.faa//g')
	do
		# declare CDS and PROT files 
		seq1_CDS=$annots_nucls/${I}.fna
		seq1_PROT=$annots_prots/${I}.faa
		
		# check for special characters in annots
		bad_CDS=$(cat $seq1_CDS | grep -v ">" | grep -f $bad_chars  | head -n 1)
		bad_PROT=$(cat $seq1_PROT | grep -v ">" | grep -f $bad_chars  | head -n 1)
		if [ $(echo $bad_CDS) | wc -l) -gt 0 ]
		then
			echo "${I%.*} CDS file has a messed characters (either . or *)." >> $char_check_failed
		fi
		if [ $(echo $bad_PROT) | wc -l) -gt 0 ]
		then
			echo "${I%.*} PROT file has a messed characters (either . or *)." >> $char_check_failed
		fi

		# check for identical sequence names
		cat $seq1_CDS | grep ">" | awk '{print $1}' | sort > CDS_names
		cat CDS_names >> all_CDS_names
		cat $seq1_PROT | grep ">" | awk '{print $1}' | sort > PROT_names
		cat PROT_names >> all_PROT_names
		# compare files 
		if [ ! $(comm -12 $seq1_CDS_names $seq1_PROT_names | wc -l) -eq $(cat $seq1_CDS_names | wc -l) ]
		then
			echo $seq1_CDS " CDS and PROT names dont match" >> $name_check_failed
			exit 1
		fi
		
	if [ ! -f name_check.failed ]
	then
		touch $annot_name_check
	else
		echo "PANIC: CDS and PROT names do not match"
		echo "  For details check $name_check_failed"
	fi
	if [ ! -f char_check.failed ]
	then
		touch $annot_char_check
	else
		echo "PANIC: Provided CDS or PROT files have illegal character (* or .)"
	fi
	 
	 # check that all CDS/PROT names are uniq...
	if [ $(cat all_names | wc -l ) -eq $(cat all_names | sort | uniq | wc -l) ]
	then
		echo "All CDS names look uniq, awesome"
	else
		echo "PANIC: CDS seq names are not uniq."
		echo "	This will wreck stuff. Aborting "
	fi

	if [ -f $name_check_failed ] || [ -f $char_check_failed ] || [ -f $uniq_check_failed ]
	then
		echo "EXITING: Annotations check failed"
		exit 1
	done
	fi
	# delete this
	echo "got to end of annot check"
	exit
}
COMMENT

CHECK_GENOME_QUALITY () {
	echo "########################"
	echo "#### RUNNING CHECKM ####"
	echo "########################"

	max_genomes=200
	checkM_out=$store/checkM_out
	mkdir $checkM_out
	cd $checkM_out || exit
	# split assemblies into different directories
	J=0
	K=0
	for I in $(ls $prots/*.faa)
	do
		if [ $((J % max_genomes)) -eq 0 ]
		then
			K=$((K+1))
			mkdir $checkM_out/prots_${K}
			mkdir $checkM_out/prots_${K}_out
			cd $checkM_out/prots_${K}
		fi
		# make simlinks for assembly/proteome subsets
		ln -s $I ./
		J=$((J+1))
	done

	# run checkM on each proteome subset
	cd $checkM_out
	J=1
	while [ $J -le $K ]
	do
		in=$checkM_out/prots_${J}
		out=$checkM_out/prots_${J}_out
		checkm taxonomy_wf -t $threads -x faa $in $out
		J=$((J+1))
	done

	# Some kinda filtering....
}

DEDUP_annot_trans () {
        echo "
	###################################
	#### RUNNING annotation dedupe ####
	###################################
	"

	local store=$1
	local trans=$2
	local prots=$3
	identity=99.9
	
	cd $store || exit
	if [ -f $store/dedupe.complete ]
       	then
    		echo "Dedupe previously completed."
        	echo "if rerun is desired, delete $store/dedupe.complete"
	else
		mkdir $trans.preDedup || exit 1
		mv $trans/*.f*a $trans.preDedup/
		mkdir $prots.preDedup || exit 1
        mv $prots/*.f*a $prots.preDedup/
		:> dedupe.stats
		:> dedupe.stats_long
		for I in $(ls $trans.preDedup/*.f*a)
		do
			# cleaner multithreading. check if the number of running jobs is greater than $threads
			# if so wait for one to finish
			if test "$(jobs | wc -l)" -ge $threads
			then
			        wait -n
			fi
			{
			base=$(basename ${I%.*})
			echo $base >> dedupe.stats
			bash dedupe.sh -Xmx1g -Xms1g in=$I out=$trans/$base.fna minidentity=99.9 2>&1 | \
				grep 'Input\|Result'  \
				>> dedupe.stats || { echo "dedupe failed for some reason :("; exit 1;}
			cat $trans/$base.fna | grep ">" | sed 's/>//g' > $trans/$base.deduped.names
			filterbyname.sh -Xmx1g -Xms1g --amino include=true \
				in=$prots.preDedup/$base.faa \
				out=$prots/$base.faa \
				names=$trans/$base.deduped.names \
				>> dedupe.stats_long 2>&1 || { echo "dedupes filterbyname failed for some reason :("; exit 1;}
			} 
		done && touch dedupe.complete
		#wait
	fi
	# remove the files telling filterbyname which prots to keep
	rm $trans/*.deduped.names
}

FIX_TRANS_NAMES () {
	echo '
	################################################
	############# fixing trans names ###############
	############## and catting to $wd ##############
	################################################
	'
	date
	func_timing_start
	echo "Fixing names and catting trans for later filterbyname"
	local_trans=$1
	cd $local_trans/ || exit
	for I in *.fna
	do
		base=${I%.*}
		cat "$I" | sed "s/>/>$base@/g" | sed 's/|.*//g'
	done \
	> $wd/all_trans.nm.fa
}


FIX_PROTS_NAMES () {
	echo '
	########################################
	######## Fixing prot names for #########
	######## for later filterbyname ########
	########################################
	'
	date
	func_timing_start
	local_prots=$1

	if [ -d $prots.fixed ]
	then
	        if [ -d $prots.fixed.bk ]
			then
				rm -r $prots.fixed.bk
			fi
	    mv $prots.fixed $prots.fixed.bk
	fi

	cd $local_prots/ || exit
	mkdir $prots.fixed
	for I in *.faa
	do
		base=${I%.*}
	        cat "$I" | sed "s/>/>$base@/g" | sed 's/|.*//g' \
		> $prots.fixed/$base.faa

		cat "$I" | sed "s/>/>$base@/g" | sed 's/|.*//g'
	done \
	> $wd/all_prots.nm.fa
}


ANI_species_shortlist () {
	echo '
	###################################################
	######## Use Average Nucleotide Identity ##########
	######## To Create Genome shortlist for  ##########
	################### Orthofinder ###################
	###################################################
	'
	date
	func_timing_start

	# set variables
	#   the mygenome_list stuff is to deal with passing array variable to function
	local genome_dir=$1
	local number_shortlist=$2
	local ANI_complete=ANI_complete
	local mygenome_list=($(ls $genome_dir))
	
	echo "ANI clustering is being done on sequences in ${genome_dir}"
	echo "Representative genomes from ${number_shortlist} clusters will be used to generate HMM profiles for each orthogroup"
	cd $store || exit 1
	mkdir ANI_working_dir
	cd ANI_working_dir

	ANI_out=ANI_out

	mkdir $prots.shortlist/
	# use fastANI directly for a list vs list comparison
	# it might do all pairwise comparisons (self and other half of matrix)
	# but testing showed it was way faster then doing it single threaded
	# and making my script parallel would be a huge pain in the arse
	if [ -f genome_names ]
	then
		rm genome_names
	fi
       	for genome1 in ${mygenome_list[@]}
        do
          	echo "$genome1" >> genome_names
        done
	cd $genome_dir || exit
	if [ -f ../ANI_working_dir/$ANI_complete ]
	then
		echo "fastANI already run"
		echo "Skipping ahead"
	else
		echo "NOTE: All fastANI output sent to file (too verbose)." 
		echo "If you think an error occured with FastANI, please check ANI_working_dir/fastANI.stdout"
		$FASTANI \
			--rl ../ANI_working_dir/genome_names \
			--ql ../ANI_working_dir/genome_names \
			-o ../ANI_working_dir/$ANI_out --matrix \
			-t $threads > ../ANI_working_dir/fastANI.stdout 2>&1 \
			&& touch ../ANI_working_dir/$ANI_complete \
			|| exit
	fi
       	cd ../ANI_working_dir || exit
	# makes a rough NJ tree to find the X most distantly related representetives
	# writes to Species_shortlist
	python $ANI_genome_picking $ANI_out genome_names $number_shortlist > ANI_picking.stdout
	for I in $(cat Species_shortlist)
	do
		cp $prots/${I%.*}.faa $prots.shortlist/
	done
}

MASH_species_shortlist () {
	echo '
	#############################################
	######## Use MASH To Create Genome ##########
	######## Shortlist for OrthoFinder ##########
	#############################################
	'
	date
	func_timing_start

	# set variables
	#   the mygenome_list stuff is to deal with passing array variable to function
	local genome_dir=$1
	local number_shortlist=$2
	local threads=$3
	local prots=$4
	local ANI_genome_picking=$5

	local MASH_complete="MASH_complete"
	local mygenome_list=($(ls $genome_dir))
	
	echo "MASH clustering is being done on sequences in ${genome_dir}"
	echo "Representative genomes from ${number_shortlist} clusters will be used to generate HMM profiles for each orthogroup"
	
	cd $store || exit 1
	MASH_working_dir=$store/MASH_working_dir
	prots_shortlist=${prots}.shortlist/
	
	mkdir $MASH_working_dir
	cd $MASH_working_dir

	MASH_out=$MASH_working_dir/MASH_out
	mkdir $prots_shortlist

	# use MASH directly for an all-vs-all assembly comparison
	if [ -f genome_names ]
	then
		rm genome_names
	fi
       	for genome1 in ${mygenome_list[@]}
        do
          	echo "$genome1" >> genome_names
        done
	cd $genome_dir || exit
	if [ -f $MASH_working_dir/$MASH_complete ]
	then
		echo "MASH already run"
		echo "Skipping ahead"
	else
		echo "NOTE: All MASH output sent to file (too verbose)." 
		echo "If you think an error occured with MASH, please check MASH_working_dir/MASH.stdout"
		mash triangle \
			-k 17 -s 5000 \
			-E -p $threads \
			*.f*a \
			> $MASH_out	&& touch $MASH_working_dir/$MASH_complete \
			|| exit
	fi
       	cd $MASH_working_dir || exit
	# makes a rough NJ tree to find the X most distantly related representetives
	# writes to Species_shortlist
	cat $MASH_out | awk '{$3=100*(1-$3) ; print $0 }' | sed 's/ /\t/g' > $MASH_out.similarity
	python $ANI_genome_picking $MASH_out.similarity genome_names $number_shortlist > MASH_picking.stdout
	for I in $(cat Species_shortlist)
	do
		cp $prots/${I%.*}.faa $prots_shortlist/
	done
}

ORTHO_RUN () {
	echo '
	###################################################
	############# Run Orthofinder For #################
	############### Predicted Prots ###################
	###################################################
	'
	date
	func_timing_start
	local_prots_fixed=$1
	cd $store || exit
	if [ -f $store/orthofinder.complete ]
	then
		echo "Orthofinder previously completed."
		echo "if rerun is desired, delete $store/ortho.complete"
		export orthodir=$(echo ${local_prots_fixed}/OrthoFinder/Results_${ortho_trial})
	else
		# this avoids getting a ${orthodir}_X for 1 restart...
		export orthodir=$(echo ${local_prots_fixed}/OrthoFinder/Results_${ortho_trial})
		if [ -d $orthodir ]
		then
			mv $orthodir $orthodir.bk
		fi
		cd $store || exit
		orthofinder -y -oa -t $threads -S diamond -n $ortho_trial \
		-M msa -A mafft -X \
		-f "$local_prots_fixed" > $wd/logs/orthofinder.stdout \
		&& touch orthofinder.complete \
		|| { echo "Orthofinder failed, check $wd/logs/orthofinder.stdout for details"; exit; }
		# make a symlink to OF dirƒ
		cd $wd
		ln -rs $orthodir ./
	fi
}

ANI_ORTHOFINDER_TO_ALL_SEQS () {
	echo '
	###################################################
	#######  Find representetives Of Orthogroups ######
	##### In annotated protiens from All Genomes  #####
	###################################################
	'
	date
	func_timing_start
	local gene_counts=$1
	local OF_alignments=$2
	local all_prots=$3
	local OG_alignmentsToHMM=$4

	local alignments_TMP="$OG_alignmentsToHMM/alignments/"

	# run this block ass a function to use multiple cores and local variables
	OG_hmm_search () {
		local local_wd=$1
		local aligned_dir=$2
		local all_prots=$3
		cd $local_wd || exit
		# build HMM model and run search
		hmmbuild hmms/${I}.hmm ../alignments/${I}.fa > hmmbuild_out.tmp
		# search all annotated prots from $I hmm.
		hmmsearch -T 25 \
		-o ../ORTHOFINDER_TO_ALL_SEQS.out \
		--tblout hmmout/${I}.hmmout \
		hmms/${I}.hmm $all_prots
		##################################################
		# This block maximizes hmm score filter, to eliminate paraogs
		# while keeping all the orthos possible
		# PROBLEM!?! if the model is biased to one clade, the scores will be too
		#   This means that a distantly related ortholog will be arteficially thrown out
		#   However, if a paralog fits this restricted  model very well, it is likey a "true" paralog in that clade
		#	Conclusion? - it is a self limiting problem?
		local no_para_score=$(cat hmmout/${I}.hmmout | grep -ve "^#" | sed 's/@/ /g' | awk '{print $1,$7}' |\
			sort -k1,1 -k2,2rn |\
			awk '{if ($1==prev) print $1,$2} {prev=$1}' |\
			sort -rnk2 |\
			head -n 1 | \
			awk '{print $2}')
		cat hmmout/${I}.hmmout |\
			grep -ve "^#" | \
			awk -v no_para_score=$no_para_score '{if ($6>no_para_score) print $1,$6}' |\
			sort -k1 \
			> hmmout/${I}.list_filter
		# get number of putative paralogs per taxa
		cat hmmout/${I}.hmmout | \
			sed 's/@/\t/g' | awk '{print $1}' | \
			sort | uniq -c | sed 's/^ *//g' | sort -nk1,1 \
			> ${I}.paralogs
		#under construction
		#echo "${I} ${seqs} ${score}" >> OG_fullset_scores
		###################################################
		# pull names from hmmout filtering and make multifasta with filterbyname
		cat hmmout/${I}.list_filter | awk '{print $1}' > hmmout/${I}.names
		filterbyname.sh -Xmx60m -Xms60m overwrite=True include=True ignorejunk=True \
		names=hmmout/${I}.names \
		in=$wd/all_prots.nm.fa \
		out=SequencesProts/${I}.faa >> filterbyname.so_verbose.out 2>&1
		mafft --quiet SequencesProts/${I}.faa > \
            AlignmentsProts/${I}.fa
	}


	dirtbag_multithreading () {
	### Run $threads number of jobs at a time
	###   Waits for all $threads jobs to finish, then starts a new round
	###   some squishy test made this seem faster. Perhaps because file I/O was the bottleneck?
	local local_wd=$1
	local input_list=$local_wd/$2
	local alignments_TMP=$3
	local outdir=$4
	local all_prots=$5
	local threads=$6
    local percent=$(( $(cat $input_list | wc -l) / 10))
	echo "Expanding OrthoFinder OGs to full genome set if found in at least $ANI_shortlist_min_OGs taxa" | tee $run_notes
	cd $outdir || exit

	J=0
	K=0
		for I in $(cat $input_list)
		do
			if test "$(jobs | wc -l)" -ge $threads
			then
			        wait -n
					trap control_c INT
			fi
			
			# Progress sent to stdout
			if [ $(((J + K) % percent)) -eq 0 ]
			then
				echo $(((J+K)/percent*10))" percent of the way through the hmm search"
			fi
			# run the OG_hmm_search module for each alignments found in $alignments_TMP
			#  if number of taxa in alignment is gt $ANI_shortlist_min_OGs
			# outer if statement avoids errors from reps>1
			#  where some cat $input_list are not in alignments because of filtering
			if [ -f $alignments_TMP/${I}.fa ]
			then
				if [ $(cat $alignments_TMP/${I}.fa | grep -e ">" | sort | uniq | wc -l) -ge $ANI_shortlist_min_OGs ]
				then
					OG_hmm_search $outdir $alignments_TMP $all_prots &
					# count actual jobs runnin gin background
					J=$((J+1))
				else
					# count inputs that are skipped (for % complete counter)
					K=$((K+1))
				fi
			else
				# count inputs that are skipped (for % complete counter)
				K=$((K+1))
			fi
		done
		wait

		# copy all HMM search aligned proteins to the alignments_TMP
		for I in $(cat $input_list)
		do
			# move file
			if [ -f $outdir/AlignmentsProts/${I}.fa ]
            then
            	rm $alignments_TMP/${I}.fa
			    cp $outdir/AlignmentsProts/${I}.fa $alignments_TMP/
            fi
			if [ -f $alignments_TMP/${I}.fa ]
			then
				echo -e "Rep_${rep}\t${I}\t"$(cat $alignments_TMP/${I}.fa| grep -c -e ">") \
					>> ${local_wd}/hmm_reps.tsv
			else
				echo -e "Rep_${rep}\t${I}\t0" \
                    >> ${local_wd}/hmm_reps.tsv
			fi
		done
	}

	echo $orthodir
	cd $orthodir || exit
	reps=2
	if [ -f $store/ANI_ORTHOFINDER_TO_ALL_SEQS.complete ]
        then
            	echo "ANI_ORTHOFINDER_TO_ALL_SEQS previously completed."
                echo "if rerun is desired, delete $store/ANI_ORTHOFINDER_TO_ALL_SEQS.complete"

        else
		mkdir $OG_alignmentsToHMM
		cd $OG_alignmentsToHMM || exit
		mkdir alignments
		touch OG_fullset_scores
		# Pick OGs of interest from Orthofinder genecounts
		#    Writes file SCO_$ANI_shortlist_min_OGs which has one OG ID per line
		echo "#### Filtering OGs for those with seqs from at least $ANI_shortlist_min_OGs taxa with no paralogs ####"
		python $OG_sco_filter $gene_counts $ANI_shortlist_min_OGs
		### Identify orthologs in the full data set from hmms derived
		###    from SCO_$ANI_shortlist_min_OGs
		### making progress indicator
		num_OGs=$(cat SCO_$ANI_shortlist_min_OGs | wc -l)
		percent=$(( num_OGs / 10))
		### Run HMM search with updataed models form full dataset
		###
		cd $OG_alignmentsToHMM
		percent=$(( num_OGs / 10))
		J=0
		for I in $(cat SCO_$ANI_shortlist_min_OGs)
		do
			cp $OF_alignments/${I}.fa $alignments_TMP/
			echo -e "Rep_0\t${I}\t"$(cat $alignments_TMP/${I}.fa| grep -c -e ">") \
				>> $OG_alignmentsToHMM/hmm_reps.tsv
		done

		for rep in $(seq $reps)
		do
			echo "##############  LIST OF ALIGNMENTS TO START Rep $rep  ############"
			cd $OG_alignmentsToHMM
			mkdir hmm_round$rep
			cd hmm_round$rep
			mkdir alignments
               		mkdir hmms
                	mkdir prots
                	mkdir hmmout
                	mkdir SequencesProts
                	mkdir AlignmentsProts
			echo "Running OG expantion rep ${rep}"
			dirtbag_multithreading $OG_alignmentsToHMM \
				SCO_$ANI_shortlist_min_OGs \
				$alignments_TMP \
				$OG_alignmentsToHMM/hmm_round$rep \
				$all_prots \
				$threads
		done

		#alignments copied to the final module output ($wd/AlignmentsProts/)
		for I in $(cat $OG_alignmentsToHMM/SCO_$ANI_shortlist_min_OGs)
		do
			if [ -f $alignments_TMP/${I}.fa ]
			then
				cp $alignments_TMP/${I}.fa $wd/AlignmentsProts/${I}.faa
			fi
		done
	fi \
	&& touch $store/ANI_ORTHOFINDER_TO_ALL_SEQS.complete
}

REALIGN_ORTHOGROUP_PROTS () {
	echo '
	###################################################
	############# Realign Orthogroup Prots ############
	#################### With MAFFT ###################
	###################################################
	'
	date
	func_timing_start
	#realigning prots because Orthofinder does some trimming which messes up prot 2 trans alignments
	#..could turn that behavior off in OF
	#..but already wrote this bit and it allows customization down the line
	#..also trimming durring Orthofinder run might help with OG accuracy?
	cd $wd || exit
	# making progress indicator
	num_OGs=$(ls "$orthodir"/MultipleSequenceAlignments/OG*.fa | wc -l)
	percent=$(( num_OGs / 10))
	J=0
	for i in "$orthodir"/MultipleSequenceAlignments/OG*.fa
	do
		if test "$(jobs | wc -l)" -ge $threads
		then
			wait -n
			trap control_c INT
		fi
		# Progress sent to stdout
		if [ $((J % percent)) -eq 0 ]
		then
			echo $((J/percent*10))" percent of the way through the extraction and realignment of OG prots"
		fi
		filter_and_Align_subfunc () {
			local base=$(basename "${i%.*}")
		    cat $i | grep ">" | sed 's/>//g' | sed 's/|.*//g' \
		        > ./OG_names/${base}.names
			# pulling out prot sequences based on names in orthofinder OGs
		    filterbyname.sh -Xmx60m -Xms60m include=t \
		        names=./OG_names/${base}.names ignorejunk=t \
		        in=$wd/all_prots.nm.fa out=./SequencesProts/${base}.faa \
				>> $wd/logs/filterbyname.realign_prots 2>&1
			# realign prots with mafft
			mafft --quiet $wd/SequencesProts/${base}.faa > \
                	$wd/AlignmentsProts/${base}.faa
		}
		# run the above subfunction multithreaded (sort-of)
		filter_and_Align_subfunc &
		J=$((J+1))
		
	done
	wait
}

GET_OG_NAMES () {
	echo '
	##############################################
	############# Parsing Names From  ############
	############## Orthogroup Prots ##############
	##############################################
	'
	OG_names_dir=$1
	OG_dir=$2
	
	local num_OGs=$(ls "$OG_dir"/OG*.f*a | wc -l)
	local percent=$(( num_OGs / 10))
	J=0
	for i in "$OG_dir"/OG*.f*a
	do
		if test "$(jobs | wc -l)" -ge $threads
		then
			wait -n
			trap control_c INT
		fi
		Get_OG_Names () {
			local base=$(basename "${i%.*}")
			#echo $base
		    cat $i | grep ">" | sed 's/>//g' | sed 's/|.*//g' \
		        > $OG_names_dir/${base}.names
		}
		Get_OG_Names &
		# Progress sent to stdout
		if [ $((J % percent)) -eq 0 ]
		then
			echo $((J/percent*10))" percent of the way through the extracting OG_names from HMM OG hits"
		fi
		J=$((J+1))
	done 
	wait
	echo "Get_OG_Names 100% done!!!" 
}




filterbyname_subfunc () {
	echo '
	#########################
	### Pull seqs for OGs ###
	#########################
	'
	names_dir=$1
	all_seqs=$2
	out_dir=$3
	J=0
	num_OGs=$(ls ${names_dir}/OG* | wc -l)
	percent=$(( num_OGs / 10))
	for i in $(ls $names_dir)
	do
	file=${i##*/}
	base=${file%%.*}
	# pulling out CDS sequences based on names in orthofinder OGs
	filterbyname.sh -Xmx60m -Xms60m include=t \
		names=$i ignorejunk=t \
		in=$all_seqs out=$out_dir/${base}.fa \
		>> $logs/filterbyname_CDS 2>&1
	J=$((J+1))
	done

}

ALIGN_PROT_n_CDS () {
	echo '
	###################################################
	############# Realign Orthogroup Prots ############
	############### And Nucs With MACSE ###############
	###################################################
	'
	# 
	date
	func_timing_start
	local wd=$1
	local OG_names=$2
	local all_CDS=$3
	local SequencesCDS=$4
	local code=$5
	local AlignmentsCDS=$6
	local AlignmentsProts=$7
	local threads=$8
	local logs=$9
	#macse_parameters="-optim 1 -max_refine_iter 3 -local_realign_init 0.3 -local_realign_dec 0.2 "
	macse_parameters="-optim 0"
	cd $wd || exit
	# making progress indicator
	num_OGs=$(ls "$OG_names"/OG*.names | wc -l)
	percent=$(( num_OGs / 10))
	J=0
	for i in "$OG_names"/OG*.names
	do
		if test "$(jobs | wc -l)" -ge $threads
		then
			wait -n
			trap control_c INT
		fi
		
		base=$(basename ${i%.*})

		# !!! redundant with new filterbyname_subfunc
		filter_and_Align_subfunc () {
			# pulling out CDS sequences based on names in orthofinder OGs
			filterbyname.sh -Xmx60m -Xms60m include=t \
		        names=$i ignorejunk=t \
		        in=$all_CDS out=$SequencesCDS/${base}.fa \
				>> $logs/filterbyname.realign_CDS 2>&1
			# realign prots with MACSE
			macse -prog alignSequences -seq $SequencesCDS/${base}.fa ${macse_parameters} >> $logs/macse.out 2>&1 
			mv $SequencesCDS/${base}_NT.fa $AlignmentsCDS/${base}.fa  || (echo "NT Alignment file not created for $base" && exit 1)
			mv $SequencesCDS/${base}_AA.fa $AlignmentsProts/${base}.faa || (echo "AA Alignment files not created for $base" && exit 1)
			
			#### TEST FIX FOR ReLeaf partition file not match prot length (* get converted to - by mafft...ugh)
			# might move to prot specific part (not needed for only CDS)
			sed -i 's/*/-/g' $AlignmentsProts/${base}.faa

		}
		# run the above subfunction multithreaded (sort-of)
		filter_and_Align_subfunc &
		J=$((J+1))
		# Progress sent to stdout
		if [ $((J % percent)) -eq 0 ]
		then
			echo $((J/percent*10-10))" percent of the way through the extraction and realignment of OG prots"
		fi
	done
	echo "Waiting for the last " $(jobs | wc -l) " alignments to finish..."
	wait
}

TRIM () {
	echo "
	##############################################
	#### trim all ${2} alignments with Trimal ####
	####### and remove gene specific names #######
	##############################################
	"
	date
	func_timing_start
	local alignment_dir=${1}
	local alignment_type=${2}
	local checkpoint_file=${3}

	if [ -f $checkpoint_file ]
    then
            	echo "Protein alignment trimming already completed."
                echo "if rerun is desired, delete " $checkpoint_file

    else
		echo -e "\nTrimming $alignment_type alignments in $alignment_dir with: 'trimal $trimal_parameter'\n"

		# make output dir
		mkdir $alignment_dir.trm

		# set up percent complete updater 
		cd $wd/ || exit
		num_OGs=$(ls $alignment_dir/OG* | wc -l) # this is 1+ the real num
			percent=$(( num_OGs / 10))
			J=0
		
		# iterate over alignment dir OG files works in chunks of %thread 
		#   Need to change it to start next job if there are < $threads trimming operations running
		for I in $alignment_dir/OG*
		do
			if test "$(jobs | wc -l)" -ge $threads
			then
				wait -n
				trap control_c INT
			fi
			# Progress sent to stdout
			if [ $((J % percent)) -eq 0 ]
			then
					echo $((J/percent*10))" percent of the way through trimming"
			fi
			TRIM_subfunc () {
				base=$(basename ${I%.*})
				trimal \
					-in $I -fasta \
					-out $alignment_dir.trm/${base}.trm.fa \
					$trimal_parameter \
					-colnumbering 1> trimmed_columns/${base}.cols 2> $wd/logs/trimal_log
				#remove all gene info from header for concat
				#this is really janky: relies on adding the @ during renaming...
				# Removes everything after the @ to leave just the sample name
				#cat $alignment_dir.trm/${base}.trm.fa \
				#| sed 's/@.*$//g' > $alignment_dir.trm.nm/${base}.trm.nm.fa
				# moved to TRIMAL_backtrans
			}
			TRIM_subfunc &
			J=$((J+1))
		done
		if test "$(jobs | wc -l)" -gt 0
			then
				wait -n
				trap control_c INT
			fi
	fi
}

GET_TRIMMED_COLS () {
    echo "
    ############################################
    Getting trimmed columns from original OP Run
    ############################################
    "
			# requires a directory full of trimal output for the -colnumbering argument (ending in .cols) 
			cols_dir=$1
            cd $cols_dir
            for I in $(ls OG*.cols); do
                base=${I%%.*}
				#get protein colums to keep through trimming
                trimal_prot_cols=$base.trimal_prot_cols
                cols_line=$(cat $I | sed 's/#ColumnsMap/{ /g;s/\t//g;s/$/ }/g;s/, /,/g')
                echo -n "-selectcols" $cols_line " -complementary "> $trimal_prot_cols
				# get CDS columns to keep based on prot columns  
				trimal_CDS_cols=$base.trimal_CDS_cols
				cols_line=$(cat $I | sed 's/#ColumnsMap\t//g;s/ //g' |\
					tr ',' '\n' |\
					awk '{print $0*3","$0*3+1","$0*3+2}' |\
					tr '\n' ',' |\
					sed 's/,$//g')
				echo -n "-selectcols { " $cols_line " } -complementary "> $trimal_CDS_cols
			done    
}

TRIM_selectcols () {
	echo "
	##############################################
	#### trim all ${2} alignments with Trimal ####
	####### and remove gene specific names #######
	##############################################
	"
	date
	func_timing_start
	local alignment_dir=${1}
	local alignment_type=${2}
	local cols_dir=${3}
	local checkpoint_file=${4}

	if [ -f $checkpoint_file ]
    then
            	echo "CDS codon alignment trimming already completed."
                echo "if rerun is desired, delete " $checkpoint_file

    else
		echo -e "\nTrimming $alignment_type alignments in $alignment_dir with columns in $cols_dir\n"

		# make output dir
		mkdir $alignment_dir.trm

		# set up percent complete updater 
		cd $wd/ || exit
		num_OGs=$(ls $alignment_dir/OG* | wc -l) # this is 1+ the real num
			percent=$(( num_OGs / 10))
			J=0
		
		# iterate over alignment dir OG files - works in chunks of #threads 
		#  starts next job if there are < $threads trimming operations running
		for I in $alignment_dir/OG*
		do
			base_tmp=${I##*/}
			base=${base_tmp%%.*}
			if test "$(jobs | wc -l)" -ge $threads
			then
				wait -n
				trap control_c INT
			fi
			# Progress sent to stdout
			if [ $((J % percent)) -eq 0 ]
			then
					echo $((J/percent*10))" percent of the way through trimming"
			fi
			trimal_parameter=$(cat $cols_dir/${base}.trimal_${alignment_type}_cols)
			TRIM_subfunc () {
				trimal \
					-in $I -fasta \
					-out $alignment_dir.trm/${base}.codon_aln.trm.fa \
					$trimal_parameter \
					> $wd/logs/trimal_log
				
				
			#remove all gene info from header for concat
			#if	[[ $alignment_type == "CDS" ]] 
			#then
				#this is really janky: relies on adding the @ during renaming...
				# Removes everything after the @ to leave just the sample name
				#	Moved to its own func
			#	cat $alignment_dir.trm/${base}.trm.fa \
			#	| sed 's/@.*$//g' > $alignment_dir.trm.nm/${base}.trm.nm.fa
			#fi
			}
			TRIM_subfunc &
			J=$((J+1))
		done
		wait
	fi
}

ALIGNMENT_STATS () {
	echo '
       	####################################################
       	########## Run alignment assessment 4 OGs  #########
       	####### from the provided directory (below) ########
       	####################################################
        '
	date
	func_timing_start
	alignment_dir=$1
    alignment_dir_vis=$(basename $1).vis
	echo "making alignemnt assesment figures for $alignment_dir"
	cd $alignment_dir/ || exit

	# make phylip formatted alignments for Alignment_Assessment
	for I in *.fa
	do
		cat $I | awk -vRS=">" -vFS="\n" -vOFS="" \
			'$0!=""{$1=substr($1,1,15);$1=sprintf ("%-17s",$1)}$0!=""' \
			> TMP.phy1
		if [ $(cat TMP.phy1 | wc -l) -gt 1 ]
		then
			num=$(cat TMP.phy1 | wc -l)
			len=$(cat TMP.phy1 | head -n 1 | awk '{print length($2)}')
			echo -e "\t"$num"\t"$len > ${I%.*}.phy
			cat TMP.phy1 >> ${I%.*}.phy
		fi
		rm TMP.phy1
	done

	# generate Master_Alignment_Assessment.txt with dportik's Alignment_Assessment GH repo
	echo "running Alignment_Assessment_v2.py to get general per OG stats"
	python $Alignment_Assessment \
		./ \
		> TMP.alignment_assessment
	cd ../
	mkdir $alignment_dir_vis
	cd $alignment_dir_vis || exit
	cp $alignment_dir/Alignment_Assessment/Master_Alignment_Assessment.txt ./

	#get number of genes per taxa
	echo "finding the number of OG genes found per taxa"
	cat $alignment_dir/*.fa | grep -e ">" | sort | uniq -c | sed 's/^ *//g;s/>//g' | sort -rnk1,1 \
		> num_genes_per_taxa.tsv

	# run modified R script to auto-generate pdf figures
	Rscript $script_home/Rscripts/Alignment_Assessment_vis.R \
		Master_Alignment_Assessment.txt \
		num_genes_per_taxa.tsv > verbose_log.txt
}

SCO_MIN_ALIGN () {
	echo "
	###################################################
	########## Identify single copy orthologs #########
	###### with at least ${2} taxon representation ####
	###################################################
	"
	date
	func_timing_start
	# takes folder of fasta alignments and will pull alingment names with >= $min_frac_orthos constituents
	local alignment_dir=${1}
	local min_frac_orthos=${2}
	local OG_sco_filter=$OG_sco_filter

	if [[ $min_frac_orthos == $(cat $store/all_input_list | wc -l) ]]
	then
		outstring="strict"
	else
		outstring=$min_frac_orthos
	fi

	if [ "$ANI" = true ]
	then
		echo "Finding all ANI OGs with $min_frac_orthos representetives in $alignment_dir"
		# make list of SCO $min_frac_orthos directly from fasta
		cd $wd || exit
		# Iterate over alignments and pull out the number of samples per multifasta
		# this has the assumption that there is only one sequence per genome
		# should be enforced by ORTHOFINDER_TO_ALL_SEQS
		# could add an error || exit if it becomes a problem
		for I in $alignment_dir/OG0*
		do
			# test transcript alignment for number of homologs (should have no paralogs)
			# if greater than min_frac_orthos (samples*min_frac_ortho) 
			# add to SCO_$min_frac_orthos list
			final_seqs=$(cat $I | grep -e ">" | wc -l)
			if [ $final_seqs -ge $min_frac_orthos ]
			then
				# this basename call  is too specific and easy to break...should clean up
				base=$(basename ${I%.*.*})
				echo $base >> $wd/SCO_$outstring
			fi
		done
	else
		# going to ignore this for now...
		cd $orthodir/ || exit
		gene_counts="$orthodir/Orthogroups/Orthogroups.GeneCount.tsv"
		echo "Finding all OGs directly from OrthoFinder with $min_frac_orthos representetives"
		# finds OGs with at least $min_orthologs SCOs
		# and writes to $orthodir/SCO_$min_frac_orthos
		python $OG_sco_filter $gene_counts $min_frac_orthos
		mv SCO_$min_frac_orthos $wd/SCO_$outstring
	fi
}


TRIMAL_backtrans () {
	echo '
	###################################################
	#####  Use TRIMAL_backtrans to produce Nuc alignments ######
	############### from AA alignemtns ################
	###################################################
	#loops of all Orthogroups to
	#	extract CDS sequences from a contatentated CDS file
	#	then do codon alignments with TRIMAL_backtrans
	'
	date
	func_timing_start
	local wd=${1}
	local prot_alignment=${2}
	local CDS_file=${3}
	local CDS_codon_alignment=${4}
	local OG_names=$5
	local SequencesCDS=$6
	local addem=$7

	if [ ! -d $CDS_codon_alignment.nm ]
	then
		mkdir $CDS_codon_alignment.nm
	fi
	
	#aligning transcripts based on orthogroup protien alignments
	cd $wd || exit
	num_OGs=$(ls $prot_alignment/OG*.faa | wc -l) # this is 1+ the real num
	percent=$(( num_OGs / 10))
        J=0
	for i in $prot_alignment/OG*.faa
	do
		if test "$(jobs | wc -l)" -ge 1
		then
			wait -n
			trap control_c INT
		fi
		# Progress sent to stdout
		if [ $((J % percent)) -eq 0 ]
		then
			echo $((J/percent*10))" percent of the way through the TRIMAL_backtrans"
		fi
		# why the hell did I declare this function with in the loop...probably because of the "local base=...". this should be passed in with a variable. 
		TRIMAL_backtrans_subfunc () {
			#get basename (OG000????)
			local file=${i##*/}
			local base=${file%%.*}
			
			#only run if in the addem workflow
			if [ "$addem" = "true" ]
			then
				trimal_parameter="$(cat $OLD_trim_cols/$base.cols_to_keep)"
			fi
			
			cat $i | grep ">" | sed 's/>//g' | sed 's/|.*//g' \
	        	> $OG_names/${base}.names
			#pull out trans sequences for each OG
			filterbyname.sh -Xmx60m -Xms60m include=t \
	        	names=$OG_names/${base}.names \
				in=$CDS_file out=$SequencesCDS/${base}.CDS.fa \
				>> $wd/logs/filterbyname.TRIMAL_backtrans 2>&1 # changed a "/" to ">" not sure how the code ran before...
			# protein to nuc alignments
			#pal2nal.pl $i ./SequencesCDS/${base}.CDS.fa -codontable 11 -output fasta \
			#> $CDS_codon_alignment/${base}.codon_aln.fa 2> TRIMAL_backtrans.stderr.tmp
			
			# backtraslation
			trimal -in $i \
				-ignorestopcodon \
				-backtrans $SequencesCDS/${base}.CDS.fa \
				-out $CDS_codon_alignment/${base}.codon_aln.fa \
				>> $wd/logs/TRIMAL_backtrans 2>&1

		}
		TRIMAL_backtrans_subfunc &
        J=$((J+1))
	done
	wait
}


prot_rename () {
	echo "
	#################################
	### Removing gene specific ###
	#### info from prot fasta ####
	##############################
	"
	local_wd=$1
	mkdir $local_wd.nm
	for I in $local_wd/OG*
	do
		file=${I##*/}
		base=${file%%.*}
		cat $I \
	    	| sed 's/@.*$//g' > $local_wd.nm/${base}.trm.nm.fa
	done
}

CDS_rename () {
	echo "
	#########################################
	######## Removing gene specific #########
	#### info from CodonAlignment fastas ####
	#########################################
	"
	local_wd=$1
	mkdir $local_wd.nm
	for I in $local_wd/OG*
	do
		local file=$(basename $I)
		local base=${file%.*}
		cat $I \
	       	| sed 's/@.*$//g' > $local_wd.nm/${base}.nm.fa 
	done
}

cat_alignments () {
	echo '
		##############################################
		####### Place SCOs in dedicated folder #######
		############ SCO_*.align and make ############
		####### concattenated alignment for ML #######
		##############################################
		'
	date
	func_timing_start
	# function inputs
	local SCO_list=${1}
	local alignment_dir=${2}
	local output_dir=${3}
	local alignment_type=${4}

	if [ ! -d $output_dir ]
	then 
		mkdir $output_dir
	fi

	# output valiable
	local SCO_dir=$SCO_list.$alignment_type.align
	local cat_alignment_file=$output_dir/$(basename $SCO_list.$alignment_type.trm.sco.nm)
	
	cd $wd || exit
	if [ ! -d $SCO_dir ]
	then	
		mkdir $SCO_dir
	fi
	echo "Creating directory with alignments of OGs from $SCO_list"
	for I in $(cat $SCO_list)
	do
		cp $alignment_dir/${I}.* $SCO_dir/
	done
 
	echo "Concatenating fasta alignments to phylip format"
	cd $output_dir || exit
	perl $catfasta2phyml_cmd -c $SCO_dir/*.fa \
		1> $cat_alignment_file.phy 2> $cat_alignment_file.partitions.tmp
	echo "Concatenating fasta alignments to fasta format"
	perl $catfasta2phyml_cmd -f -c $SCO_dir/*.fa \
        1> $cat_alignment_file.fa 2>> $store/verbose_log.txt
	# fixing partion file to make compatable with iqtree
	if [ $alignment_type == "CDS" ]
	then
		type=DNA
	elif [ $alignment_type == "PROT" ]
	then
		type=PROTEIN
	fi
	cat $cat_alignment_file.partitions.tmp \
		| awk -v type=$type '{print type",\t"$0}' \
		> $cat_alignment_file.partitions #&& rm $cat_alignment_file.partitions.tmp
}

#####################################
#####################################



TREE_BUILD () {
	echo '
	##########################################
	###### Build trees and bootstrap #########
	####### from SCO_$min_frac_orthos #########
	########### or SCO_stricts ###############
	##########################################
	'
	date
	func_timing_start
	echo "Building trees from...."
	echo "$2"
	output_dir=$1
	input_alignment=$2
	threads=$3
	alignment_type=$4
	output_name=$(basename ${input_alignment%.*.*.*.*})
	partition_file=${input_alignment%.*}.partitions
	
	# grab tree options based on what type of data (CDS vs PROT)
	if [[ $alignment_type == "CDS" ]]
	then
		# load CDS specific options
		RAxML_speciestree_options=$RAxML_CDS_speciestree_options
		fasttree_speciestree_options=$fasttree_CDS_speciestree_options
		IQtree_speciestree_options=$IQtree_CDS_speciestree_options
		IQtree_partition_options=$IQtree_CDS_partition_options
	
	elif [[ $alignment_type == "PROT" ]]
	then
		# load PROT specific parameters
		RAxML_speciestree_options=$RAxML_PROT_speciestree_options
		fasttree_speciestree_options=$fasttree_PROT_speciestree_options
		IQtree_speciestree_options=$IQtree_PROT_speciestree_options
		IQtree_partition_options=$IQtree_PROT_partition_options
	fi

	# if using partitions set up IQtree and RAxML options
	if [ $use_partitions == "true" ]
	then
		# edge-unlinked partition merging for IQTREE
		IQtree_partitions="-Q $partition_file $IQtree_partition_options"
		RAxML_partitions="-q $partition_file"
	fi

	cd $output_dir || exit
	# subfunction to run RAxml
	RAxML_run () {
		mkdir raxml
		cd raxml
		raxmlHPC-PTHREADS $RAxML_speciestree_options \
		-s $input_alignment \
		$RAxML_partitions \
		-n ${output_name} \
		> ./RAxML_output.1 2>&1
		treefile=$(ls RAxML_bipartitionsBranchLabels*)
		mv $treefile ../${treefile}.tree && touch ./${treefile%.*}.complete || echo -e "\n!!!!!!!!\nWARNING\n!!!!!!!!!!\n\tRAxML failed to run on "$input_alignment
		cd $output_dir || exit
	}

	# subfuncton to run FastTreeMP
	FASTTREE_run () {
		# fastTree doesnt like any of the phylip formats I have tried. using .fa
		#    A bit clunky...
		export OMP_NUM_THREADS=$threads
		FastTreeMP $fasttree_speciestree_options \
		-out ./fastTree.${output_name}.tree_working \
		-nt ${input_alignment%.*}.fa \
		> ./fastTree_log.${output_name} 2>&1
		treefile=fastTree.${output_name}.tree_working
		mv $treefile ./${treefile%.*}.tree && touch ./${treefile%.*}.complete || echo -e "\n!!!!!!!!\nWARNING\n!!!!!!!!!!\n\tFastTree failed to run on "$input_alignment
		cd $output_dir || exit
	}
	IQTREE_run () {
		if [ ! -d "iqtree" ]; then mkdir iqtree ; fi
		cd iqtree
		# deal with numerical overflow on large datasets. 
		if [[ $(cat $store/genome_list | wc -l) -gt 400 ]]
		then
			safe="--safe"
		fi
		iqtree -s $input_alignment \
		--prefix iqtree.${output_name} $IQtree_partitions $IQtree_speciestree_options $safe\
		-T AUTO --threads-max $threads --seed 1234 > iqtree.${output_name}.long_log
		treefile=$(ls iqtree.${output_name}.treefile)
		cp $treefile ../${treefile%.*}.tree && touch ../${treefile%.*}.complete || echo -e "\n!!!!!!!!\nWARNING\n!!!!!!!!!!\n\tIQTree failed to run on "$input_alignment
		cd $output_dir || exit

	}

	# decide which tree method(s) to use for the cancatenate gene nuc matrix
	

	if [[ " ${tree_method[*]} " =~ " raxml " ]]
	then
		echo "################################"
		echo "Building species tree with RAxML"
		echo "################################"
		RAxML_run
	fi
	if [[ " ${tree_method[*]} " =~ " fasttree " ]]
	then
		echo "###################################"
		echo "Building species tree with FastTree"
		echo "###################################"
		func_timing_start
		FASTTREE_run
	fi
	if [[ " ${tree_method[*]} " =~ " iqtree " ]]
	then
		echo "#################################"
		echo "Building species tree with IQtree"
		echo "#################################"
		func_timing_start
		IQTREE_run
	fi
	
}

orthofinderGENE2SPECIES_TREE () {
	echo '
	#############################################################
	Species tree reconstruction from all Orthogroup protien trees
	#############################################################
	Running Astral on all orthologs from orthofinder
	'

	# This currently doesnt work because Astal P, which is required for dealing with paralogs
	#   is incompatable with our cluster
	#  Jarfile needs to be rebuild...low priority
	mkdir $wd/astral_trees
	cd $wd/astral_trees || exit
	#make names compatable with Astral and concatenate
	# remove everything but species/strain identifyer with
	# non-greedy sed (name@blah.blah.blah:BranchLength) > name:Branchlength
	cat $orthodir/Gene_Trees/OG000*_tree.txt | \
	sed 's/;/;\n/g' | sed 's/@[^:]*:/:/g' \
	> all_genes.orthofinder.tre
	java -Djava.library.path=$ASTRAL_P_lib -jar $ASTRAL_P  -i all_genes.orthofinder.tre \
	-o all_genes.orthofinder.astral.tre \
	-t $threads
}

allGENE_TREEs () {
	echo '
	##################################################
	####### Build Gene trees for CDS or PROT  ########
	########## Alignments to use for Astral ##########
	############# Species Tree Estemation ############
	##################################################
	'
	date
	func_timing_start
	local alignment_type=$1
	local gene_alignment_dir=$2
	local out_dir=$3
	local logs=$wd/logs/gene_tree_logs

	echo "Building gene trees with "$alignment_type" alignments using "${gene_tree_methods[@]}
	if [ -d $out_dir/ ]
	then
	   rm -r $out_dir.bk
	   mv $out_dir/ $out_dir.bk
	fi
	mkdir $out_dir/
	mkdir $out_dir.nm/
	cd $out_dir || exit
	num_OGs=$(ls $gene_alignment_dir/OG*.fa | wc -l) # this is 1+ the real num
        percent=$(( num_OGs / 10))
        J=0

	#set some tree building parameters for CDS vs Prot
	if [[ $alignment_type == "CDS" ]]
	then
		local fasttree_options=$fasttree_CDS_genetree_options
		local iqtree_options=$iqtree_CDS_genetree_options
	elif [[ $alignment_type == "PROT" ]]
	then
		local fasttree_options=$fasttree_PROT_genetree_options
		local iqtree_options=$iqtree_PROT_genetree_options
	fi

	# loops through all alignments
	num_OG_lt_4_OGs=0
	for i in $gene_alignment_dir/OG*.fa
	do
		if test "$(jobs | wc -l)" -ge $threads
		then
			wait -n
			trap control_c INT
		fi
		# Progress sent to stdout
		if [ $((J % percent)) -eq 0 ]
		then
				echo $((J/percent*10))" percent of the way through building gene trees"
		fi
		TRANS_TREE_subfunc () {
			local file=${i##*/}
			local base=${file%%.*}
			# filter out sequences if they have lt 4 sequences. Not informative. 
			#  Still throws an error if there are 3 seqs that are identical
			if [[ $(cat $gene_alignment_dir/${base}.*.fa | grep ">" | wc -l) -lt 4 ]]
			then 
				echo -e "$gene_alignment_dir/${base}.*.fa has less than 4 samples in it...Skipping" >> $logs
				export num_OG_lt_4_OGs=$((num_OG_lt_4_OGs+1))
			else
				if [[ " ${gene_tree_methods[*]} " =~ " fasttree " ]]
				then
					fasttree $fasttree_options \
					$gene_alignment_dir/${base}.*.fa \
					> ./${base}.fasttree.tree 
					cat ./${base}.fasttree.tree | sed 's/@[^:]*:/:/g' \
					> $out_dir.nm/${base}.fasttree.tree
				fi
				if [[ " ${gene_tree_methods[*]} " =~ " iqtree " ]]
				then
					iqtree $iqtree_options \
					-s $gene_alignment_dir/${base}.*.fa \
					--prefix ${base} -T 1 > ${base}.iqtree.log
					mv ./${base}.treefile ./${base}.iqtree.tree
					cat ./${base}.iqtree.tree | sed 's/@[^:]*:/:/g' \
					> $out_dir.nm/${base}.iqtree.tree
				fi
			fi
		}

		TRANS_TREE_subfunc &
		J=$((J+1))
	done
	wait
	echo -e "\n$num_OG_lt_4_OGs Ortholog alignments were not used to build gene trees because they had less than 4 sequences\n"
}



astral_allCDSGENE2SPECIES_TREE () {
	echo '
	#######################################################
	Species tree reconstruction from CDS or PROT alignments
	#######################################################
	'

	#astral all genes
	#TODO
	echo "Does nothing!
	TODO
	Need to set it up to run astral_p to handle paralogs
	I think this will only be useful in limited cases:
	back burner
	ALSO, functionality could be put into astral_CDSGENE2SPECIES_TREE
	just need to change astral_CMD
	"
}


astral_GENE2SPECIES_TREE () {
	echo '
	#######################################################
	Species tree reconstruction from CDS or PROT gene trees 
	#######################################################
	'
	date
	func_timing_start
	local alignment_type=$1
	local in_dir=$2
	local genelist=$3

	echo -e "\nRunning Astral on "${genelist}" genelist for "$alignment_type" using "${gene_tree_methods[@]}" gene trees.\n"

	local base_list=$(basename ${genelist})
	local ASTRAL_cmd=$ASTRAL_cmd
	#outputs
	local out_dir=astral_trees
	local excluded=$wd/$out_dir/trimmed_from_${base_list}_gene_list
	local concat_trees=$wd/$out_dir/${base_list}

	cd $wd || exit
	mkdir $out_dir
	# use iqtree or fasttree (or both)
	for method in ${gene_tree_methods[@]}
	do
		concat_tree_method=$concat_trees.${method}.trees
		astral_tree=astral.${base_list}.${alignment_type}.${method}.tree
		#remove concatenated tree file if it exists.
		if [ -f $concat_tree_method ]
		then
		   rm -r $concat_tree_method.bk
		   mv $concat_tree_method concat_tree_method.bk
		fi

		# Delete previous file of genes ommited because trimming reduced the number of leaves to > 4
		# (in case of rerun)
		if [ -f $excluded ]
		then 
			rm $excluded
		fi

		# concatenate gene trees in $in_dir with names in $genelist
		for I in $(cat $genelist)
		do
            base=$(basename ${I%.*})
			OG_tree=$in_dir/${base}.${method}.tree
			if [ -f $OG_tree ]
			then
				# only include trees with > 3 leaves
				# some genetrees are reduced from min_num_align because of trimming... I think...
				if [ $(awk -F',' '{print NF-1}' $OG_tree) -gt 3 ]
				then
					cat $OG_tree \
					>> $concat_tree_method
				fi
			else
				# add ${base} to list of removed genes (because of trimming)
				echo ${base} >> $excluded
			fi
		done
		cd $wd/$out_dir || exit
		java -jar $ASTRAL_cmd \
		-i $concat_tree_method \
		-o $astral_tree 2>&1 >> $astral_tree.log
	done
}


WRAP_UP () {
        echo '
       	########################################################
        ########  Cleaning up some intermediate files   ########
       	########################################################
        '
	date
	func_timing_start
	# clean up files
	if [[ $cleanup == "TRUE" ]]
	then
		# remove files and dirs. $clean_me declared in control files
		mkdir $store/tmp_trash || exit
		for I in ${clean_me[*]}
		do
			if [ -f $I ]
			then
				mv $store/$I $store/tmp_trash
			fi
		done
		rm -r $store/tmp_trash
	fi
	# but keep AlignmentCDS.trm.nm ans AlignemntProt
	# gather all SpeciesTrees in one place?
	echo '
	#######################################################
	########## Moving final species trees to  #############
	###########  $store/FINAL_SPECIES_TREES ###############
	#####   And running a Generalized RF comparison  ######
	#######################################################
	'

	mkdir $store/FINAL_SPECIES_TREES
	cp $wd/SpeciesTree/*.tree $store/FINAL_SPECIES_TREES
	cp $wd/astral_trees/*.tree $store/FINAL_SPECIES_TREES
	ete3 compare --unrooted --taboutput \
		-t $store/FINAL_SPECIES_TREES/*tree -r $store/FINAL_SPECIES_TREES/*tree \
		>> $store/FINAL_SPECIES_TREES/trees.compare
}

TESTER_compare () {
	echo "
	#############################################
	######   Compare Trees generated with  ######
    #######  TESTER to prepackaged Trees  #######
	#############################################
	"
	# TODO:
	## generalize for full pipeline
	cd $store/FINAL_SPECIES_TREES || exit
	if [ -f $script_home/TESTER.fail ]
	then
		rm TESTER.fail
	fi
	for I in ./*tree
	do
		ete3 compare --unrooted --taboutput \
		-t $I -r $script_home/TESTER/REFERENCE_TESTER_TREES/$I \
		> ${I%.*}.compare
		RF_test=$(cat ${I%.*}.compare | tail -n 1 | awk -F "\t" '{print $5}')
		if (( $(echo "$RF_test == 0" | bc -l) ))
		then
			echo $I " was generated correctly"
		else
			echo $I " did not get generated correctly"
			touch $script_home/TESTER.fail
		fi
		if [ ! -f $script_home/TESTER.fail ]
		then
			echo '
			###### It looks like the install was successful ######
			'
			touch $script_home/TESTER.pass
		else
			echo '
               	        ###### Unsuccessful install######
                       	'
		fi
	done
}
