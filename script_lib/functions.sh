#!/bin/bash

# functions used in orthophylo.sh
# i.e.
# source $script_home/scipt_lib/functions.sh

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
		echo "   Resetting min_num_orthos in case OrthoPhyl crashed between Setup and SCO_min_num_orthos"
		export min_num_orthos=$(printf %.0f $(echo "$(cat $store/all_input_list | wc -l)*$min_frac_orthos" | bc))
		echo " If this is not correct please rm $store/setup.complete and restart OrthoPhyl"
	else
		# set up annotation output directories
		mkdir $store
		mkdir $trans
		mkdir $prots
		mkdir $annots

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
		for I in $(ls ./ | grep \)); do   mv $I ${I//\)/}; done
		for I in $(ls ./ | grep \(); do   mv $I ${I//\(/}; done
		for I in $(ls ./ | grep "_\=_" ); do mv $I ${I//_\=_} ; done

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
		export min_num_orthos=$(printf %.0f $(echo "$(cat $store/all_input_list | wc -l)*$min_frac_orthos" | bc))
		echo "All Single copy orthologs and SCOs represented in $min_num_orthos will be used to create species trees with ASTRAL and ML method chosen" | tee -a $run_notes
		###########################################
		# just setting up the directory structure #
		###########################################
		mkdir $wd
		cd $wd || exit
		mkdir AlignmentsProts/
		mkdir AlignmentsTrans/
		mkdir AlignmentsTrans.trm/
		mkdir AlignmentsTrans.trm.nm/
		mkdir SequencesTrans/
		mkdir SequencesProts/
		mkdir OG_names/
		mkdir RAxMLtrees/
		mkdir SpeciesTree
		mkdir logs

		touch $store/setup.complete
	fi
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
		if test "$(jobs | wc -l)" -ge $threads; then
		        wait -n
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

	trans=$trans
	prots=$prots
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
			if test "$(jobs | wc -l)" -ge $threads; then
			        wait -n
			fi
			{
			base=$(basename ${I%.*})
			echo $base >> dedupe.stats
			bash dedupe.sh -Xmx1g -Xms1g in=$I out=$trans/$base.fna minidentity=99.9 2>&1 | grep 'Input\|Result'  \
				>> dedupe.stats || { echo "dedupe failed for some reason :("; exit 1;}
			cat $trans/$base.fna | grep ">" | sed 's/>//g' > $trans/$base.deduped.names
			filterbyname.sh -Xmx1g -Xms1g --amino include=true in=$prots.preDedup/$base.faa out=$prots/$base.faa names=$trans/$base.deduped.names \
				>> dedupe.stats_long 2>&1
			} &
		done && touch dedupe.complete
		wait
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
			rm $prots.fixed.bk
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
	cd $store || exit
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
	gene_counts="$orthodir/Orthogroups/Orthogroups.GeneCount.tsv"
	OF_alignments=$orthodir/MultipleSequenceAlignments
	all_prots=$wd/all_prots.nm.fa


	alignments_TMP="$orthodir/OG_alignmentsToHMM/alignments/"

	# run this block ass a function to use multiple cores and local variables
	OG_hmm_search () {
		local_wd=$1
		aligned_dir=$2
		all_prots=$3
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
		cat hmmout/${I}.hmmout | grep -ve "^#" |awk '{print $1,$6}' \
			> hmmout/${I}.list
		cat hmmout/${I}.list \
			> hmmout/${I}.list_filter
		# get number of putative paralogs per taxa
		cat hmmout/${I}.list | \
			sed 's/@/\t/g' | awk '{print $1}' | \
			sort | uniq -c | sed 's/^ *//g' | sort -nk1,1 \
			> ${I}.paralogs
		# Init file for HMM-filtering output
		echo -e "OG\thmm_score\tnum_species\tnum_seqs" > hmmout/${I}.filtering_out
		# get number of represented species in hmmout
		local species=$(cat hmmout/${I}.list_filter | awk 'BEGIN { FS="@" } {print $1}' | sort | uniq | wc -l)
		# get total number of seqs in hmm output
		local seqs=$(cat hmmout/${I}.list_filter | wc -l)
		local score=25
		while [ ! $seqs -eq $species ]
		do
	            	cat hmmout/${I}.list_filter | awk -v score="$score" '{if ($2 > score) print $1,$2}' > hmmout/${I}.list_filter_iter
			mv hmmout/${I}.list_filter_iter hmmout/${I}.list_filter
			local species=$(cat hmmout/${I}.list_filter | awk 'BEGIN { FS="@" } {print $1}' | sort | uniq | wc -l)
			local seqs=$(cat hmmout/${I}.list_filter | wc -l)
			echo -e "$I\t$score\t$species\t$seqs" >> hmmout/${I}.filtering_out
			local score=$((score+25))
		done
		echo "${I} ${seqs} ${score}" >> OG_fullset_scores
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
	local_wd=$1
	input_list=$local_wd/$2
	alignments_TMP=$3
	outdir=$4
	all_prots=$5
	threads=$6
        percent=$(( $(cat $input_list | wc -l) / 10))
	echo "Expanding OrthoFinder OGs to full genome set if found in at least $ANI_shortlist_min_OGs taxa" | tee $run_notes
	cd $outdir || exit

	J=0
	K=0
		for I in $(cat $input_list)
		do
			if test "$(jobs | wc -l)" -ge $threads; then
			        wait -n
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
		mkdir OG_alignmentsToHMM
		cd OG_alignmentsToHMM || exit
		mkdir alignments
		touch OG_fullset_scores
		# Pick OGs of interest from Orthofinder genecounts
		#    Writes file OG_SCO_$ANI_shortlist_min_OGs which has one OG ID per line
		echo "#### Filtering OGs for those with seqs from at least $ANI_shortlist_min_OGs taxa with no paralogs ####"
		python $OG_sco_filter $gene_counts $ANI_shortlist_min_OGs
		### Identify orthologs in the full data set from hmms derived
		###    from OG_SCO_$ANI_shortlist_min_OGs
		### making progress indicator
		num_OGs=$(cat OG_SCO_$ANI_shortlist_min_OGs | wc -l)
		percent=$(( num_OGs / 10))
		### Run HMM search with updataed models form full dataset
		###
		cd $orthodir/OG_alignmentsToHMM
		percent=$(( num_OGs / 10))
		J=0
		for I in $(cat OG_SCO_$ANI_shortlist_min_OGs)
		do
			cp $OF_alignments/${I}.fa $alignments_TMP/
			echo -e "Rep_0\t${I}\t"$(cat $alignments_TMP/${I}.fa| grep -c -e ">") \
				>> $orthodir/OG_alignmentsToHMM/hmm_reps.tsv
		done

		for rep in $(seq $reps)
		do
			echo "##############  LIST OF ALIGNMENTS TO START Rep $rep  ############"
			cd $orthodir/OG_alignmentsToHMM
			mkdir hmm_round$rep
			cd hmm_round$rep
			mkdir alignments
               		mkdir hmms
                	mkdir prots
                	mkdir hmmout
                	mkdir SequencesProts
                	mkdir AlignmentsProts
			echo "Running OG expantion rep ${rep}"
			dirtbag_multithreading $orthodir/OG_alignmentsToHMM OG_SCO_$ANI_shortlist_min_OGs $alignments_TMP $orthodir/OG_alignmentsToHMM/hmm_round$rep $all_prots $threads
		done

		#alignments copied to the final module output ($wd/AlignmentsProts/)
		for I in $(cat $orthodir/OG_alignmentsToHMM/OG_SCO_$ANI_shortlist_min_OGs)
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
		# Progress sent to stdout
               	if [ $((J % percent)) -eq 0 ]
                then
			echo $((J/percent*10))" percent of the way through the hmm search"
               	fi
		filter_and_Align_subfunc () {
			local base=$(basename "${i%.*}")
		        cat $i | grep ">" | sed 's/>//g' | sed 's/|.*//g' \
		        > ./OG_names/${base}.names
			# pulling out prot sequences based on names in orthofinder OGs
		        filterbyname.sh -XmX60m -Xms60m include=t \
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
                if [ $(($J % $threads)) -eq 0 ]
                then
                    	wait
                fi
	done
	wait


}


PAL2NAL () {
	echo '
	###################################################
	#####  Use PAL2NAL to produce Nuc alignments ######
	############### from AA alignemtns ################
	###################################################
	#loops of all Orthogroups to#	make file with OG gene names
	#	extract transcripts sequences from a contatentated transcript file
	#	do codon alignments with pal2nal
	'
	date
	func_timing_start
	#aligning transcripts based on orthogroup protien alignments
	cd $wd || exit
	num_OGs=$(ls $wd/AlignmentsProts/OG*.faa | wc -l) # this is 1+ the real num
	percent=$(( num_OGs / 10))
        J=0
	for i in $wd/AlignmentsProts/OG*.faa
	do
		# Progress sent to stdout
                if [ $((J % percent)) -eq 0 ]
                then
			echo $((J/percent*10))" percent of the way through the PAL2NAL"
                fi
		# why the hell did I declare this function with in the loop...probably because of the "local bass=...". this should be passed in with a variable. 
		PAL2NAL_subfunc () {
			#get basename (OG000????)
			local base=$(basename ${i%.*})
			cat $i | grep ">" | sed 's/>//g' | sed 's/|.*//g' \
	        	> ./OG_names/${base}.names
			#pull out trans sequences for each OG
			filterbyname.sh -Xmx60m -Xms60m include=t \
	        	names=./OG_names/${base}.names \
			in=all_trans.nm.fa out=./SequencesTrans/${base}.trans.fa \
			>> $wd/logs/filterbyname.PAL2NAL 2>&1 # changed a "/" to ">" not sure how the code ran before...
			# protein to nuc alignments
			pal2nal.pl $i ./SequencesTrans/${base}.trans.fa -codontable 11 -output fasta \
			> ./AlignmentsTrans/${base}.codon_aln.fa 2> pal2nal.stderr.tmp
		}
		PAL2NAL_subfunc &
                J=$((J+1))
                if [ $(($J % $threads)) -eq 0 ]
                then
                    	wait
		fi
	done
	wait
}


TRIM_TRANS () {
	echo '
	##############################################
	#### trim all Tran alignments with Trimal ####
	####### and remove gene specific names #######
	##############################################
	'
	date
	func_timing_start
	cd $wd/ || exit
	num_OGs=$(ls $wd/AlignmentsTrans/OG*.fa | wc -l) # this is 1+ the real num
        percent=$(( num_OGs / 10))
        J=0
	for I in AlignmentsTrans/OG0*.fa
	do
		# Progress sent to stdout
                if [ $((J % percent)) -eq 0 ]
                then
                    	echo $((J/percent*10))" percent of the way through trimming"
                fi
                TRIM_TRANS_subfunc () {
			base=$(basename ${I%.*.*})
			trimal \
			-in $I -fasta \
			-out AlignmentsTrans.trm/${base}.codon_aln.trm.fa \
			$trimal_parameter > $wd/logs/trimmal_logs 2>&1
			#remove all gene info from header for concat
			#this is really janky: relies on adding the @ during renaming...
			# Removes everything after the @ to leave just the sample name
			cat AlignmentsTrans.trm/${base}.codon_aln.trm.fa \
			| sed 's/@.*$//g' > AlignmentsTrans.trm.nm/${base}.codon_aln.trm.nm.fa
		}
                TRIM_TRANS_subfunc &
               	J=$((J+1))
                if [ $(($J % $threads)) -eq 0 ]
                then
                        wait
                fi
	done
	wait
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
	echo '
	#######################################################
	###### cat alignments 4 genes that are single copy ####
	############# and in at least X samples ###############
	#######################################################
	'
	date
	func_timing_start
	min_num_orthos=$1
	if [ "$ANI" = true ]
	then
		echo "Finding all ANI OGs with $min_num_orthos representetives"
		# make list of SCO $min_num_orthos directly from fasta
		cd $wd || exit
		# Iterate over trimmed alignments and pull out the number of samples per multifasta
		# this has the assumption that there is only one sequence per genome
		# should be enforced by ORTHOFINDER_TO_ALL_SEQS
		# could add an error || exit if it becomes a problem
		for I in AlignmentsTrans.trm.nm/OG0*.codon_aln.trm.nm.fa
		do
			# test transcript alignment for number of homologs (should have no paralogs)
			# if greater than min_num_orthos (samples*min_frac_ortho) 
			# add to OG_SCO_$min_num_orthos list
			final_seqs=$(cat $I | grep -e ">" | wc -l)
			if [ $final_seqs -ge $min_num_orthos ]
			then
				base=$(basename ${I%.*.*.*.*})
				echo $base >> $wd/OG_SCO_$min_num_orthos
			fi
		done
	else
		cd $orthodir/ || exit
		gene_counts="$orthodir/Orthogroups/Orthogroups.GeneCount.tsv"
		echo "Finding all OGs directly from OrthoFinder with $min_num_orthos representetives"
		# finds OGs with at least $min_orthologs SCOs
		# and writes to $orthodir/OG_SCO_$min_num_orthos
		python $OG_sco_filter $gene_counts $min_num_orthos
		mv OG_SCO_$min_num_orthos $wd/OG_SCO_$min_num_orthos
	fi

	cd $wd || exit
	mkdir OG_SCO_$min_num_orthos.align
	echo "Creating directory with alignments of OGs with at least $min_num_orthos orthos"
	for I in $(cat $wd/OG_SCO_$min_num_orthos):
	do
		cp AlignmentsTrans.trm.nm/${I}.codon_aln.trm.nm.fa ./OG_SCO_$min_num_orthos.align/
	done
	echo "Concatenating fasta alignments to phylip format"
	cd $wd/SpeciesTree/ || exit
	perl $catfasta2phyml_cmd -c ../OG_SCO_$min_num_orthos.align/*.fa \
		1> SCO_$min_num_orthos.codon_aln.trm.sco.nm.phy 2> SCO_$min_num_orthos.codon_aln.trm.sco.nm.partitions
	echo "Concatenating fasta alignments to fasta format"
	perl $catfasta2phyml_cmd -f -c ../OG_SCO_$min_num_orthos.align/*.fa \
        1> SCO_$min_num_orthos.codon_aln.trm.sco.nm.fa 2>> $store/verbose_log.txt
}

#####################################
#####################################

SCO_strict () {
	echo '
        ####################################################
        ############ Identify Single Copy Orthos ###########
        ########### And Place in SCO_strict.align ##########
        #### and make concattenated alignment for RAxML ####
        ####################################################
        '
	date
	func_timing_start
	OG_SCO_strict=$wd/OG_SCO_strict
	if [ -f $OG_SCO_strict ]
	then
		rm $OG_SCO_strict
	fi
	cd $wd || exit
	mkdir OG_SCO_strict.align
	if [ $ANI == "true" ]
	then
		echo "Finding all proper SCOs from ANI Orthogroups"
		# make list of all proper SCOs directly from fasta
		cd $wd || exit
		# Iterate over trimmed alignments and pull out the number of samples per multifasta
		# this has the assumption that there is only one sequence per genome
		# should be enforced by ORTHOFINDER_TO_ALL_SEQS
		# could add an error || exit if it becomes a problem
		for I in AlignmentsTrans.trm.nm/OG0*.codon_aln.trm.nm.fa
		do
			final_seqs=$(cat $I | grep -e ">" | wc -l)
			if [ $final_seqs -ge $(cat $store/all_input_list | wc -l) ]
			then
				# add SCO_strict alignemtns to a specific DIR
				cp $I ./OG_SCO_strict.align/
				base=$(basename ${I%.*.*.*.*})
				echo $base >> $OG_SCO_strict
			fi
		done
	else
		for I in $orthodir/Single_Copy_Orthologue_Sequences/OG0*.fa
		do
		  	base=$(basename ${I%.*})
		  	cp AlignmentsTrans.trm.nm/${base}.codon_aln.trm.nm.fa ./OG_SCO_strict.align/
			#used later for ASTRAL gene tree selection
			echo $base >> $OG_SCO_strict
		done
	fi
	#cat SCO_strict nuc alignments
	cd $wd/SpeciesTree/ || exit
	perl $catfasta2phyml_cmd -c $wd/OG_SCO_strict.align/*.fa \
		1> SCO_strict.codon_aln.trm.sco.nm.phy 2> SCO_strict.codon_aln.trm.sco.nm.partitions.tmp
	perl $catfasta2phyml_cmd -f -c $wd/OG_SCO_strict.align/*.fa \
        1> SCO_strict.codon_aln.trm.sco.nm.fa 2>> $store/verbose_log.txt
	cat SCO_strict.codon_aln.trm.sco.nm.partitions.tmp | awk '{print "DNA,\t"$0}' \
		> SCO_strict.codon_aln.trm.sco.nm.partitions && rm SCO_strict.codon_aln.trm.sco.nm.partitions.tmp
}


TREE_BUILD () {
	echo '
	##########################################
	###### Build trees and bootstrap #########
	####### from SCO_$min_num_orthos #########
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
	output_name=$(basename ${input_alignment%.*.*.*.*})
	partition_file=${input_alignment%.*}.partitions
	cd $output_dir || exit
	# subfunction to run RAxml
	RAxML_run () {
		mkdir raxml
		cd raxml
		raxmlHPC-PTHREADS $RAxML_speciestree_options \
		-s $input_alignment \
		-q $partition_file \
		-n ${output_name} \
		> ./RAxML_output.1 2>&1
		treefile=$(ls RAxML_bipartitionsBranchLabels*)
		mv $treefile ../${treefile}.tree
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
		mv $treefile ./${treefile%.*}.tree
		cd $output_dir || exit
	}
	IQTREE_run () {
		mkdir iqtree
		cd iqtree
		iqtree2 -s $input_alignment \
		--prefix iqtree.${output_name} \
		-Q $partition_file -m MFP+MERGE \
		-T $threads --seed 1234 > iqtree.long_log
		treefile=$(ls *.treefile)
		mv $treefile ../${treefile%.*}.tree
		cd $output_dir || exit

	}

	# decide which tree method to use for the cancatenate gene nuc matrix
	

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
		FASTTREE_run
	fi
	if [[ " ${tree_method[*]} " =~ " iqtree " ]]
	then
		echo "#################################"
		echo "Building species tree with IQtree"
		echo "#################################"
		IQTREE_run
	fi
	
}

orthofinderGENE2SPECIES_TREE () {
	echo '
	#############################################################
	Species tree reconstruction from all Orthogroup protien trees
	#############################################################
	Running Astral on all ortholog from orthofinder
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

allTransGENE_TREEs () {
	echo '
	##################################################
	###### Build Gene trees for all Transcript #######
	########## Alignments to use for Astral ##########
	############# Species Tree Estemation ############
	##################################################
	'
	date
	func_timing_start
	if [ -d $wd/trans_gene_trees/ ]
	then
	   rm -r ${wd}/trans_gene_trees.bk
	   mv $wd/trans_gene_trees/ ${wd}/trans_gene_trees.bk
	fi
	mkdir $wd/trans_gene_trees/
	mkdir $wd/trans_gene_trees.nm/
	cd $wd/trans_gene_trees || exit
	num_OGs=$(ls $wd/AlignmentsTrans.trm/OG00*.codon_aln.trm.fa | wc -l) # this is 1+ the real num
        percent=$(( num_OGs / 10))
        J=0
	for i in $wd/AlignmentsTrans.trm/OG00*.codon_aln.trm.fa
	do
		# Progress sent to stdout
                if [ $((J % percent)) -eq 0 ]
                then
                    	echo $((J/percent*10))" percent of the way through building trans trees"
                fi
		TRANS_TREE_subfunc () {
		        base=$(basename ${i%.*.*.*})
		        #add partition finder?
		        #raxmlHPC-PTHREADS -T $threads -s $wd/AlignmentsTrans.trm/${base}.codon_aln.trm.fa \
			#-n ${base}.tree -m GTRGAMMA
		        #>> ./RAxML_output.tmp
        		if [[ " ${gene_tree_methods[*]} " =~ " fasttree " ]]
        		then
				fasttree -gtr -quiet -gamma \
				-nt $wd/AlignmentsTrans.trm/${base}.codon_aln.trm.fa \
				> ./${base}.fasttree.tree
				cat ./${base}.fasttree.tree | sed 's/@[^:]*:/:/g' \
				> $wd/trans_gene_trees.nm/${base}.fasttree.tree
			elif [[ " ${gene_tree_methods[*]} " =~ " iqtree " ]]
			then
				iqtree -s $wd/AlignmentsTrans.trm/${base}.codon_aln.trm.fa \
				--prefix ${base} -T 1 > ${base}.iqtree.log
				mv ./${base}.treefile ./${base}.iqtree.tree
				cat ./${base}.iqtree.tree | sed 's/@[^:]*:/:/g' \
				> $wd/trans_gene_trees.nm/${base}.iqtree.tree
			else
			echo "Gene tree estemation from gene alignemtns not done; gene_tree_method not set to either fasttree, raxml or iqtree" 
			exit 1
        		fi
		}

                TRANS_TREE_subfunc &
                J=$((J+1))
                if [ $(($J % $threads)) -eq 0 ]
                then
                        wait
                fi

	done
}



astral_allTransGENE2SPECIES_TREE () {
	echo '
	##########################################################
	Species tree reconstruction from all trans alignments
	##########################################################
	'

	#astral all genes
	#TODO
	echo "Does nothing!
	TODO
	Need to set it up to run astral_p to handle paralogs
	I think this will only be useful in limited cases:
	back burner
	ALSO, functionality could be put into astral_TransGENE2SPECIES_TREE
	just need to change astral_CMD
	"
}


astral_TransGENE2SPECIES_TREE () {
	echo '
	###############################################################
	Species tree reconstruction from SCO_min_align trans alignments
	###############################################################
	'
	date
	func_timing_start
	local genelist=$1
	echo "
	Running Astral on ${genelist} genelist
	"
	#inputs
	# in_dir=$2
	#genelist=$orthodir/OG_SCO_$min_num_orthos
	local in_dir=$wd/trans_gene_trees.nm
	local base_list=$(basename ${genelist})
	local ASTRAL_cmd=$ASTRAL_cmd
	#outputs
	local out_dir=astral_trees
	local excluded=$wd/$out_dir/trimmed_from_${base_list}_gene_list
	local concat_trees=$wd/$out_dir/${base_list}

	cd $wd || exit
	mkdir $out_dir
	# use iqtree or fasttree (or both)
	for method in $gene_tree_methods
	do
		concat_tree_method=$concat_trees.${method}.trees
		astral_tree=astral.${base_list}.${method}.tree
		#remove concatenated tree file if it exists.
		if [ -f $concat_tree_method ]
		then
		   rm -r $concat_tree_method.bk
		   mv $concat_tree_method concat_tree_method.bk
		fi

		# remove list of genes ommited because trimming reduced the number of leaves to > 4
		rm $excluded

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
			mv $store/$I $store/tmp_trash
		done
		rm -r $store/tmp_trash
	fi
	# but keep AlignmentTrans.trm.nm ans AlignemntProt
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
