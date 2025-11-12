
SET_UP_ADDASM_DIRS () {
    cd $store || exit
    mkdir $addasm_dir   
    cd $addasm_dir || exit
    mkdir $genome_dir
    mkdir $prots
    mkdir $trans
    mkdir $annots
    mkdir $new_hmm_output
    mkdir $new_OG_prots
    mkdir $new_OG_CDS
    mkdir $new_prot_alignments
    mkdir $new_prot_alignments.trm
    mkdir $new_prot_alignments.trm.nm
    mkdir $new_CDS_alignments
    mkdir $new_CDS_alignments.trm
    mkdir $new_CDS_alignments.trm.nm
    mkdir $OG_names
    mkdir $OLD_TREE_INFO
    mkdir $SpeciesTree
    mkdir $new_trees
    mkdir $logs
    mkdir trimmed_columns # need to make a variable somewhere...
}

GET_OLD_TREE_INFO () {
    echo "
    ##################################################
    ######## Setting variables based on old ##########
    ###### files if not declared on CMD line #########
    ##################################################"
    ################################
    #### set phylogenetic tool #####
    ###### For tree building #######
    ################################
    OLD_tree_method=( $(ls $store/FINAL_SPECIES_TREES/ | \
			grep -v compare | \
			awk -F"." '{print $1}' | \
			sort | uniq) )
    if [[ -z ${tree_method+x} ]]
    then
        echo "Setting tree methods based on old trees"
        tree_method=(${OLD_tree_method[@]})
    else
        diff=$(echo ${tree_method[@]} ${OLD_tree_method[@]} | tr ' ' '\n' | sort | uniq -u )
        if [ $(echo ${tree_method[@]} ${diff[@]} | tr ' ' '\n' | sort | uniq -D | uniq | wc -l ) -gt 0 ]
        then
            echo "WARNING: Tree methods declared on the command line with \"-p " ${tree_method[@]}"\""
            echo "Does not match available trees from previous OrthoPhyl Run!!!"
            exit 1  
        fi
    fi
    
    # test is FastTree was used.
    if [[ " ${tree_method[*]} " =~ " fastTree " ]] || [[ " ${tree_method[*]} " =~ " astral " ]]
    then
        echo "#######################################################################################"
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "STILL IN DEVELOPMENT"
        echo "Adding sequences to FastTree trees from original Orthophyl Run is still in development."
        echo "   This should be functional now, but some more testing is needed (big trees, low support and such)"
        echo "Tree(s) built by iqtree will be added to."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "#######################################################################################"
        sleep 10
    fi

    ################################
    #### set the SCO set to use ####
    ###### For tree building #######
    ################################
    if [[ -z ${SCO_sets+x} ]]
    then
        echo "Setting SCO sets to use based on old trees"
		SCO_sets=( $(ls $store/FINAL_SPECIES_TREES/ | \
			grep -v compare | \
			awk -F"." '{print $2}' | \
			sort | uniq) )
	fi	
	
    ################################
    #### set which data to use #####
    ###### For tree building #######
    ################################
    OLD_TREE_DATA_list=( $(ls $store/FINAL_SPECIES_TREES/ | \
			grep -v compare | \
			awk -F"." '{print $3}' | \
			sort | uniq) )
    if [[ -z ${TREE_DATA+x} ]]
    then
        echo "Setting to use CDS, PROTs, or BOTH based on old trees"
		TREE_DATA_list=(${OLD_TREE_DATA_list[@]})
    else 
        if [ "$TREE_DATA" = "BOTH" ]
        then
            TREE_DATA_list=("CDS" "PROT")
        else
            TREE_DATA_list=$TREE_DATA
        fi
        diff=$(echo ${TREE_DATA_list[@]} ${OLD_TREE_DATA_list[@]} | tr ' ' '\n' | sort | uniq -u )
        if [ $(echo ${TREE_DATA_list[@]} ${diff[@]} | tr ' ' '\n' | sort | uniq -D | uniq | wc -l ) -gt 0 ]
        then
            echo "WARNING: Tree methods declared on the command line with \"-o " ${TREE_DATA}"\""
            echo "Does not match available trees from previous OrthoPhyl Run!!!"
            exit 1  
        fi
    fi	

    

    # print v ariable telling addasm which trees to add to
    echo "tree_method = ${tree_method[@]}"
    echo "SCO_sets = ${SCO_sets[@]}"
    echo "TREE_DATA_list = ${TREE_DATA_list[@]}"

    echo "
    #########################################
    # check that expected trees are present #
    #########################################
    "
    for I in ${tree_method[@]}
    do
        if [ $I == "astral" ] # fix the ASTRAL tree test...not important until we add seqs to ASTAL trees
        then   
            echo "Skipping ASTRAL tree test for now"
        else
            for J in ${SCO_sets[@]}
            do
                for K in ${TREE_DATA_list[@]}
                do
                    if [ ! -f $store/FINAL_SPECIES_TREES/$I.$J.$K*.tree ] # * is to deal with gene tree methods used with astral 
                    then    
                        echo -e "Tree expected at:\n$store/FINAL_SPECIES_TREES/$I.$J.$K.*.tree\n was not found!!!\nIt's hard to add sequences to a tree that doesn't exist"
                        exit 1
                    fi
                done
            done
        fi
    done
}

ANNOTATIONS (){
    # runs PRODIGAL_PREDICT from main functions.sh 
    if [[ ${genomes_provided+x} ]]
    then
        # defined in functions.sh
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
    # defined in functions.sh
    DEDUP_annot_trans $addasm_dir $trans $prots 
    FIX_TRANS_NAMES $trans
    FIX_PROTS_NAMES $prots
}

HMM_search () {
    echo "
    ##########################
    ### Running HMM_search ###
    ##########################
    "
    echo -e "\nRun hmmsearch for each precomputed OG in:" 
    echo -e "$hmm_dir\n"
    new_hmm_output=$1
    hmm_dir=$2
    all_prots=$3
    new_OG_prots=$4
    cd ${new_hmm_output} || exit 
    num_OGs=$(ls $hmm_dir/*.hmm | wc -l) # this is 1+ the real num
    percent=$(( num_OGs / 10))
    J=0
    for I in $hmm_dir/*.hmm
    do
        # If number of multithreaded jobs is greater than $threads, wait.
        #  Test again if a job finished
    	if test "$(jobs | wc -l)" -ge $threads; then
			wait -n
            trap control_c INT
		fi
        # Progress sent to stdout
        if [ $((J % percent )) -eq 0 ]
        then
			echo $((J/percent*10))" percent of the way through the $FUNCNAME"
        fi
        hmm_search_subfunc () {
            base=$(basename ${I%.*})
            hmmsearch -T 25 \
            -o ../ORTHOFINDER_TO_ALL_SEQS.out \
            --tblout ${base}.hmmout \
            $I $all_prots
            ##################################################
            # This block maximizes hmm score filter, to eliminate paraogs
            # while keeping all the orthos possible
            local no_para_score=$(cat ${base}.hmmout | grep -ve "^#" | sed 's/@/ /g' | awk '{print $1,$7}' |\
                sort -k1,1 -k2,2rn |\
                awk '{if ($1==prev) print $1,$2} {prev=$1}' |\
                sort -rnk2 |\
                head -n 1 | \
                awk '{print $2}')
            cat ${base}.hmmout |\
                grep -ve "^#" | \
                awk -v no_para_score=$no_para_score '{if ($6>no_para_score) print $1,$6}' |\
                sort -k1 \
                > ${base}.list_filter
            cat ${base}.list_filter | awk '{print $1}' > ${base}.names
            filterbyname.sh -Xmx60m -Xms60m overwrite=True include=True ignorejunk=True \
                names=${base}.names \
                in=$all_prots \
                out=${new_OG_prots}/${base}.faa >> filterbyname.so_verbose.out 2>&1
        }
        hmm_search_subfunc &
        # Add one to J for tracking multithreading
        J=$((J+1))
    done
    wait
}

ADD_2_ALIGNMENTS () {
    ################################################
    # take old alignment and add new seqs
    ##################################################
    echo "
        ###############################
        ### Running ADD_2_alignment ###
        ###############################
        "
    local OLD_alignments=$1
    local new_OG_seqs=$2
    local new_alignments=$3
    local threads=$4
    
    num_alignments=$(ls $OLD_alignments/OG* | wc -l) # this is 1+ the real num
    percent=$(( num_alignments / 10))

    J=0
    for OG_alignment in $OLD_alignments/OG*
    do
        # If number of multithreaded jobs is greater than $threads, wait.
        #  Test again if a job finishes
    	if test "$(jobs | wc -l)" -ge $threads; then
			wait -n
            trap control_c INT
		fi
        # Progress sent to stdout
        if [ $((J % percent)) -eq 0 ]
        then
			echo $((J/percent*10))" percent of the way through the $FUNCNAME"
        fi
        base=$(basename $OG_alignment)
        OG=${base%%.*}
        if [ -s $new_OG_seqs/${OG}.fa* ]
        then
            #echo "adding new sequences to $OG"
            sed 's/!/N/g' -i $OG_alignment
            mafft --add $new_OG_seqs/${OG}.fa* --keeplength $OG_alignment \
                > ${new_alignments}/${OG}.new_aligned.fa \
                2>> $addasm_dir/logs/${OG}.mafft_out &
        else
            #echo "NOT adding to $OG"
            ln -sr $OG_alignment ${new_alignments}/${OG}.new_aligned.fa
        fi
    J=$((J+1))
    done
    wait
}


ADD_2_IQTREE () {

    echo "
    ###############################
    ### RUNNING ADD_2_IQTREE!!! ###
    ###############################
    "
    new_alignment=$1 # new concatenated alignment with same length as one that generated $old_tree
    old_tree=$2 # old tree to contrain the new one
    old_scheme=$3 # tree_file.best_scheme from old iqtree run
    new_prefix=$4 # new tree location/name
    threads=$5

    if [ -f $old_scheme ]
    then
        partition_scheme="-p $old_scheme"
    fi
    
    iqtree -s $new_alignment \
        -g $old_tree \
        $partition_scheme \
        --prefix $new_prefix \
        -T AUTO --threads-max $threads
}
ADD_2_fasttree () {
    echo "#################################"
    echo "### RUNNING ADD_2_FASTTREE!!! ###"
    echo "#################################"
    local new_alignment=$1 # new concatenated alignment with same length as one that generated $old_tree
    local old_tree=$2 # old tree to contrain the new one
    local new_prefix=$3 # new tree location/name
    local tree_data=$4
    local threads=$5

    # lots to do before this works
    #need to get bioperl to work in my mamba environment
     #works in a clean env...
    #need to collaps bad splits...dont think I can do it during constraint conversion
    
    if [[ $tree_data == "CDS" ]]
    then
        fasttree_speciestree_options=$fasttree_CDS_speciestree_options
        alignemnt="-nt "$new_alignment
    elif [[ $J == "PROT" ]]
    then
        fasttree_speciestree_options=$fasttree_PROT_speciestree_options
        alignment=$new_alignment
    fi
    
    #need to get bioperl to work in my mamba environment
     #works in a clean env...
    #need to collaps bad splits...dont think I can do it during constraint conversion
    echo " Creating contraint tree from Earlier OP run"
    constraint_tree=$addasm_dir/$(basename $old_tree).contraints
    python $script_home/python_scripts/Newick2FastTreeConstraints.py $old_tree > $constraint_tree.fa
    cat $I | awk -vRS=">" -vFS="\n" -vOFS="" \
        '$0!=""{$1=substr($1,1,15);$1=sprintf ("%-17s",$1)}$0!=""' \
        > TMP.phy1
    local num=$(cat TMP.phy1 | wc -l)
    local len=$(cat TMP.phy1 | head -n 1 | awk '{print length($2)}')
    echo -e "\t"$num"\t"$len > $constraint_tree
    cat TMP.phy1 >> $constraint_tree
    rm TMP.phy1
    echo "Running FastTree..."
    FastTreeMP $fasttree_speciestree_options \
		-out ./fastTree.${output_name}.tree_working \
         -constraints $constraint_tree \
        $alignemnt
}


RUN_ML_SUPERMATRIX_WORKFLOW () {
    ###################################################
    # Run ML 
    ####################################################
    # Run Prot workflow
    
    for SCO_set in ${SCO_sets[@]}
    do
        ln -s $store/phylo_current/$SCO_set $addasm_dir/$SCO_set
        if [ ! -f $store/phylo_current/$SCO_set ]
        then
            echo "WARNING: SCO_sets indicate a file at $store/phylo_current/$SCO_set that doesnt exist"
            exit 1
        fi
    done

    for data in ${TREE_DATA_list[@]}
    do
        
    
        if [ $data == "PROT" ]
        then
            new_alignments=$new_prot_alignments.trm.nm
        elif [ $data == "CDS" ]
        then   
            new_alignments=$new_CDS_alignments.trm.nm
        fi


        for SCO_set in ${SCO_sets[@]}
        do
        cat_alignments \
            $addasm_dir/$SCO_set \
            $new_alignments \
            $AlignmentsConcatenated \
            $data
        #ALIGNMENT_STATS $wd/SCO_${min_frac_orthos}.align
        done
    done

    if [ $total_tree_rerun == "True" ]
    then
        for I in $AlignmentsConcatenated/*.$data.*.phy
        do
            TREE_BUILD $SpeciesTree $I $threads $data
        done
    else
        if [[ " ${tree_method[*]} " =~ " iqtree " ]]
        then
            for I in ${SCO_sets[@]}
            do
                for J in ${TREE_DATA_list[@]}
                do
                ADD_2_IQTREE \
                    $AlignmentsConcatenated/$I.$J.trm.sco.nm.fa \
                    $OLD_TREES/iqtree.$I.$J.tree \
                    $OLD_iqtree_output/iqtree.$I.$J.best_scheme \
                    $new_trees/iqtree.$I.$J.addasm \
                    $threads
                done
            done
        fi
        if [[ " ${tree_method[*]} " =~ " fastTree " ]]
        then
            for I in ${SCO_sets[@]}
            do
                for J in ${TREE_DATA_list[@]}
                do
                
                ADD_2_fasttree \
                    $AlignmentsConcatenated/$I.$J.trm.sco.nm.fa \
                    $OLD_TREES/fastTree.$I.$J.tree \
                    $new_trees/fastTree.$I.$J.addasm \
                    $J \
                    $threads
                done
            done
        fi 
    fi
 
}

ASTRAL_WORKFLOW () {
    
    if [[ "$TREE_DATA" = "CDS" || "$TREE_DATA" = "BOTH" ]]
    then        
        tree_data=CDS
        # build gene trees for all CDS alignments
        allGENE_TREEs $tree_data $wd/AlignmentsCDS.trm.nm $wd/${tree_data}_gene_trees

        # astral_allCDSGENE2SPECIES_TREE #Still not written (needs a different ASTRAL )
        
        # build ASTRAL tree from SCO_relaxed gene trees
        if [[ "$relaxed" != false ]]
            then
            astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_$min_frac_orthos
        fi

        # Build ASTRAL tree from SCO_strict
        astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_strict
    fi
    if [[ "$TREE_DATA" = "PROT" || "$TREE_DATA" = "BOTH" ]]
    then
        tree_data=PROT
        # build gene trees for all CDS alignments
        allGENE_TREEs $tree_data $wd/AlignmentsProts.trm.nm $wd/${tree_data}_gene_trees 

        # astral_allCDSGENE2SPECIES_TREE #Still not written (needs a different ASTRAL )
        
        # build ASTRAL tree from SCO_relaxed gene trees
        if [[ "$relaxed" != false ]]
        then
            astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_$min_frac_orthos
        fi

        # Build ASTRAL tree from SCO_strict
        astral_GENE2SPECIES_TREE $tree_data $wd/${tree_data}_gene_trees $wd/SCO_strict
    fi
}