#!/bin/bash

control_c () {
    echo "Keybourd interupt!!!"
    #kill 
    exit 1
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
        echo "Args comming in through file set with -c, other args seem good"
    fi

    

}

SET_RIGOR () {
    echo "
    #############################################
    ##### Setting methods for tree building #####
    #############################################
	"

    echo "Setting trees to build based of provided \"-R,--Rigor\" declared on the command line"
    echo "rigor = ${rigor}"
    echo "*********************************************"
    if [[ $rigor = "full" ]]
    then    
        echo "Building the biggest/baddest trees with IQTree and ASTRAL*"
        echo "  Good for: Supporting taxonomy changes, arguing with reviewer 2, and confusing yourself"
        echo "  Use both CDS and PROT alignments* (*helps* with saturation issues)"
        echo "  Strict and relaxed single copy orthologs*"
        echo "  Full partition model merging for IQTree to reducing overfitting (takes for-ever)"
        echo "  *(unless overridden on command line)"
        export tree_method=("iqtree" "astral")
        export IQtree_CDS_partition_options="-m MFP+MERGE --rclusterf 10" # cluster partitions with FreeRate models using a fast$
        export IQtree_PROT_partition_options="-m MFP+MERGE --rclusterf 10" # cluster partitions with FreeRate models using a fas$
        export gene_tree_methods=("iqtree")
        export iqtree_CDS_genetree_options="-m MFP --lmap 2000 --symtest -B 1000 -t PARS --ninit 2 "
        export iqtree_PROT_genetree_options="-m MFP --lmap 2000 --symtest -B 1000 -t PARS --ninit 2 "

    elif [[ $rigor = "fast" ]]
    then    
        echo "Building trees with fast methods (fasttree and ASTRAL)"
        echo "  Good for: exploritory work, quick turnarounds, arguing with committee members" 
        echo "  Will use only PROT or CDS alignments (default PROT)*"
        echo "  Only \"strict\" single copy orthologs will be used*"
        echo "  *(unless overridden on command line)"
        export tree_method=("fasttree" "astral")
        export fasttree_PROT_speciestree_options=""
        export gene_tree_methods=("fasttree")
        if [[ "$ARGS_SET" == *o* ]]
        then 
            echo "  Using "$TREE_DATA" alignments for tree generation!" 
        else
            TREE_DATA=PROT
            echo "  Using "$TREE_DATA" alignments for tree generation!" 
        fi
        if [[ "$ARGS_SET" == *m* ]]
        then 
            echo "  Using OGs found in at least"$min_num_orthos" fracton on samples" 
        else
            min_num_orthos=1
            echo "  Only using OGs found in "$min_num_orthos" all samples" 
        fi

    elif [[ $rigor = "medium" ]]
    then
        echo "Building trees with iqtree and ASTRAL"
        echo "  Good for: "
        echo "  iqtree will test GTR and FreeRate models for CDS and standard models for Prots"
        echo "     but will not merge patitions (genes)"
        export tree_method=("iqtree" "astral")
        export IQtree_CDS_partition_options="-m MFP" # test GTR and FreeRate models
        export IQtree_PROT_partition_options="" # tests standard PROT models
        export gene_tree_methods=("iqtree")
        export iqtree_CDS_genetree_options="--lmap 2000 --symtest -B 1000 -t PARS --ninit 2 "
        export iqtree_PROT_genetree_options="--lmap 2000 --symtest -B 1000 -t PARS --ninit 2 "
    else
        echo "You used a -R/--rigor option that is not recognized! \n options are full, medium, or fast"
        exit 1
    fi
    
}