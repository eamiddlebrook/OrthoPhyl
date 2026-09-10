#!/bin/bash
# USAGE:
# utils/gather_filter_local_asms.sh taxon local_db_dir output_dir threads
#
# This script searches a local genome database instead of downloading from NCBI.
# It mirrors the workflow of gather_filter_asms.sh but uses pre-downloaded genomes.
#
# Local database structure expected:
# local_db_dir/
#   ├── genomes/              # All genome FASTA files
#   │   ├── GCF_*.fna.gz
#   │   └── GCA_*.fna.gz
#   ├── assembly_summary.txt  # NCBI assembly summary (optional)
#   └── metadata/             # Assembly metadata (optional)
#       └── All-*-info.assembly.xml

# Inside a Singularity/Docker container the gather_genomes env is already on
#   PATH (see the container %environment section), and $HOME maps to the host
#   home -- sourcing the host .bash_profile and running `conda activate` there
#   would (wrongly) put the *host* checkm2 on PATH, which cannot execute inside
#   the container ("required file not found"). So only activate the env when
#   NOT in a container, matching the guard in OrthoPhyl.sh.
if [[ -z ${SINGULARITY_CONTAINER+x} ]] && [[ -z ${DOCKER+x} ]]
then
	source ~/.bash_profile
	echo "To create conda env use
conda create -n gather_genomes \\
-c bioconda -c conda-forge \\
checkm2 bbmap entrez-direct ncbi-datasets-cli
then download the CheckM2 DIAMOND database once with:
checkm2 database --download
"
	conda activate gather_genomes || exit
fi

export taxon="$1"
export local_db_dir="$2"
export wd="$3"
threads="$4"

if [ -z "$taxon" ] || [ -z "$local_db_dir" ] || [ -z "$wd" ] || [ -z "$threads" ]; then
    echo "Usage: $0 taxon local_db_dir output_dir threads"
    echo "Example: $0 Rhizobium /data/local_genomes/ my_output/ 8"
    exit 1
fi

# Convert paths to absolute
if [[ "$local_db_dir" = ~* ]]; then
    local_db_dir="${local_db_dir/#\~/$HOME}"
fi
if [[ "$local_db_dir" != /* ]]; then
    local_db_dir="$(cd "$local_db_dir" 2>/dev/null && pwd)"
fi

if [[ "$wd" = ~* ]]; then
    wd="${wd/#\~/$HOME}"
fi
if [[ "$wd" != /* ]]; then
    wd="$(cd "$(dirname "$wd")" 2>/dev/null && pwd)/$(basename "$wd")"
fi

echo "======================================"
echo "Local Genome Database Search"
echo "======================================"
echo "Taxon: $taxon"
echo "Local DB: $local_db_dir"
echo "Output: $wd"
echo "Threads: $threads"
echo ""

# Validate local database
if [ ! -d "$local_db_dir" ]; then
    echo "ERROR: Local database directory not found: $local_db_dir"
    exit 1
fi

genomes_dir="$local_db_dir/genomes"
if [ ! -d "$genomes_dir" ]; then
    echo "ERROR: genomes/ subdirectory not found in $local_db_dir"
    echo "Expected structure:"
    echo "  $local_db_dir/genomes/GCF_*.fna.gz"
    exit 1
fi

# Create output directory
mkdir -p "$wd"
cd "$wd" || exit

# Pick a temp dir that is BOTH writable AND short, then point every temp-dir
#   var at it. Two container constraints collide here:
#     * Singularity's image /tmp is read-only, which breaks edirect's
#       nquire/mktemp ("mktemp: ... Read-only file system") -> needs writable.
#     * CheckM's multiprocessing.Manager binds an AF_UNIX socket under $TMPDIR,
#       and the socket path has a hard 108-char kernel limit. A temp dir deep
#       under $wd overflows it ("OSError: AF_UNIX path too long") -> needs short.
#   /dev/shm is writable (tmpfs) and short in a container; fall back through
#   /tmp, $HOME, and finally $wd for non-container / unusual setups.
op_tmp=""
for base in /dev/shm /tmp "$HOME" "$wd"; do
	if [ -d "$base" ] && [ -w "$base" ]; then
		op_tmp=$(mktemp -d "$base/op_tmp.XXXXXX" 2>/dev/null) && break
	fi
done
if [ -z "$op_tmp" ]; then
	echo "ERROR: could not create a writable temp dir (tried /dev/shm /tmp \$HOME \$wd)" >&2
	exit 1
fi
export TMPDIR="$op_tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
# clean up the scratch dir when the script exits (any reason)
trap 'rm -rf "$op_tmp"' EXIT
echo "Using temp dir (writable + short for CheckM AF_UNIX sockets): $TMPDIR"

echo "Searching for $taxon genomes in local database..."

#######################################################
#### Main function for local genome gathering ####
#######################################################

main () {
    get_local_metadata
    filter_by_taxon
    filter_NCBI_genomes
    link_local_assemblies
    get_all_asm_list
    get_biosample_GEOdata
    merge_metadata_geoloc
    all_sample_metadata
    aggregate_assemblies "$wd"/assemblies_all.TMP
    get_stats_with_checkM "$wd"/assemblies_all.TMP/ $wd/checkM_out genome fna
    filter_asm_by_stats default default default default default default 95 1.0
    get_all_asms_to_remove
    filter_for_redundancy
}

get_local_metadata () {
    echo ""
    echo "#########################################"
    echo "##### Get metadata for local genomes ####"
    echo "#########################################"
    echo ""
    
    # Check if we have pre-downloaded metadata
    metadata_file="$local_db_dir/metadata/All-$taxon-info.assembly.xml"
    
    if [ -f "$metadata_file" ]; then
        echo "Using pre-downloaded metadata: $metadata_file"
        cp "$metadata_file" All-$taxon-info.assembly.xml
    else
        echo "Downloading metadata from NCBI (genomes will be used from local DB)..."
        esearch -db assembly -query "txid${taxon}[organism]" | esummary \
            > All-$taxon-info.assembly.xml
    fi
    
    # Parse metadata
    cat All-$taxon-info.assembly.xml \
        | xtract -pattern DocumentSummary \
            -def "NA" -element BioSampleAccn RefSeq Genbank SpeciesName Sub_value FtpPath_GenBank FtpPath_RefSeq Taxid taxonomy-check-status ExclFromRefSeq | \
            sed 's/ /_/g' \
            > All-$taxon.assembly.BS_to_meta
    
    n_assemblies=$(cat All-$taxon.assembly.BS_to_meta | wc -l)
    echo "Found $n_assemblies assemblies in NCBI database for $taxon"
}

filter_by_taxon () {
    echo ""
    echo "#########################################"
    echo "#### Find matching genomes in local DB ##"
    echo "#########################################"
    echo ""
    
    mkdir -p assemblies_datasets_uniq
    
    # Get list of accessions from metadata
    cat All-$taxon.assembly.BS_to_meta | awk '{print $2}' | grep GCF > assemblies_expected_refseq.names
    cat All-$taxon.assembly.BS_to_meta | awk '{print $3}' | grep GCA > assemblies_expected_genbank.names
    cat assemblies_expected_refseq.names assemblies_expected_genbank.names > assemblies_expected.names
    
    n_expected=$(cat assemblies_expected.names | wc -l)
    echo "Expected $n_expected assemblies based on metadata"
    
    # Find which ones exist in local database
    :> assemblies_found.names
    
    while read acc; do
        # Try to find this accession in local database
        # Look for exact match or version variants (GCF_000123.1, GCF_000123.2, etc.)
        base_acc=$(echo $acc | sed 's/\.[0-9]*$//')
        
        found=false
        for genome_file in "$genomes_dir"/${base_acc}*.fna.gz "$genomes_dir"/${base_acc}*.fna; do
            if [ -f "$genome_file" ]; then
                # Extract accession from filename
                found_acc=$(basename "$genome_file" | sed 's/\.[^.]*$//' | sed 's/_genomic$//')
                echo "$found_acc" >> assemblies_found.names
                
                # Create symlink in working directory
                ln -sf "$genome_file" assemblies_datasets_uniq/${found_acc}.fna.gz
                found=true
                break
            fi
        done
        
        if [ "$found" = "false" ]; then
            echo "  Missing: $acc" >> assemblies_missing.names
        fi
    done < assemblies_expected.names
    
    n_found=$(cat assemblies_found.names | wc -l)
    n_missing=$(cat assemblies_missing.names 2>/dev/null | wc -l || echo 0)
    
    echo "✓ Found $n_found genomes in local database"
    
    if [ $n_missing -gt 0 ]; then
        echo "⚠ Missing $n_missing genomes (listed in assemblies_missing.names)"
        echo "  You may want to download these separately"
    fi
    
    if [ $n_found -eq 0 ]; then
        echo "ERROR: No matching genomes found in local database"
        echo "Check that:"
        echo "  1. Genomes are in: $genomes_dir"
        echo "  2. Filenames match NCBI accessions (GCF_*.fna.gz or GCA_*.fna.gz)"
        exit 1
    fi
    
    # Use found assemblies as our dataset
    cp assemblies_found.names assemblies_datasets.names
}

filter_NCBI_genomes () {
    echo ""
    echo "#################################################"
    echo "####### Remove RefSeq/Genbank redundancy ########"
    echo "####### Prefering to keep RefSeq entries ########"
    echo "#################################################"
    echo ""
    
    # Same logic as original script
    cat assemblies_datasets.names | grep GCF > assemblies_datasets_refseq.names
    cat assemblies_datasets.names | grep GCA > assemblies_datasets_genbank.names
    cat assemblies_datasets_refseq.names | sed 's/GCF_//g' > assemblies_datasets_refseq.edit.names
    cat assemblies_datasets_genbank.names | sed 's/GCA_//g' > assemblies_datasets_genbank.edit.names
    comm -13 assemblies_datasets_refseq.edit.names assemblies_datasets_genbank.edit.names \
        | sed 's/^/GCA_/g' > assemblies_datasets_genbank.uniq.names
    comm -23 assemblies_datasets_refseq.edit.names assemblies_datasets_genbank.edit.names \
        | sed 's/^/GCF_/g' > assemblies_datasets_refseq.uniq.names
    comm -12 assemblies_datasets_refseq.edit.names assemblies_datasets_genbank.edit.names \
        | sed 's/^/GCF_/g' > assemblies_datasets_refseq.union.names
    cat assemblies_datasets_refseq.union.names \
        assemblies_datasets_refseq.uniq.names \
        assemblies_datasets_genbank.uniq.names \
        >  assemblies_datasets_uniq.names
    
    n_uniq=$(cat assemblies_datasets_uniq.names | wc -l)
    echo "After removing redundancy: $n_uniq unique assemblies"
}

link_local_assemblies () {
    echo ""
    echo "#########################################"
    echo "#### Link assemblies from local DB ######"
    echo "#########################################"
    echo ""
    
    # Assemblies are already linked in assemblies_datasets_uniq/ from filter_by_taxon
    echo "✓ Using $(ls assemblies_datasets_uniq/*.gz 2>/dev/null | wc -l) assemblies from local database"
}

get_all_asm_list () {
    echo "###########################################"
    echo "##### Make a list of all assemblies  ######"
    echo "###########################################"
    :> all_asm_acc
    for I in $(ls assemblies_datasets_uniq/)
    do
        J=$(basename ${I%.*.*})
        echo ${J} >> all_asm_acc
    done
}

get_biosample_GEOdata () {
    echo "################################################"
    echo "##### mapping biosample to geo locations  ######"
    echo "###########   or Isolation country   ###########"
    echo "################################################"
    
    # Check for pre-downloaded biosample data
    biosample_file="$local_db_dir/metadata/All-$taxon-info.biosample.xml"
    
    if [ -f "$biosample_file" ]; then
        echo "Using pre-downloaded biosample metadata: $biosample_file"
        cp "$biosample_file" All-$taxon-info.biosample.xml
    else
        echo "Downloading biosample metadata from NCBI..."
        esearch -db biosample -query "txid${taxon}[organism]" | esummary \
            > All-$taxon-info.biosample.xml
    fi
    
    cat All-$taxon-info.biosample.xml |\
        xtract -pattern DocumentSummary -def "NA" \
            -element Accession  \
            -def "NA" \
            -block Attribute  \
            -if Attribute@harmonized_name \
            -equals geo_loc_name \
            -element Attribute | \
            sed 's/ /_/g' | \
            sed 's/Missing.*$/Unknown/g'|\
            sed 's/missing.*$/Unknown/g' |\
            sed 's/unknown.*$/Unknown/g' |\
            sed 's/not_collected.*$/Unknown/g' |\
            sed 's/not_applicable.*$/Unknown/g' |\
            sed 's/NONE.*$/Unknown/g' |\
            awk  '{if ($2 == "") print $1,"NA" ;else print $0}' | \
            sed 's/ /\t/g' | \
            sed 's/:/\t/g' | \
            awk '{print $1"\t"$2}'  \
            > All-$taxon.biosample.BS_to_Geoloc
}

merge_metadata_geoloc () {
    echo "################################################"
    echo "#####  add Geoloc data to metadata file   ######"
    echo "################################################"
    :> All-$taxon.BS_to_all_meta
    cat All-$taxon.assembly.BS_to_meta |
    while read BS BLAH
    do
        cat All-$taxon.biosample.BS_to_Geoloc |
        while read BS1 BLAH1
        do
            if [ $BS = $BS1 ]
            then
                echo -e $BS'\t'$BLAH'\t'$BLAH1
            fi
        done
    done >> All-$taxon.BS_to_all_meta
}

all_sample_metadata () {
    echo "################################################"
    echo "##### Filter all_metadata for assemblies  ######"
    echo "################################################"
    cat all_asm_acc | grep GCF > all_asm_acc.GCF
    cat all_asm_acc | grep GCA > all_asm_acc.GCA

    :> all_asm_acc_metadata
    cat All-$taxon.BS_to_all_meta |\
    grep -f all_asm_acc.GCF |\
    awk '{print $2,$4,$4"."$5,$1,$8,$9,$10,$11}' \
    >> all_asm_acc_metadata

    cat All-$taxon.BS_to_all_meta |\
    grep -f all_asm_acc.GCA |\
    awk '{print $3,$4,$4"."$5,$1,$8,$9,$10,$11}' \
    >> all_asm_acc_metadata
}

aggregate_assemblies () {
    echo "#####################################"
    echo "##### Create directory with all #####"
    echo "####### assemblies gunzip'd  ########"
    echo "#####################################"
    cd $wd
    out_dir=$1
    mkdir $out_dir
    
    for I in $(ls $wd/assemblies_datasets_uniq/*.gz)
    do
        base=$(basename ${I%.*.*})
        if [[ $I == *.gz ]]; then
            cat $I | gunzip > $out_dir/${base}.fna
        else
            cp $I $out_dir/${base}.fna
        fi
    done
}

get_stats_with_checkM () {
    echo "############################################################"
    echo "###### Run CheckM2 to get completeness, contamination ######"
    echo "####### and general assembly  metrics for filtering  #######"
    echo "############################################################"
    # CheckM2 (DIAMOND + pretrained ML models) replaces legacy CheckM1's
    #   pplacer/reference-tree placement, which OOM-killed even with
    #   --reduced_tree. No per-genome memory blow-up, so run the whole folder in
    #   a single `checkm2 predict` call. Kept in sync with
    #   utils/gather_filter_asms.sh:get_stats_with_checkM.
    threads=$threads
    checkM_input=$1
    checkM_dir=$2
    checkM_type=$3
    suffix=$4

    # $lowmem_flag may be unset in this script (no flag parser); default empty.
    lowmem_flag="${lowmem_flag:-}"
    # --force lets a re-run overwrite a non-empty output dir (CheckM2 aborts
    #   otherwise); the dir is recreated fresh each QC run anyway.
    if [[ $checkM_type == "protien" ]]
    then
        checkM_args="--threads $threads --genes -x $suffix --force $lowmem_flag"
    elif [[ $checkM_type == "genome" ]]
    then
        checkM_args="--threads $threads -x $suffix --force $lowmem_flag"
    else
        echo "Unknown checkM input type" && exit
    fi

    mkdir $checkM_dir
    cd $checkM_dir || exit

    echo "  Running CheckM2 predict with args: ${checkM_args}"
    checkm2 predict ${checkM_args} --input "$checkM_input" --output-directory "$checkM_dir"
    checkm_rc=$?
    if [ $checkm_rc -ne 0 ]
    then
        echo "ERROR: CheckM2 predict exited with code $checkm_rc." >&2
        echo "       If this is an out-of-memory kill, allocate more memory or" >&2
        echo "       add --lowmem. Aborting so unfiltered assemblies are NOT" >&2
        echo "       passed to OrthoPhyl." >&2
        exit 1
    fi
    quality_report="${checkM_dir}/quality_report.tsv"
    if [ ! -s "$quality_report" ]
    then
        echo "ERROR: CheckM2 produced no report output: missing/empty $quality_report" >&2
        echo "       CheckM2 likely crashed (possibly OOM). Aborting." >&2
        exit 1
    fi

    # Aggregate CheckM2 output into the legacy 18-column stats layout that
    #   filter_asm_by_stats consumes. duplication_ratio has no CheckM2 analogue
    #   (placeholder 0.00). Columns matched by HEADER NAME (order is not stable).
    echo "acc lineage #markerGenes #genomes_based_on missing 1copy 2copy 3copy 4copy 5+copy duplication_ratio completeness contamination GC GC_std Genome-size #scaffs scaff_N50" \
        > $wd/assemblies_all.stats.txt
    awk -F'\t' '
        NR==1 {
            for (i=1; i<=NF; i++) col[$i]=i
            next
        }
        {
            name=$(col["Name"]); comp=$(col["Completeness"]); cont=$(col["Contamination"])
            gc=$(col["GC_Content"]); size=$(col["Genome_Size"]); n50=$(col["Contig_N50"])
            nctg=(col["Total_Contigs"] ? $(col["Total_Contigs"]) : "NA")
            print name,"checkm2","NA","NA","NA","NA","NA","NA","NA","NA","0.00",comp,cont,gc,"NA",size,nctg,n50
        }
    ' "$quality_report" \
        >> $wd/assemblies_all.stats.txt
}

filter_asm_by_stats () {
    echo "###############################################"
    echo "#### filter out assemblies with low N50s, ####"
    echo "#### low total length and extra length  ######"
    echo "##############################################"
    cd $wd
    
    MIN_LEN=$1
    MAX_LEN=$2
    MIN_N50=$3
    MIN_GC=$4
    MAX_GC=$5
    MAX_dup=$6
    MIN_completeness=$7
    MAX_contam=$8
    
    # Same filtering logic as original script
    # (abbreviated here for space - copy from gather_filter_asms.sh)
    
    cat assemblies_all.stats.txt | \
    awk '{print $16,$18,$14,$11,$12,$13}' |\
    awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}
        END {for (i=1;i<=NF;i++) {
         printf "%f %f ", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' > assemblies_all.stats.ave_stddev

    while read scaf_bp scaf_bp_stddev scaf_N50 scaf_N50_stddev gc_avg gc_avg_stddev dup_ave dup_stddev complete_ave complete_stddev contam_ave contam_ave_stddev
    do
        if [ $MIN_LEN == "default" ]; then
            MIN_LEN=$(bc -l <<< "$scaf_bp-($scaf_bp_stddev*3)")
        fi
        if [ $MAX_LEN == "default" ]; then
            MAX_LEN=$(bc -l <<< "$scaf_bp+($scaf_bp_stddev*3)")
        fi
        if [ $MIN_N50 == "default" ]; then
            MIN_N50=$(bc -l <<< "$scaf_N50-($scaf_N50_stddev*3)")
        fi
        if [ $MIN_GC == "default" ]; then
            MIN_GC=$(bc -l <<< "$gc_avg-($gc_avg_stddev*3)")
        fi
        if [ $MAX_GC == "default" ]; then
            MAX_GC=$(bc -l <<< "$gc_avg+($gc_avg_stddev*3)")
        fi
        if [ $MAX_dup == "default" ]; then
            MAX_dup="0.02"
        fi
    done <<< $(cat assemblies_all.stats.ave_stddev)

    echo "Filtering assemblies based on:"
    echo "MIN_LEN=$MIN_LEN"
    echo "MAX_LEN=$MAX_LEN"
    echo "MIN_N50=$MIN_N50"
    echo "MIN_GC=$MIN_GC"
    echo "MAX_GC=$MAX_GC"
    echo "MAX_dup=$MAX_dup"
    echo "MIN_completeness=$MIN_completeness"
    echo "MAX_contam=$MAX_contam"

    :> assemblies_to_remove.stats
    cat assemblies_all.stats.txt |\
    tail -n +2 |\
    awk '{print $1,$16,$18,$14,$11,$12,$13}' |\
    while read acc scaf_bp ctg_N50 gc_avg dup complete contam
    do
        if (( $(echo "$scaf_bp < $MIN_LEN" |bc -l) )) || \
           (( $(echo "$scaf_bp > $MAX_LEN" |bc -l) )) || \
           (( $(echo "$ctg_N50 < $MIN_N50" |bc -l) )) || \
           (( $(echo "$gc_avg < $MIN_GC" |bc -l) )) || \
           (( $(echo "$gc_avg > $MAX_GC" |bc -l) )) || \
           (( $(echo "$dup > $MAX_dup" |bc -l) ))  || \
           (( $(echo "$complete < $MIN_completeness" |bc -l) )) || \
           (( $(echo "$contam > $MAX_contam" |bc -l) ))
        then
            echo ${acc} >> assemblies_to_remove.stats
        fi
    done
}

get_all_asms_to_remove () {
    echo ""
    echo "#######################################"
    echo "####### Get final list of ASMs ########"
    echo "#############  to remove  #############"
    echo "#######################################"
    echo ""
    cd $wd || exit
    cat assemblies_to_remove.* 2>/dev/null | sort | uniq > final_assemblies_to_remove
}

filter_for_redundancy () {
    echo ""
    echo "##############################"
    echo "####### Filter ASMs ##########"
    echo "####### For redundancy #######"
    echo "##############################"
    echo ""

    cd $wd || exit
    mkdir genomes_to_keep/
    ls assemblies_all.TMP > genome_list
    
    cat genome_list | sed 's/GC._//g' | sed 's/\..*fna//g' | sort | uniq | sed 's/^ *//g' > genome_list.accNum
    for I in $(cat genome_list.accNum)
    do 
        cat genome_list | grep $I | sort | tail -n 1 | sed 's/.fna//g'
    done | sort > genome_list.nunRedundant
    
    for I in $(comm -23 genome_list.nunRedundant final_assemblies_to_remove)
    do
        cp assemblies_all.TMP/$I.fna genomes_to_keep/
    done
    
    echo ""
    echo "✓ Final genome count: $(ls genomes_to_keep/*.fna | wc -l)"
    echo "✓ Output: $wd/genomes_to_keep/"
}

# Run main pipeline
main

echo ""
echo "======================================"
echo "✓ Local genome gathering complete"
echo "======================================"
echo "Output directory: $wd"
echo "Filtered genomes: $wd/genomes_to_keep/"
