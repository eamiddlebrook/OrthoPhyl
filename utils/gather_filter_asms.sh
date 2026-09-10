#!/bin/bash
# USAGE:
# utils/gather_genomes.X.sh taxID_# PATH_TO/output_dir_name

#checkM2 needs its own conda ENV (DIAMOND/lightgbm deps conflict with the base)
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
conda create -n gather_genomes \
-c bioconda -c conda-forge \
checkm2 bbmap entrez-direct ncbi-datasets-cli
then download the CheckM2 DIAMOND database once with:
checkm2 database --download
"
	conda activate gather_genomes || exit
fi


export taxon="$1"
# deal with WD below
threads=$3

# NCBI Entrez query term used by the esearch metadata steps below.
#   The "txid" prefix is only valid for a numeric taxID; a taxon *name*
#   must be used bare (e.g. "Gloeotrichia[organism]"). The wrapper may pass
#   either form, so pick the right one here instead of hardcoding "txid".
if [[ "$taxon" =~ ^[0-9]+$ ]]; then
	export entrez_query="txid${taxon}[organism]"
else
	export entrez_query="${taxon}[organism]"
fi

# Low-memory / stats-engine flags.
#   CheckM2 replaces the old CheckM1 pplacer/reference-tree design (which
#   OOM-killed even with --reduced_tree) with DIAMOND + pretrained ML models.
#   Its --lowmem flag halves DIAMOND's RAM use at the cost of runtime. The
#   legacy --reduced_tree flag is accepted as an alias for --lowmem so old
#   invocations (and the wrapper) keep working.
lowmem_flag=""
use_bbmap=false

# Parse optional flags
for arg in "$@"; do
    if [[ "$arg" == "--lowmem" || "$arg" == "--reduced_tree" ]]; then
        lowmem_flag="--lowmem"
        echo "Low RAM mode enabled: passing --lowmem to CheckM2"
    elif [[ "$arg" == "--use-bbmap" ]]; then
        use_bbmap=true
        echo "Using bbmap statswrapper instead of CheckM2 for genome statistics"
    fi
done

## Convert wd to absolute path (handles relative paths, ~, and absolute paths)
if [[ "$2" = ~* ]]; then
    # Expand tilde
    wd="${2/#\~/$HOME}"
elif [[ "$2" = /* ]]; then
	# Already absolute path
    export wd="$2"
else
    # Relative path - make absolute
	export wd="$(cd "$(dirname "$2")" 2>/dev/null && pwd)/$(basename "$2")"
fi


echo "Working directory (absolute): $wd"
# Set these if you hav a good idea of the expected values
# if not, set to "default" and they will be set as the average+-(3*stddev)
# of course stddev is a bit weird because the values are almost surely not normal..
# NOTE: N50 is checkM's definition (length of contig, not number)
MIN_LEN="default"
MAX_LEN="default"
MIN_N50="default"
MIN_GC="default"
MAX_GC="default"
MAX_dup="default"
MIN_completeness="95"
MAX_contam="1.0"

MAX_dup_default="0.02"
MIN_completness_default="98"
MAX_contam_default="0.2"

#taxonomy_check=FALSE

mkdir $wd
cd $wd || (echo $wd "doesnt exist...exiting" ; exit)
mkdir ./assemblies_datasets_uniq
echo "Output will be in $wd"

# Pick a temp dir that is BOTH writable AND short, then point every temp-dir
#   var at it. Two container constraints collide here:
#     * Singularity's image /tmp is read-only, which breaks edirect's
#       nquire/mktemp ("mktemp: ... Read-only file system") -> needs writable.
#     * Tools that bind an AF_UNIX socket under $TMPDIR (edirect helpers, and
#       any multiprocessing.Manager) hit a hard 108-char kernel path limit. A
#       temp dir deep under $wd overflows it ("OSError: AF_UNIX path too long")
#       -> needs short.
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

#######################################################
#### Declare main function for gather genomes pipe ####
#######################################################
# comment out function calls ass needed
main () {
	get_NCBI_genomes
	filter_NCBI_genomes
	get_asm_metadata
	get_non_datasets_assemblies
	get_all_asm_list
	get_biosample_GEOdata
	merge_metadata_geoloc
	all_sample_metadata
	aggregate_assemblies "$wd"/assemblies_all.TMP
	#filter_asm_by_taxCheck
	
	# Choose stats method based on flag
	if [ "$use_bbmap" = true ]; then
		get_asm_stats "$wd"/assemblies_all.TMP/
	else
		get_stats_with_checkM "$wd"/assemblies_all.TMP/ $wd/checkM_out genome fna
	fi
	filter_asm_by_stats $MIN_LEN $MAX_LEN $MIN_N50 $MIN_GC $MAX_GC $MAX_dup $MIN_completeness $MAX_contam
	get_all_asms_to_remove
	filter_for_redundancy
}

get_NCBI_genomes () {
echo "
#########################################
####### Get taxon $taxon genomes ########
####### From NCBI using datasets ########
#########################################
"
	# this might be a little redundant if I can grab GCF/GCAs from assemblies DB
	#   However, DLing with assembly's ftp path using wget
	#   led to a few corrupt gz files (~6/150)
	#   while 0/880 DL'd with datasets were corrupt
	#   Plus datasets is natively multithreaded
	
	# Download with retry logic for NCBI server failures
	max_retries=3
	attempt=1
	success=false
	
	while [ $attempt -le $max_retries ] && [ "$success" = "false" ]; do
		echo "Download attempt $attempt of $max_retries..."
		
		# Clean up any partial downloads from previous attempts
		if [ $attempt -gt 1 ]; then
			echo "Cleaning up partial downloads from previous attempt..."
			rm -rf ncbi_dataset ncbi_dataset.zip README.md rehydrate.log
		fi
		
		# Download and dehydrate
		if datasets download genome taxon $taxon --dehydrated; then
			if unzip -o -q ncbi_dataset.zip; then
				# Rehydrate with error checking
				if datasets rehydrate --gzip --directory ./ 2>&1 | tee rehydrate.log; then
					# Check for 503 errors in output
					if grep -q "503 Service Unavailable" rehydrate.log; then
						unavailable_count=$(grep -c "503 Service Unavailable" rehydrate.log)
						echo "WARNING: $unavailable_count files unavailable due to NCBI server errors (503)"
						
						# Check if we got at least some files
						if [ -d "ncbi_dataset/data" ]; then
							retrieved_count=$(ls ncbi_dataset/data/ | grep -c GC)
							echo "Retrieved $retrieved_count assemblies (some failed)"
							
							if [ $retrieved_count -gt 0 ]; then
								if [ $attempt -lt $max_retries ]; then
									echo "Will retry to get missing files..."
								else
									echo "WARNING: Proceeding with $retrieved_count assemblies (some downloads failed)"
									echo "You may want to re-run this later to get complete dataset"
									success=true
								fi
							else
								echo "ERROR: No assemblies retrieved"
							fi
						else
							echo "ERROR: ncbi_dataset/data/ directory not created"
						fi
					else
						# No errors in rehydration
						echo "✓ Download completed successfully"
						success=true
					fi
				else
					echo "ERROR: Rehydration failed"
				fi
			else
				echo "ERROR: Failed to unzip ncbi_dataset.zip"
			fi
		else
			echo "ERROR: datasets download command failed"
		fi
		
		# If not successful and more retries available, wait and retry
		if [ "$success" = "false" ] && [ $attempt -lt $max_retries ]; then
			# Exponential backoff: 30s, 60s, 120s
			wait_time=$((30 * 2**(attempt-1)))
			echo "Waiting ${wait_time}s before retry (NCBI servers may be overloaded)..."
			sleep $wait_time
		fi
		
		attempt=$((attempt + 1))
	done
	
	if [ "$success" = "false" ]; then
		echo "ERROR: Failed to download genomes after $max_retries attempts"
		echo "NCBI servers may be experiencing issues. Try again later."
		exit 1
	fi
	
	# generate a file with all accessions grabbed by datasets
	ls ncbi_dataset/data/ | grep GC > assemblies_datasets.names
}


filter_NCBI_genomes () {
echo "
#################################################
####### Remove RefSeq/Genbank redundancy ########
####### Prefering to keep RefSeq entries ########
#################################################
"
	# make a list of nonredundant accessions (using GCF if both are present)
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

	# move all nonredundant assemblies to another dir (./assemblies_datasets_uniq/)
	#   errors when there are only unplaced scaffold files
	cat assemblies_datasets_uniq.names | \
	while read -r I ; do cat ncbi_dataset/data/${I}/*fna.gz \
		> ./assemblies_datasets_uniq/$I.fna.gz ; done
	#pretty messy. errors for all asseemblies without "unplaced*" in it (could fix - too lazy)
	#cat assemblies_datasets_uniq.names | \
	#while read -r I ; do cp ncbi_dataset/data/${I}/unplace* ./assemblies_datasets_uniq/$I.fna.gz ; done
}

# Run `esearch -db DB -query QUERY | esummary > OUT` with bounded retries.
#   NCBI's eutils endpoint drops connections intermittently (e.g.
#   "curl (56) SSL_ERROR_SYSCALL"), and edirect does not always recover on its
#   own. Retry a few times and verify the result actually contains records.
#   Returns 0 on success (OUT has >=1 DocumentSummary), 1 on persistent failure.
run_esearch_with_retry () {
	local db="$1"
	local query="$2"
	local out="$3"
	local max_retries=3
	local attempt=1

	while [ $attempt -le $max_retries ]; do
		echo "esearch attempt $attempt of $max_retries (db=$db, query=$query)..."
		if esearch -db "$db" -query "$query" | esummary > "$out" 2>/dev/null \
			&& grep -q "<DocumentSummary" "$out"; then
			echo "✓ esearch on $db returned records"
			return 0
		fi
		echo "WARNING: esearch on $db failed or returned no records (attempt $attempt)"
		attempt=$((attempt + 1))
	done

	echo "ERROR: esearch on $db failed after $max_retries attempts (query=$query)"
	return 1
}

# need to fix the problem with "Sub_value" having multiple values screwing up the column numbers
get_asm_metadata () {
        echo "################################################"
        echo "##### mapping biosample to GB accesstions ######"
        echo "################################################"
	if ! run_esearch_with_retry assembly "$entrez_query" All-$taxon-info.assembly.xml; then
		# datasets (get_NCBI_genomes) already fetched the bulk assemblies;
		# this metadata only drives supplementary (non-datasets) downloads.
		# Degrade gracefully with an empty table so the downstream grep -vf
		# in get_non_datasets_assemblies has something to read.
		echo "WARNING: proceeding without assembly metadata; supplementary (non-datasets) assemblies will be skipped."
		:> All-$taxon.assembly.BS_to_meta
		return 0
	fi
	cat All-$taxon-info.assembly.xml \
		| xtract -pattern DocumentSummary \
			-def "NA" -element BioSampleAccn RefSeq Genbank SpeciesName Sub_value FtpPath_GenBank FtpPath_RefSeq Taxid taxonomy-check-status ExclFromRefSeq | \
			sed 's/ /_/g' \
			> All-$taxon.assembly.BS_to_meta
}

get_non_datasets_assemblies () {
	echo "############################################"
	echo "###### Identify and DL all assemblies ######"
	echo "#### in NCBI asm DB but not in datasets ####"
	echo "############################################"

	#get info for assemblies not in datasets
	#  Uses All-$taxon.assembly.BS_to_meta to identify asm accessions
	cat All-$taxon.assembly.BS_to_meta | \
	grep -vf assemblies_datasets_uniq.names \
	> assemblies_not_in_datasets.BS_to_meta

	mkdir assemblies_additional
	cd assemblies_additional
	cat ../assemblies_not_in_datasets.BS_to_meta |
	sed 's/ /_/g' |
	while read BS GCF GCA species strain GB_path RS_path TaxID
	do
		if [ $RS_path != "NA" ]
		then
			RS_file="${RS_path##*/}_genomic.fna.gz"
                        RS_path_full="${RS_path}/${RS_file}"
            wget "$RS_path_full"
			mv $RS_file $GCF.fna.gz
			# if gz file corrupt, try again
			if gunzip -t $GCF.fna.gz
			then
				echo "$GCF.fna.gz DL'd without error"
			else
				wget "$RS_path_full"
				mv $RS_file $GCF.fna.gz
			fi
		elif [ $GB_path != "NA" ]
		then
			GB_file="${GB_path##*/}_genomic.fna.gz"
			GB_path_full="${GB_path}/${GB_file}"
			wget "$GB_path_full"
			mv $GB_file $GCA.fna.gz
			# if gz file corrupt, try again
			if gunzip -t $GCA.fna.gz
			then
					echo "$GCA.fna.gz DL'd without error"
			else
					wget "$GB_path_full"
					mv $GB_file $GCA.fna.gz
			fi

		else
			echo "No Path for $GCF or $GCA assembly of $BS"
		fi
	done
	cd $wd
}

get_all_asm_list () {
	echo "###########################################"
        echo "##### Make a list of all assemblies  ######"
        echo "###########################################"
        :> all_asm_acc
        for I in $(ls assemblies_additional/)
        do
                J=$(basename ${I%.*.*})
                echo ${J} >> all_asm_acc
        done
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
	# grab XML file with all brucalla info from BioSample DB
	if ! run_esearch_with_retry biosample "$entrez_query" All-$taxon-info.biosample.xml; then
		# Geolocation is decorative metadata; degrade gracefully to an empty
		# table so merge_metadata_geoloc still has a file to join against.
		echo "WARNING: proceeding without biosample geolocation metadata."
		:> All-$taxon.biosample.BS_to_Geoloc
		return 0
	fi
	#Tried to capture most cases of unknown value with sed. Super messy and dumb
	#   Rewrote to use the ATTR@subATTR syntax
	#   The "if" statement stuff is really dumb, cant figure out how to use "def" value
	#   still need sed cmds to change stuff to "Unknown"
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
	# this is super duper dumb
	#   there is probably a builtin way to merge files by column falue
	#   oh well, files arnt huge
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
	for I in $(ls $wd/assemblies_additional/*.gz)
	do
                base=$(basename	${I%.*.*})
                cat $I | gunzip > $out_dir/${base}.fna
	done
	for I in $(ls $wd/assemblies_datasets_uniq/*.gz)
        do
                base=$(basename	${I%.*.*})
                cat $I | gunzip > $out_dir/${base}.fna
        done
}

filter_asm_by_taxCheck () {
        echo "#####################################################"
        echo "###### filter out assemblies with inconclusive ######"
        echo "##### taxonomy checks or other metadata values  #####"
        echo "#####################################################"
	cat all_asm_acc_metadata |\
		awk '{if ($6 != "OK") print $1}' \
		> assemblies_to_remove.taxCheck
}

get_stats_with_checkM () {
	echo "############################################################"
        echo "###### Run CheckM2 to get completeness, contamination ######"
        echo "####### and general assembly  metrics for filtering  #######"
        echo "############################################################"
	# CheckM2 (DIAMOND + pretrained ML models) replaces legacy CheckM1's
	#   pplacer/reference-tree placement, which OOM-killed even with
	#   --reduced_tree. There is no per-genome memory blow-up, so we run the
	#   whole assembly folder in a single `checkm2 predict` call instead of
	#   splitting into 200-genome pplacer batches.
	threads=$threads
	checkM_input=$1
    checkM_dir=$2
	checkM_type=$3
	suffix=$4

	# CheckM2 takes a folder + extension. Predicted-protein input uses --genes
	#   (the gather scripts only ever pass "genome", but keep the branch).
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

	if [[ -n "$lowmem_flag" ]]; then
		echo "  Using CheckM2 --lowmem (low RAM mode)"
	fi
	mkdir $checkM_dir
	cd $checkM_dir || exit

	echo "  Running CheckM2 predict with args: ${checkM_args}"
	checkm2 predict ${checkM_args} --input "$checkM_input" --output-directory "$checkM_dir"
	checkm_rc=$?
	# CheckM2 can still be OOM-killed (SIGKILL -> exit 137) or fail silently. If
	# we do not abort here, aggregation below produces a header-only stats file,
	# the stats filter matches nothing, and EVERY raw assembly is passed
	# downstream unfiltered. Fail loudly instead.
	if [ $checkm_rc -ne 0 ]
	then
		echo "ERROR: CheckM2 predict exited with code $checkm_rc." >&2
		echo "       If this is an out-of-memory kill, re-run with --lowmem (low RAM" >&2
		echo "       mode) or --use-bbmap, or allocate more memory. Aborting so" >&2
		echo "       unfiltered assemblies are NOT passed to OrthoPhyl." >&2
		exit 1
	fi
	# Even on a 0 exit code, verify CheckM2 actually wrote its report table.
	quality_report="${checkM_dir}/quality_report.tsv"
	if [ ! -s "$quality_report" ]
	then
		echo "ERROR: CheckM2 produced no report output:" >&2
		echo "       missing or empty $quality_report" >&2
		echo "       CheckM2 likely crashed (possibly OOM). Aborting so unfiltered" >&2
		echo "       assemblies are NOT passed to OrthoPhyl." >&2
		exit 1
	fi

	# Aggregate CheckM2 output into the legacy 18-column stats layout that
	#   filter_asm_by_stats consumes (it reads: acc[1], dup[11], completeness[12],
	#   contamination[13], GC[14], Genome-size[16], scaff_N50[18]).
	#   CheckM2 has no marker-copy duplication metric, so duplication_ratio is a
	#   placeholder (0.00), exactly like the bbmap path. Columns are matched by
	#   HEADER NAME because CheckM2's column order changes with mode/--genes.
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
			# acc lineage mark genomes missing 1c 2c 3c 4c 5c dup comp contam GC GCstd size scaffs N50
			print name,"checkm2","NA","NA","NA","NA","NA","NA","NA","NA","0.00",comp,cont,gc,"NA",size,nctg,n50
		}
	' "$quality_report" \
		>> $wd/assemblies_all.stats.txt

	# Sanity check: every input assembly must have a stats row. A mismatch means
	# CheckM2 silently dropped genomes (partial crash) -- abort rather than let the
	# stats filter under-report and pass unfiltered assemblies downstream.
	n_input=$(ls $checkM_input/*.$suffix 2>/dev/null | wc -l)
	n_stats=$(tail -n +2 $wd/assemblies_all.stats.txt | wc -l)
	echo "  CheckM2 stats: $n_stats rows for $n_input input assemblies"
	if [ "$n_stats" -eq 0 ]
	then
		echo "ERROR: CheckM2 produced no per-assembly stats ($wd/assemblies_all.stats.txt" >&2
		echo "       contains only a header). Aborting so unfiltered assemblies are NOT" >&2
		echo "       passed to OrthoPhyl." >&2
		exit 1
	fi
	if [ "$n_stats" -ne "$n_input" ]
	then
		echo "ERROR: CheckM2 stats row count ($n_stats) does not match the number of" >&2
		echo "       input assemblies ($n_input). CheckM2 likely crashed on some genomes." >&2
		echo "       Aborting so partially-filtered assemblies are NOT passed to OrthoPhyl." >&2
		exit 1
	fi
}

get_asm_stats () {
	echo "###############################################"
        echo "###### get stats (BBmap) for assemblies  ######"
	echo "###############################################"
	local input_dir=$1
	cd $wd

	echo "Running statswrapper.sh on assemblies..."
	# Run bbmap statswrapper on all assemblies
	statswrapper.sh ${input_dir}*.fna > assemblies_all.stats.bbmap.txt
	
	# Convert bbmap output to format expected by filter function
	# bbmap columns: #file n_scaffolds scaf_bp ... n50 ... gc_avg
	# Expected format matches checkM: acc lineage ... duplication_ratio completeness contamination GC GC_std Genome-size #scaffs scaff_N50
	# For bbmap: set dup=0, completeness=100, contamination=0 as placeholders
	echo "Converting bbmap stats to expected format..."
	echo "acc lineage #markerGenes #genomes_based_on missing 1copy 2copy 3copy 4copy 5+copy duplication_ratio completeness contamination GC GC_std Genome-size #scaffs scaff_N50" \
		> $wd/assemblies_all.stats.txt
	
	# Parse bbmap output and reformat
	# bbmap columns: n_scaffolds(1) scaf_bp(3) scaf_N50(6) gc_avg(18) filename(20)
	# Target CheckM format: acc(1) ... GC(14) GC_std(15) Genome-size(16) #scaffs(17) scaff_N50(18)
	cat assemblies_all.stats.bbmap.txt | \
		grep -v "^#" | \
		awk 'NR>1 {print $20,$1,$3,$6,$18}' | \
		sed 's|.*/||; s/.fna//g' | \
		awk '{print $1,"bbmap","NA","NA","NA","NA","NA","NA","NA","NA","0.00","100","0",$5,"NA",$3,$2,$4}' \
		>> $wd/assemblies_all.stats.txt
	
	echo "Stats collection complete using bbmap"
	echo "Note: Completeness, contamination, and duplication set to placeholders (100, 0, 0)"
	echo "      Filtering will only use: genome size, N50, and GC content"
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
	if [ $MAX_contam == '' ]
        then
            	echo "!!!!! Not enough args given...needs 8 !!!!!!!!"
                echo "please look at function filter_asm_by_stats for details"
                echo "Set variables....
                MIN_LEN=$1
                MAX_LEN=$2
                MIN_N50=$3
                MIN_GC=$4
                MAX_GC=$5
                MAX_dup=$6
                MIN_completeness=$7
                MAX_contam=$8
                "
                exit
        fi
	echo "USER FILTERING input"
	echo "MIN_LEN=$MIN_LEN"
        echo "MAX_LEN=$MAX_LEN"
        echo "MIN_N50=$MIN_N50"
        echo "MIN_GC=$MIN_GC"
        echo "MAX_GC=$MAX_GC"
	echo "MAX_dup=$MAX_dup"
        echo "MIN_completeness=$MIN_completeness"
        echo "MAX_contam=$MAX_contam"

	#create a file with averages and stdDevs for columns
	#   in checkM_stats_aggregated
	#   scaf_bp  scaf_N50 gc_avg dup complete contam
	cat assemblies_all.stats.txt | \
	awk '{print $16,$18,$14,$11,$12,$13}' |\
	awk '{for(i=1;i<=NF;i++) {sum[i] += $i; sumsq[i] += ($i)^2}}
        	END {for (i=1;i<=NF;i++) {
         	printf "%f %f ", sum[i]/NR, sqrt((sumsq[i]-sum[i]^2/NR)/NR)}
         }' > assemblies_all.stats.ave_stddev

	# set filtering defaults if not set at start of script
	#   This is too complicated for bash....should have writen in python.
	while read scaf_bp scaf_bp_stddev scaf_N50 scaf_N50_stddev gc_avg gc_avg_stddev dup_ave dup_stddev complete_ave complete_stddev contam_ave contam_ave_stddev
	do
		# if [ "default" in variables ]
		if [ $MIN_LEN == "default" ]
		then
                        MIN_LEN=$(bc -l <<< "$scaf_bp-($scaf_bp_stddev*3)")
		fi
		if [ $MAX_LEN == "default" ]
                then
                        MAX_LEN=$(bc -l <<< "$scaf_bp+($scaf_bp_stddev*3)")
		fi
		if [ $MIN_N50 == "default" ]
		then
			MIN_N50=$(bc -l <<< "$scaf_N50-($scaf_N50_stddev*3)")
		fi
		if [ $MIN_GC == "default" ]
                then
			MIN_GC=$(bc -l <<< "$gc_avg-($gc_avg_stddev*3)")
                fi
                if [ $MAX_GC == "default" ]
                then
			MAX_GC=$(bc -l <<< "$gc_avg+($gc_avg_stddev*3)")
                fi
		if [ $MAX_dup == "default" ]
                then
                    	MAX_dup=$MAX_dup_default
                fi
		if [ $MIN_completeness == "default" ]
                then
                    	MAX_complete=$MAX_complete_default
                fi
		if [ $MAX_contam == "default" ]
                then
                    	MIN_contam=$MIN_contam_default
                fi
	done <<< $(cat assemblies_all.stats.ave_stddev)


	# Filter assemblies
	echo "Filtering assemblies based on:"
        echo "MIN_LEN=$MIN_LEN"
        echo "MAX_LEN=$MAX_LEN"
        echo "MIN_N50=$MIN_N50"
        echo "MIN_GC=$MIN_GC"
        echo "MAX_GC=$MAX_GC"
	echo "MAX_dup=$MAX_dup"
        echo "MIN_completeness=$MIN_completeness"
        echo "MAX_contam=$MAX_contam"
	echo "Paths to assemblies being filtered out are found in assemblies_to_remove.stats"
	# Emply output file
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
			#echo "$scaf_bp $ctg_N50 $gc_avg $dup $complete $contam"
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
	cat assemblies_to_remove.* | sort | uniq > final_assemblies_to_remove
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
	# this removes multiple versions of assemblies (i.e. GCA_XXXX.1 and GCA_XXXX.2)
	cat genome_list | sed 's/GC._//g' | sed 's/\..*fna//g' | sort | uniq | sed 's/^ *//g' > genome_list.accNum
	for I in $(cat genome_list.accNum)
	do 
		cat genome_list | grep $I | sort | tail -n 1 | sed 's/.fna//g'
	done | sort > genome_list.nunRedundant
	# iterate over acc that are in genome_list.nunRedundant but not final_assemblies_to_remove
	for I in $(comm -23 genome_list.nunRedundant final_assemblies_to_remove)
	do
		cp assemblies_all.TMP/$I.fna genomes_to_keep/

	done
}

# run main pipe
main


######################
##### notes ##########
#####################
# change file names
rename_asm_files () {
	while IFS= read -r line 
		do new_name=$(echo -e $line | awk '{print $3_$1}') 
		name=$(echo -e $line | awk '{print $1}') ; echo $new_name 
		echo $name
		cp Glutamicibacter_genomes.8.28.23/assemblies_all.TMP/$name* Paeniglutamicibacter_genomes4manuscript_tree/$new_name.fna
		done < Glutamicibacter_genomes.8.28.23/all_asm_acc_metadata
}
