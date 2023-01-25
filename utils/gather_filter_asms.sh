#!/bin/bash
# USAGE:
# utils/gather_genomes.X.sh taxID_# output_dir_name


source ~/.bash_profile
#checkM and R are incompatable in conda, so this script needs its own conda ENV
echo "To create conda env use
conda create -n gather_genomes \
-c bioconda -c conda-forge \
checkm-genome bbmap entrez-direct ncbi-datasets-cli 
"

conda activate gather_genomes

#proxy-work

export taxon="$1"
export NCBI_assemblies=$2
wd=$(pwd)/$NCBI_assemblies
threads=2
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
MIN_completeness="98"
MAX_contam="1.0"

MAX_dup_default="0.02"
MIN_completness_default="98"
MAX_contam_default="0.2"

#taxonomy_check=FALSE

mkdir $NCBI_assemblies
cd $NCBI_assemblies
mkdir ./assemblies_datasets_uniq
echo "Output will be in $wd"

#######################################################
#### Declare main function for gather genomes pipe ####
#######################################################
# comment out function calls ass needed
main () {
	#get_NCBI_genomes
	#filter_NCBI_genomes
	#get_asm_metadata
	#get_non_datasets_assemblies
	#get_all_asm_list
	#get_biosample_GEOdata
	#merge_metadata_geoloc
	#all_sample_metadata
	#aggregate_assemblies "$wd"/assemblies_all.TMP
	#filter_asm_by_taxCheck
	#get_stats_with_checkM "$wd"/assemblies_all.TMP/ $wd/checkM_out genome fna
	###get_asm_stats # bbmap stats superseded by checkm
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
	datasets download genome taxon $taxon --dehydrated --exclude-genomic-cds --exclude-gff3 --exclude-protein --exclude-rna
	unzip ncbi_dataset.zip
	datasets rehydrate --gzip --directory ./
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

get_asm_metadata () {
        echo "################################################"
        echo "##### mapping biosample to GB accesstions ######"
        echo "################################################"
	esearch -db assembly -query "txid${taxon}[organism]" | esummary \
		> All-Brucella-info.assembly.xml
	cat All-Brucella-info.assembly.xml \
		| xtract -pattern DocumentSummary \
			-def "NA" -element BioSampleAccn RefSeq Genbank SpeciesName Sub_value FtpPath_GenBank FtpPath_RefSeq Taxid taxonomy-check-status ExclFromRefSeq | \
			sed 's/ /_/g' \
			> All-Brucella.assembly.BS_to_meta
}

get_non_datasets_assemblies () {
	echo "############################################"
        echo "###### Identify and DL all assemblies ######"
        echo "#### in NCBI asm DB but not in datasets ####"
        echo "############################################"

	#get info for assemblies not in datasets
	#  Uses All-Brucella.assembly.BS_to_meta to identify asm accessions
        cat All-Brucella.assembly.BS_to_meta | \
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
	esearch -db biosample -query "txid${taxon}[organism]" | esummary \
                > All-Brucella-info.biosample.xml
	#Tried to capture most cases of unknown value with sed. Super messy and dumb
	#   Rewrote to use the ATTR@subATTR syntax
	#   The "if" statement stuff is really dumb, cant figure out how to use "def" value
	#   still need sed cmds to change stuff to "Unknown"
        cat All-Brucella-info.biosample.xml |\
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
			> All-Brucella.biosample.BS_to_Geoloc
}

merge_metadata_geoloc () {
	echo "################################################"
	echo "#####  add Geoloc data to metadata file   ######"
	echo "################################################"
	# this is super duper dumb
	#   there is probably a builtin way to merge files by column falue
	#   oh well, files arnt huge
	:> All-Brucella.BS_to_all_meta
	cat All-Brucella.assembly.BS_to_meta |
	while read BS BLAH
	do
		cat All-Brucella.biosample.BS_to_Geoloc |
		while read BS1 BLAH1
		do
			if [ $BS = $BS1 ]
			then
				echo -e $BS'\t'$BLAH'\t'$BLAH1
			fi
		done
	done >> All-Brucella.BS_to_all_meta
}

all_sample_metadata () {
	echo "################################################"
        echo "##### Filter all_metadata for assemblies  ######"
	echo "################################################"
	cat all_asm_acc | grep GCF > all_asm_acc.GCF
        cat all_asm_acc | grep GCA > all_asm_acc.GCA

	:> all_asm_acc_metadata
	cat All-Brucella.BS_to_all_meta |\
	grep -f all_asm_acc.GCF |\
	awk '{print $2,$4,$4"."$5,$1,$8,$9,$10,$11}' \
	>> all_asm_acc_metadata

        cat All-Brucella.BS_to_all_meta |\
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
        echo "###### Run CheckM to get completeness, contamination  ######"
        echo "####### and general assembly  metrics for filtering  #######"
        echo "############################################################" 
	# checkM has the option to write stdout to file...
	#   dont know if this would make aggregating its output
	#   stdout is really messy
	max_genomes=200
	threads=$threads
	checkM_input=$1
        checkM_dir=$2
	checkM_type=$3
	suffix=$4
	if [[ $checkM_type == "protien" ]]
	then
		checkM_args="-t $threads -g -x $suffix"
	elif [[ $checkM_type == "genome" ]]
	then
		checkM_args="-t $threads -x $suffix"
	else
		echo "Unknown checkM input type" && exit
	fi
        mkdir $checkM_dir
        cd $checkM_dir || exit
        # split assemblies into different directories
        J=0
	K=0
	for I in $(ls $checkM_input/*.$suffix)
        do
          	if [ $((J % max_genomes)) -eq 0 ]
                then
                    	K=$((K+1))
                        mkdir $checkM_dir/${checkM_type}_${K}
                        mkdir $checkM_dir/${checkM_type}_${K}_out
                        cd $checkM_dir/${checkM_type}_${K}
                fi
                # make simlinks for assembly/proteome subsets
                ln -s $I ./
                J=$((J+1))
        done

        # run checkM on each proteome subset
        cd $checkM_dir
        J=1
	while [ $J -le $K ]
        do
          	in=$checkM_dir/${checkM_type}_${J}
                out=$checkM_dir/${checkM_type}_${J}_out
                checkm lineage_wf ${checkM_args} $in $out
                J=$((J+1))
        done

        # Aggregate checkM output
	echo "acc lineage #markerGenes #genomes_based_on missing 1copy 2copy 3copy 4copy 5+copy duplication_ratio completeness contamination GC GC_std Genome-size #scaffs scaff_N50" \
		> $wd/assemblies_all.stats.txt
	cat $checkM_dir/${checkM_type}_*_out/storage/bin_stats_ext.tsv | \
		sed 's/{,//g;s/,//g' | \
		awk '{print $1,$4,$10,$14,$16,$18,$20,$22,$24,$26,($26*5+$24*4+$22*3+$20*2)/$10,$28,$30,$32,$35,$38,$45,$57}' \
		>> $wd/assemblies_all.stats.txt
}

# CheckM gives all these stats, plus the completeness etc.
#   might be faster though, so if someone wants to only filter on simples stats
#   It could be useful...will save for now. 
#   Could format like CheckM output to use the same filtering function
get_asm_stats () {
	echo "###############################################"
        echo "###### get stats (BBmap) for assemblies  ######"
	echo "###############################################"
	cd $wd
	# get stats for assemblies.
	statswrapper.sh assemblies_additional/* > assemblies_additional.stats
	statswrapper.sh assemblies_datasets_uniq/* > assemblies_datasets.stats
	cat assemblies_additional.stats assemblies_datasets.stats |\
	grep -v n_scaffolds > assemblies_all.stats.txt

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
	#   This is also too complicated for bash...really slow
	#   SHOULD BE something like:
	#     more assemblies_all.stats.txt.old | grep Brucella | awk '{if ($11 <= 0.05) print $0}' | awk '{if ($12 >= 98) print $0}' | awk '{if ($13 <= 2.5) print $0}' | wc -l
	#     Must pass varable into awk explicitly (unfortunately)
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
