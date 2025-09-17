export script_home=$(dirname "$(readlink -f "$0")")
cd $script_home
threads=$1
num_for_OrthoFinder=$2
outdir=./TESTER/FULLTEST_OUT.$(date +%m-%d-%Y)

#if [ -d ${outdir} ]
#then
#	mv ${outdir} ${outdir}.old
#fi

./OrthoPhyl.sh \
	-s ${outdir} \
	-a ./TESTER/annots_nucls_fasttest,./TESTER/annots_prots_fasttest \
	-o BOTH \
	-R full \
	-g ./TESTER/genomes_fasttest \
	-t $threads \
	-c control_file.user \
	-n $num_for_OrthoFinder

./ReLeaf.sh \
	-g ~/gits/OrthoPhyl/TESTER/genomes_fasttest_addasm \
	-a ./TESTER/annots_nucls_fasttest_addasm,./TESTER/annots_prots_fasttest_addasm \
	-s ~/gits/OrthoPhyl/${outdir} \
	-t $threads \
	-p iqtree \
	-o BOTH
