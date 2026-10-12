
# Pick a writable TMPDIR and export it for the rest of the run. Singularity's
#   image /tmp is often read-only at runtime, which breaks anything that
#   hardcodes /tmp or falls back to it (OrthoFinder's mafft dependency check,
#   mafft/hmmbuild/trimal's own temp-file handling, etc). /dev/shm is a
#   writable tmpfs even in that case; fall back through /tmp, $HOME, and
#   finally $store for non-container/unusual setups. Mirrors the identical
#   fix in utils/gather_filter_asms.sh for the download/QC phase.
# op_tmp is deliberately not `local` -- the deferred `trap` must still
#   resolve it in the calling script's scope when the shell exits.
SETUP_OP_TMPDIR () {
	op_tmp=""
	for base in /dev/shm /tmp "$HOME" "$store"; do
		if [ -d "$base" ] && [ -w "$base" ]; then
			op_tmp=$(mktemp -d "$base/op_tmp.XXXXXX" 2>/dev/null) && break
		fi
	done
	if [ -z "$op_tmp" ]; then
		echo "ERROR: could not create a writable temp dir (tried /dev/shm /tmp \$HOME \$store)" >&2
		exit 1
	fi
	export TMPDIR="$op_tmp"
	export TMP="$TMPDIR"
	export TEMP="$TMPDIR"
	trap 'rm -rf "$op_tmp"' EXIT
	echo "Using temp dir (writable, container-safe): $TMPDIR"
}

# check if dir is absolute or relative to running dir
relative_absolute () {
if [[ $1 = "/"* ]] || [[ $1 = "~/"* ]]
then
    	echo $1
else
    	echo $(pwd)/"$1"
fi
}

# print timeing for running functions
func_timing_start () {
	echo "${FUNCNAME[1]} started at $(date)" >> $store/timing
}
func_timing_end () {
	echo "${FUNCNAME[1]} ended at $(date)" >> $store/timing
}
