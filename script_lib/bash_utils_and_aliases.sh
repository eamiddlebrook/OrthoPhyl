
# check if dir is absolute or relative to running dir
relative_absolute () {
if [[ $1 = ~/* ]] || [[ $1 = /* ]]
then
    	echo "$1"
else
    	echo $(pwd)/"$1"
fi
}

# print timeing for running functions
func_timing_start () {
	echo "${FUNCNAME[0]} started at $(date)" >> $store/timing
}
func_timing_end () {
	echo "${FUNCNAME[0]} ended at $(date)" >> $store/timing
}
