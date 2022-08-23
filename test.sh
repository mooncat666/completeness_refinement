while [ -n "$1" ] ; do
case "$1" in
-d) outdir=$2 
final=$2 ;;
-s | --species) species=$2 ;;
--temp) temp_outdir=$2 ;;
-v | --verbose) _verbose=True;;
-h) echo "Usage: bash -i pipe.sh -s species [-d output_directory] [--temp temp_dir] [-h]"
exit 1
shift ;;
esac
shift
done


if [[ -n ${_verbose} ]];
then
function log () {
echo "$@"
}
else
function log () { 
:
}
fi

echo start
log There is the log
echo end
