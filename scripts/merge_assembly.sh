OUTFILE=$1
shift # Move all arguments to the left, so the outfile is skipped in the for loop
mkdir -p $(dirname $OUTFILE)
truncate -s 0 $OUTFILE
for arg in "$@"; do
    cat $arg >> $OUTFILE
done