OUTFILE=$1
shift # Move all arguments to the left, so the outfile is skipped in the for loop
mkdir -p $(dirname $OUTFILE)
truncate -s 0 $OUTFILE
for arg in "$@"; do
    file_name=$(basename $arg)
    if [[ "$file_name" = *".annotation.gtf" ]]; then
        id="${file_name%.*.*}"
        sed "s/gene_id \"/gene_id \"${id}./g" $arg >> $OUTFILE
    else
        cat $arg >> $OUTFILE
    fi
done