UNMAPPED=$1
ORIGINAL=$2

zcat $UNMAPPED | while read HEADER; read SEQUENCE; read SEPARATOR; read QUALITY;
do
    ## TODO comparison
    echo $HEADER;
    _ORIG=$(zgrep $HEADER $ORIGINAL -A3;)
    echo $HEADER;
    echo $_ORIG
    echo ""
done