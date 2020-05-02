#!/bin/bash
FILE="./seqs.txt"
if [ ! -f "$FILE" ]
then
    touch $FILE
elif [ -s "$FILE" ]
then
    > "$FILE"
fi

while read -r CURRENT_LINE
    do
        /usr/local/ncbi/blast/bin/blastdbcmd -db ~/ColauttiLabScratch/16SMicrobialDB/16SMicrobial -entry $CURRENT_LINE | tee -a ./seqs.txt
        ((LINE++))
done < "./GBACC.txt"
