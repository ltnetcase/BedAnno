#!/bin/bash

if [ $# -lt 3 ]; then
    echo "Usage: $0 <test.bed> <whole_db.bed.gz> <whole_trSeq.fa.gz>"
    exit 1
fi

BED=$1
DB=$2
TRSEQ=$3
bindir=`dirname $0`
tooldir="$bindir/../tools"
$tooldir/tabix -R $1 $2 | bgzip > test_db.bed.gz
$tooldir/tabix -p bed -f test_db.bed.gz
gzip -dc test_db.bed.gz | gawk '{for (i=4;i<=NF;i++) {split($i,itms,"|"); hash[itms[1]] = 1;}} END {for (tr in hash) {print tr;}}' | sort > trans.list
gzip -dc $3 | grep -A1 -f trans.list | grep -v -- "^--" | gzip > test.fas.gz
#samtools faidx test.fas.gz
