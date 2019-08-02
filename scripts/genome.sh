#!/usr/bin/env bash

PATH="$( cd "$(dirname "$0")" || exit 1 ; pwd -P )":"$PATH"

taxon="$1"
db="$2"
slen1="$3"
slen2="$4"

mkdir -p "$taxon" && cd "$taxon" || exit 1

query="txid$taxon[PORGN] AND $slen1:$slen2[SLEN] AND biomol_genomic[PROP] NOT gbdiv_pat[PROP] NOT gbdiv_syn[PROP]"
echo "$query" > qry.txt

esearch -db "$db" -query "$query" > rec.tmp
efetch -format fasta < rec.tmp > rec.fas
efetch -format gb < rec.tmp > rec.gbk
efetch -format docsum < rec.tmp > rec.xml
rm -f rec.tmp

xtract \
	-input rec.xml -pattern DocumentSummary -def NA \
	-element AccessionVersion -element SubType -element SubName | \
	subtype.R 2> /dev/null | \
	lubridate.R - collection_date > rec.tsv 2> /dev/null
