#!/usr/bin/env bash

msa="$1"
ref="$2"

root="$(dirname "$msa")"
out="$root/vcf"

rm -rf "$out" "$root/ref.fas.fai"
mkdir -p "$out"

echo "samtools..." > "$root/gmm.log"
samtools faidx "$root/rec.fas" "$ref" > "$root/ref.fas"

echo "calling variants..." > "$root/gmm.log"
./call.py "$msa" -ref "$ref" -out "$out" | \
while read -r ele
do
	echo "bcftools $ele..." > "$root/gmm.log"
	bcftools sort "$ele" | bcftools norm -c e -d all -f "$root/ref.fas" -o "$ele.tmp"
	mv "$ele.tmp" "$ele"
	echo "$ele"
done > "$root/gmm.txt" 2> /dev/null

rm "$root/gmm.log"

exit 0
