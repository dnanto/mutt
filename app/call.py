#!/usr/bin/env python3

import re
import sys
from argparse import ArgumentParser, FileType
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


def process_indels(aln):
	gaps = np.array([[nt == "-" for nt in rec.seq] for rec in aln], np.int, order="F").T
	# print(np.array([[nt == "-" for nt in rec.seq] for rec in aln], np.int, order="F"), file=sys.stderr)
	# print(aln, file=sys.stderr)
	# print(gaps, file=sys.stderr)
	result = "".join("01"[int(ele.all())] for ele in gaps[:-1] == gaps[1:])
	if len(result) > 0:
		result += result[-1]
	else:
		result = "0"
	# print(result, file=sys.stderr)
	# print(*re.finditer(r"(0|1+0?)", result), file=sys.stderr, sep="\n")
	for match in re.finditer(r"(0|1+0?)", result):
		start, end = match.start(0), match.end(0)
		ref = aln[0, start:end]
		is_ins = "-" in ref.seq
		for rec in aln[1:]:
			alt = rec[start:end]
			if alt.seq != ref.seq and ("-" in alt.seq or is_ins):
				yield start, ref, alt


def positions(ref):
	pos = 0
	coors = []

	for i, e in enumerate(ref):
		if e != "-":
			pos += 1
		coors.append(pos)

	return coors


def parse_args(argv):
	parser = ArgumentParser()
	parser.add_argument("path", type=FileType())
	parser.add_argument("-ref", "--ref")
	parser.add_argument("-out", "--out", type=Path, default=Path())
	args = parser.parse_args(argv)
	return args


def main(argv):
	args = parse_args(argv[1:])

	msa = AlignIO.read(args.path, "fasta", alphabet=IUPACAmbiguousDNA())

	if args.ref:
		msa = MultipleSeqAlignment(sorted(msa, key=lambda rec: rec.id != args.ref))

	chrom = str(msa[0].seq).replace("-", "")
	coors = positions(msa[0].seq)

	gaps = "".join(" -"["-" in msa[:, pos]] for pos in range(len(msa[0])))

	results = defaultdict(list)

	for match in re.finditer(r"--+", gaps):
		aln = msa[:, match.start(0):match.end(0)]
		# n = len({re.sub(r"[^-]+", " ", str(rec.seq)) for rec in aln})
		# if n > 2:
		for start, ref, alt in process_indels(aln):
			pos = coors[match.start(0) + start]
			if pos == 1:
				ini = chrom[pos + len(ref) - 1]
				rseq, aseq = (ref.seq + ini, ini) if "-" in alt.seq else (ini, alt.seq + ini)
			else:
				pos -= 1
				ini = chrom[pos - 1]
				rseq, aseq = (ini + ref.seq, ini) if "-" in alt.seq else (ini, ini + alt.seq)
			# print(ref.id, pos - 1, ".", rseq.upper(), aseq.upper(), ".", ".", ".", sep="\t")
			results[alt.id].append((pos, rseq.upper(), aseq.upper()))

	for rec in msa[1:]:
		for idx, ele in enumerate(zip(msa[0].seq, rec.seq)):
			ref, alt = ele
			if ref != alt and ref != "-" and alt != "-":
				# print(msa[0].id, coors[idx], ".", ref.upper(), alt.upper(), ".", ".", ".", sep="\t")
				results[rec.id].append((coors[idx], ref.upper(), alt.upper()))

	for key, val in results.items():
		path = args.out.joinpath(key + ".vcf")
		with path.open("w") as file:
			print("##fileformat=VCFv4.0", file=file)
			print("##contig=<ID=", msa[0].id, ",length=", len(msa[0]), ">", sep="", file=file)
			print("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", sep="\t", file=file)
			for ele in val:
				print(msa[0].id, ele[0], ".", ele[1], ele[2], ".", ".", ".", sep="\t", file=file)
		print(path)

	return 0


if __name__ == "__main__":
	sys.exit(main(sys.argv))
