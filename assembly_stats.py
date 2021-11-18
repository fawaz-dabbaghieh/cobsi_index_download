import gzip
import os
import sys
import argparse
import matplotlib.pyplot as plt
from math import log


def read_fasta_gen(fasta_file_path):
	"""
	A generator function that reads one read at a time
	Can be used for big FASTA files to not keep them in memory

	:param fasta_file_path: path to fasta file
	:yield: a tuple of sequence id and sequence
	"""

	if not os.path.exists(fasta_file_path):
		print(f"file {fasta_file_path} does not exist")
		sys.exit()

	if fasta_file_path.endswith("gz"):
		fasta_file = gzip.open(fasta_file_path, "rt")
	else:
		fasta_file = open(fasta_file_path, "r")

	seqs = []
	seq_name = ""
	for line in fasta_file:
		line = line.strip()
		if not line:  # empty line
			continue

		if line.startswith(">"):
			if len(seqs) != 0:  # there was a sequence before
				yield seq_name, "".join(seqs)
				seq_name = line[1:]
				seqs = []
			else:
				seq_name = line[1:]
		else:
			seqs.append(line)

	# last sequence
	if seqs:
		yield seq_name, "".join(seqs)



parser = argparse.ArgumentParser(description='Some stats related to downloaded assemblies', add_help=True)
subparsers = parser.add_subparsers(help='Available subcommands', dest="subcommands")


creating_table = subparsers.add_parser('assemb_stats', help='Creating a table with minimum assemblies statistics')

creating_table.add_argument("--in_dir", dest="in_dir", type=str, default=None,
							help="Provide a directory that contains the assemblies in .fasta or .fasta.gz format")
creating_table.add_argument("--out_table", dest="out_table", type=str, default="assembly_info_table.tsv",
							help="The output table path/name. Default: assembly_info_table.tsv")

histograms = subparsers.add_parser('histograms', help='Plotting simple histograms using the table created in assemb_stats')
histograms.add_argument("--in_table", dest="in_table", type=str, default=None,
						help="Give the input table generated with assemb_stats to generate the histograms")
histograms.add_argument("--out_png", dest="out_png", type=str, default="dists.png",
						help="Output png file with contigs and assembly length distributions. Default: dists.png")

args = parser.parse_args()

if args.subcommands == "assemb_stats":

	if args.in_dir is None:
		print("You need to provide an input directory")
		sys.exit()

	else:
		if not os.path.exists(args.in_dir):
			print(f"The path {args.in_dir} does not exist")
			sys.exit()

	if os.path.exists(args.out_table):
		print(f"The file {args.out_table} already exists")
		sys.exit()

	assembly_files = []
	if not args.in_dir.endswith(os.sep):
		args.in_dir += os.sep

	for f in os.listdir(args.in_dir):
		if f.endswith(".fasta") or f.endswith(".fa") or f.endswith("gz"):
		  assembly_files.append(args.in_dir + f)

	with open(args.out_table, "w") as out_file:
		# writing the header
		out_file.write("\t".join(["file_name", "num_of_contigs", "sequence_len"]) + "\n")

		for f in assembly_files:
			n_contigs = 0
			seq_len = 0
			for seq_name, seq in read_fasta_gen(f):
				n_contigs += 1
				seq_len += len(seq)
			out_file.write("\t".join([f, str(n_contigs), str(seq_len)]) + "\n")


if args.subcommands == "histograms":
	if not os.path.exists(args.in_table):
		print(f"The table {args.in_table} does not exist")
		sys.exit()

	n_contigs = []
	seq_lens = []
	with open(args.in_table, "r") as in_file:
		next(in_file)  # skipping header
		for l in in_file:
			l = l.strip().split()
			n_contigs.append(int(l[1]))
			seq_lens.append(int(l[2]))


	fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
	fig.suptitle('Distributions for assemblies', fontsize=12)

	#  Sturge's rule for choosing number of bins for a histogram
	contigs_bins = int(1 + (3.22*log(len(n_contigs))))

	axs[0].hist(n_contigs, bins=contigs_bins)
	axs[0].set_title("Contigs distribution")
	axs[0].set_xlabel("# contigs")
	axs[0].set_ylabel("Frequency")

	seq_bins = int(1 + (3.22*log(len(seq_lens))))
	axs[1].hist(seq_lens, bins=seq_bins)
	axs[1].set_title("Assembly lengths distribution")
	axs[1].set_xlabel("Assembly lengths")
	axs[1].set_ylabel("Frequency")

	if os.path.exists(args.out_png):
		print(f"Path {args.out_png} already exists")
		sys.exit()

	plt.savefig(args.out_png, dpi=900)
