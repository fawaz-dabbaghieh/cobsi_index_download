import pdb
import sys
import os
import requests
import argparse
import subprocess
import pdb
import xml.etree.ElementTree as ET
import multiprocessing as mp


ACCESSION = 0
ENA_ACCESSION = 2
ALIAS = 1
TAXON_ID = 4
SAMPLE_NAME = 5
URL = 6

def download_file(info_line, out_path, req_timeout):

	out_file_name = out_path + "_".join(["_".join(info_line[SAMPLE_NAME].split()), info_line[TAXON_ID], info_line[ENA_ACCESSION], info_line[ACCESSION]]).replace(os.sep, "_").replace(".", "_").replace("__", "_") + ".fa.gz"
	if not os.path.exists(out_file_name):

		print(f"Going to download {info_line[URL]} and put it in {out_file_name}")
		try:
			r = requests.get(info_line[URL], allow_redirects=True, timeout=req_timeout)
		except Exception as e:
			print(f"Error: {e} there seems to be an error with getting {info_line[URL]}")
			r = ""

		if r:
			with open(out_file_name, "wb") as out_file:
				out_file.write(r.content)
	else:
		print(f"File {out_file_name} already exists, skipping that entry")

	# urllib.urlretrieve(url, out_file_name)

def match_name(given_name, sample_name):
	# I match to any word in the sample name in the xml file
	# because sometimes there are extra information related to strain
	# so I don't know which index after splitting sample_name will match
	# this is more of a general match, so it will download anything with given_name in it
	for c in sample_name.split():
		if c == given_name:
			return True
	return False


def taxon_or_org(type, first, second):
	if type == "taxon":
		if first == second:
			return True
		return False

	elif type == "org_name":
		return match_name(first, second)


def get_xml(accession):
	xml = requests.get(f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}?download=true")
	if xml.status_code != 200:
		return False

	with open(f"tmp/{accession}", "wb") as out_file:
		out_file.write(xml.content)
	return True


def parse_xml(accession, path, queue):

	info_dict = {"path":path, "ena_accession":accession}

	if not get_xml(accession):
		queue.put(info_dict)
		return None


	try:
		file_tree = ET.parse(f"tmp/{accession}")
		tree_root = file_tree.getroot()

		for root_item in tree_root:
			if root_item.tag == "SAMPLE":
				sample = root_item
				# the SAMPLE tag attributes has accession, alias and broker_name information

				for key, value in sample.attrib.items():
					info_dict[key] = value

				for sample_element in sample:
					if sample_element.tag == 'SAMPLE_NAME':
						for element in sample_element:
							if element.tag == "TAXON_ID":
								info_dict["taxon_id"] = int(element.text)
							if element.tag == "SCIENTIFIC_NAME":
								info_dict["scientific_name"] = element.text.lower()
		queue.put(info_dict)

	except FileNotFoundError:
		queue.put(info_dict)


def return_line(info_dict):
	keys = ["accession", "alias","ena_accession", "broker_name", "taxon_id", "scientific_name", "path"]
	output = []
	# output = [info_dict["accession"], info_dict["alias"], info_dict["broker_name"], info_dict["taxon_id"], info_dict["scientific_name"], info_dict["path"]]
	for k in keys:
		if k in info_dict:
			output.append(info_dict[k])
		else:
			output.append("NA")

	return("\t".join([str(x) for x in output]))


def yeild_accessions(sample_id_cobsi_table):
	
	with open(sample_id_cobsi_table, "r") as in_file:
		for l in in_file:
			yield(l.strip().split())


parser = argparse.ArgumentParser(description='Get COBSI database info for an organism', add_help=True)
subparsers = parser.add_subparsers(help='Available subcommands', dest="subcommands")


parser.add_argument("--cores", dest="cores", type=int, default=1,
					help="You can specify corse to make the process faster")


creating_table = subparsers.add_parser('create_table', help='Creating a table with information from accessions from COBSI database')
creating_table.add_argument("--samples", dest="sample_id", type=str, default=None,
					help="give the path to the table from COBSI database with accession and assemblies paths")
creating_table.add_argument("--output_table", dest="output_table", type=str, default="cobsi_sample_information.tsv",
							help="The output path for the table")


get_contigs = subparsers.add_parser("get_contigs", help="Get assemblies related to taxon id or organism given")
get_contigs.add_argument("--info_table", dest="table", type=str, default=None,
						help="The information table that was produce with crate_table subcommand")

get_contigs.add_argument("--taxon_id", dest="taxon_id", type=int, default=None,
						help="Give the taxon id of the organism you are interested in")

get_contigs.add_argument("--org_name", dest="org_name", type=str, default=None,
						help="Give the genus or species name, e.g. Myxoccoccus, Escherichia")

get_contigs.add_argument("--output_dir", dest="out_dir", type=str, default=".",
						help="Speicfy the output directory for the assemblies")

get_contigs.add_argument("--request_timeout", dest="req_timeout", type=int, default=20,
						help="How long before timing out the file download request, so the script can move to next assembly. Default: 20s")

args = parser.parse_args()

if args.cores > os.cpu_count():
	print("You don't have that many corse unfortunately")
	sys.exit()



if args.subcommands == "create_table":
	if args.sample_id is None:
		print("You need to provide a path to the COBSI file with accessions and path to assemblies")
		sys.exit()

	try:
		print("making tmp directory")
		os.mkdir("tmp")
	except:
		subprocess.run("rm tmp/*", shell=True)

	samples = []
	for accession, path in yeild_accessions(args.sample_id):
		samples.append((accession, path))


	out_file = open(args.output_table, "w")
	print(f"Going to build the table and put it in {out_file}")
	# writing header
	out_file.write("\t".join(["accession", "alias", "ena_accession" , "broker_name", "taxon_id", "sample_name", "path_to_assembly"]) + "\n")
	counter = 0
	processes = []
	checkpoint = int(len(samples)/10)
	if checkpoint == 0:  # avoid divide by 0
		checkpoint = 1

	queue = mp.Queue()
	processes = []

	print(f"Going to loop through the samples and get their xml information to build the table")
	for accession, path in samples:
		# fixing the path to be a valid ftp path
		new_path = path.split("/")
		new_path = "http://ftp.ebi.ac.uk/" + "/".join(new_path[3:])

		process = mp.Process(target=parse_xml, args=(accession, new_path, queue,))
		processes.append(process)
		counter += 1
		if len(processes) == args.cores:
			for p in processes:
				p.start()
			for p in processes:
				p.join()
			for p in processes:
				out_file.write(return_line(queue.get()) + "\n")
			# emptying to prepare the next batch of graphs
			processes = []
			queue = mp.Queue()
			subprocess.run("rm tmp/*", shell=True)

		if counter % checkpoint == 0:
			print(f"So far {counter} accessions have been processed")

	# processing the leftovers
	if processes:
		for p in processes:
			p.start()
		for p in processes:
			p.join()
		for p in processes:
			out_file.write(return_line(queue.get()) + "\n")

	subprocess.run("rm -r tmp/", shell=True)

	out_file.close()
	print(f"Done!")


if args.subcommands == "get_contigs":
	if args.table is None:
		print("You need to provide the path to the table produced with create_table subcommand")
		sys.exit()

	if (args.taxon_id is None) and (args.org_name is None):
		print("You need to either give a taxon id as integer or an organism name")
		sys.exit()

	elif (args.taxon_id is not None) and (args.org_name is not None):
		print("You need to either give a taxon id or an organism, not both")
		sys,exit()


	if args.taxon_id is not None:
		print(f"Will search for taxon id {args.taxon_id} and download the assemblies related to it in {args.out_dir}")
		in_type = "taxon"
	elif args.org_name is not None:
		print(f"Will search for organism {args.org_name} and download the assemblies related to it in {args.out_dir}")
		in_type = "org_name"

	processes = []

	print(f"Looping through the table to check which samples match the desired species chosen")
	with open(args.table, "r") as in_file:
		next(in_file)  # skipping header
		for l in in_file:
			l = l.strip().split("\t")
			if in_type == "taxon":
				first = l[TAXON_ID]
				second = str(args.taxon_id)
			else:
				first = args.org_name
				second = l[SAMPLE_NAME].lower()

			if taxon_or_org(in_type, first, second) is True:
				process = mp.Process(target=download_file, args=(l, args.out_dir, args.req_timeout,))
				processes.append(process)
				if len(processes) == args.cores:
					for p in processes:
						p.start()
					for p in processes:
						p.join()

					processes = []

		print(f"Finished most of the processes and checking if there are any leftovers")
		if processes:
			for p in processes:
				p.start()
			for p in processes:
				p.join()

		print(f"Done!")

	# elif args.org_name is not None:
	# 	type = "org_name"

	# 	processes = []

	# 	with open(args.table, "r") as in_file:
	# 		next(in_file)
	# 		for l in in_file:
	# 			l = l.strip().split("\t")
	# 			if match_name(args.org_name, l[SAMPLE_NAME].lower()):

	# 				process = mp.Process(target=download_file, args=(l, args.out_dir))
	# 				processes.append(process)
	# 				# pdb.set_trace()
	# 				if len(processes) == args.cores:
	# 					for p in processes:
	# 						p.start()
	# 					for p in processes:
	# 						p.join()
	# 					processes = []

	# 		if processes:
	# 			for p in processes:
	# 				p.start()
	# 			for p in processes:
	# 				p.join()
