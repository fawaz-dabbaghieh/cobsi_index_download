# Downloading ENA assemblies from COBSI index
Generating an information table and downloading files from the COBSI index

A simple python script that takes the table from the ENA2018 database [here](http://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k/sampleid_assembly_paths.txt) and produces a table with each sample's accession number and assembly URL. Second functionality is that it takes this big table and downloads all the assemblies according to a species name of Taxon ID that the user gives it. It also support multi processing to make things faster.

## Generate table
First functionality is generating the big table.

Example:
```bash
python3 get_assemblies.py --cores 10 create_table --samples sampleid_assembly_paths.txt --output_table sampleid_assembly_table.tsv
```

This table will have this information `accession alias ena_accession broker_name taxon_id sample_name path_to_assembly`.


## Download assemblies
The subcommand `get_contigs` takes the table produced in the first step and takes either a taxon id or a genus or species name (e.g. bacillus, coli) and download all the assemblies corresponding to the chosen organism.

Example:
```bash
python3 create_info_tabl.py --cores 20 get_contigs --info_table sampleid_assembly_table.tsv --taxon_id 1392 --output_dir assemblies_id_1392/ --request_timeout 15 > log.txt

python3 create_info_tabl.py --cores 10 get_contigs --info_table sampleid_assembly_table.tsv --org_name bacillus --output_dir assemblies_bacillus/ > log.txt
```

The `request_timout` argument is to cancel the download in case there was a problem and the download request wasn't being processed. This way, it reports it in the log file but also doesn't keep the thread running and can move on to other ones to download. In case this happens, one simple way to get the missing assemblies is to run the script again, if nothing of the inputs changed, then it will skip the files that were already downloaded and download the missing one.


# Statistics on the downloaded assemblies
The script `assembly_stats.py` can help generate simple stats on the assemblies

## Generate stats table
By calling `python3 assembly_stats.py assemb_stats --in_dir assemblies/ --out_table assemblies_info.tsv` you get a TSV files with three columns, the assembly file path, number of contigs and size of sequence in that assembly.

## Generate simple histogram
By calling `python3 assembly_stats.py histograms --in_table assemblies_info.tsv --out_png assemblies_hists.png` It will produces a PNG files with two histograms, one with distribution of number of contigs and the other with distribution of size of assemblies downloaded.