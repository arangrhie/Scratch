import sys
import re
import os
import argparse


parser = argparse.ArgumentParser(description="Convert GFF3 file to geneType, colors, and geneNames files")
parser.add_argument("gff3_file", help="Path to the GFF3 file to be converted")
parser.add_argument("-useName", action="store_true", help="Use Name instead of ID")
args = parser.parse_args()

gff3_file = args.gff3_file
use_name = args.useName

if not os.path.isfile(gff3_file):
	print(f"Error: {gff3_file} does not exist")
	sys.exit(1)

# Define the output file paths
base_filename = os.path.splitext(gff3_file)[0]
gene_type_file = base_filename + ".geneType.txt"
colors_file = base_filename + ".colors.txt"
gene_names_file = base_filename + ".geneNames.txt"

# Color codes
col_protein_coding = "25\t25\t112" # midnight blue
col_pseudogene = "65\t105\t225" # royal blue
col_rna = "50\t205\t50" # lime green
col_else = "128\t128\t128" # gray

with open(colors_file, 'w') as color_file:
    pass

# Function to take gene_biotype as input, write output to colors.txt
def write_rgb_codes(id, gene_biotype):

	# Open the colors.txt file for writing
	with open(colors_file, "a") as color_file:
		
		# Check the gene_biotype and write the corresponding RGB code
		if gene_biotype == "protein_coding":
			color_file.write(f"{id}\t{col_protein_coding}\n")
		elif gene_biotype.endswith("pseudogene"):
			color_file.write(f"{id}\t{col_pseudogene}\n")
		elif gene_biotype.endswith("RNA"):
			color_file.write(f"{id}\t{col_rna}\n")
		else: # gray
			color_file.write(f"{id}\t{col_else}\n")
			SystemError(f"Unknown gene_biotype: {gene_biotype=}")

# Open the input GFF3 file for reading
with open(gff3_file, "r") as file, open(gene_type_file, "w") as gene_type_out, open(gene_names_file, "w") as gene_names_out:
	# Iterate over each line in the input file
	for line in file:
		# Skip comment lines
		if line.startswith("#"):
			continue
		# Split the line into fields
		fields = line.strip().split("\t")
		
		# Extract the feature type
		feature_type = fields[2]
		
		# Extract the attributes field
		attributes = fields[8]
		
		# Extract gene_biotype, gene_name, and gene_synonym from gene feature
		if feature_type.endswith("gene"):

			gene_biotype_match = re.search(r"gene_biotype=([^;]+)", attributes)
			if gene_biotype_match:
				gene_biotype = gene_biotype_match.group(1)
			else:
				gene_biotype = ""

			if (use_name):
				# RefSeq Hack
				id_match = re.search(r"Name=([^;]+)", attributes)
				gene_name_match = id_match
			else:
				# Liftoff
				id_match = re.search(r"ID=([^;]+)", attributes)
				gene_name_match = re.search(r"gene_name=([^;]+)", attributes)
			if id_match:
				id = id_match.group(1)
			else:
				# No Name field found, use gene field
				id_match = re.search(r"gene=([^;]+)", attributes)
				id = id_match.group(1)
			
			gene_type_out.write(f"{id}\t{gene_biotype}\n")
			write_rgb_codes(id, gene_biotype)

			# Extract gene name
			if gene_name_match:
				gene_name = gene_name_match.group(1)
			else:
				gene_name = ""

			# Extract gene synonyms from gene feature
			# Duplicate gene_name to match the 3 column requirements
			# if not, genePredToBigGenePred gets confused
			gene_synonyms_match = re.search(r"gene_synonym=([^;]+)", attributes)
			if gene_synonyms_match:
				gene_synonyms = gene_synonyms_match.group(1)
				gene_names_out.write(f"{id}\t{gene_name}\t{gene_synonyms}\n")
			else:
				gene_synonyms = ""
				gene_names_out.write(f"{id}\t{gene_name}\t{gene_name}\n")
			
			if (use_name):
				parent_name = id
			else:
				parent_name = gene_name

		# Extract ID from transcript feature
		# Duplicate gene_name to match the 3 column requirements
		# if not, genePredToBigGenePred gets confused
		elif feature_type.endswith("transcript") or feature_type.endswith("RNA") or feature_type.endswith("gene_segment") or feature_type == "enhancer": 
			if (use_name):
				id_match = re.search(r"Name=([^;]+)", attributes)
			else:
				id_match = re.search(r"ID=([^;]+)", attributes)
			if id_match:
				id = id_match.group(1)
			else:
				# No Name field found, use gene field
				id_match = re.search(r"gene=([^;]+)", attributes)
				if id_match:
					id = id_match.group(1)
				else:
					# No gene field found either, use parent_name instead.
					# This is for IGH enhancers,
					# which comes with no Name nor gene field.
					id = parent_name

			# Write the transcript ID and gene_biotype to the output file
			gene_type_out.write(f"{id}\t{gene_biotype}\n")
			write_rgb_codes(id, gene_biotype)

			# Extract gene from transcript feature
			gene_name_match = re.search(r"gene=([^;]+)", attributes)
			if gene_name_match:
				gene_name = gene_name_match.group(1)
			
			if gene_name == parent_name and gene_synonyms != "":
				gene_names_out.write(f"{id}\t{gene_name}\t{gene_synonyms}\n")
			elif gene_name != parent_name and gene_synonyms != "":
				gene_names_out.write(f"{id}\t{gene_name}\t{parent_name},{gene_synonyms}\n")
			elif gene_synonyms == "":
				gene_names_out.write(f"{id}\t{gene_name}\t{parent_name}\n")

