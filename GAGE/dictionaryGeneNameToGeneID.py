translate_gene_name = {}

with open('GAGE/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.55.gff3', 'r') as gff_file:
    for line in gff_file:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'gene':
                attributes = fields[8].split(';')
                gene_name = gene_id = None
                for attr in attributes:
                    key, value = attr.strip().split('=')
                    if key == 'ID' and value.startswith('gene:'):
                        gene_id = value[5:] # Remove 'gene:' prefix
                    elif key == 'Name':
                        gene_name = value
                if gene_name and gene_id:
                    translate_gene_name[gene_id] = gene_name

#print(translate_gene_name.values())

# Read the file and rename gene names
with open('GAGE/cre_quantmerge.txt', 'r') as file:
    lines = file.readlines()

# Extract and preserve the header line
header = lines[0]

# Process and rename gene names for the remaining lines
for i in range(1, len(lines)):
    columns = lines[i].strip().split('\t')
    gene_name_with_suffix = columns[0]
    
    # Extract the gene name without the suffix (e.g., ".t1.2")
    gene_name = gene_name_with_suffix.split('.t')[0]

    # Check if the gene name is in the dictionary values (gene IDs)
    if gene_name in translate_gene_name.values():
        # Replace the gene name with the corresponding gene ID
        for key, value in translate_gene_name.items():
            if value == gene_name:
                columns[0] = key
                break  # Exit the loop once a match is found

    lines[i] = '\t'.join(columns) + '\n'  # Add back the newline character

# Write the modified lines (including the preserved header) back to the file
with open('GAGE/cre_quantmerge.translatedGeneName.txt', 'w') as file:
    file.write(header)
    file.writelines(lines[1:])











