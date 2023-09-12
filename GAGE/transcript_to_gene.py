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
                        gene_id = value[5:] # Remove 'gene:' prefix, single quotes, and extra whitespace
                    elif key == 'Name':
                        gene_name = value # Remove single quotes and extra whitespace
                if gene_name and gene_id:
                    translate_gene_name[gene_id] = gene_name

print(translate_gene_name)