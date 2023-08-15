transcript_to_gene = {}

with open('Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.55.gff3', 'r') as gff_file:
    for line in gff_file:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'mRNA':
                attributes = fields[8].split(';')
                transcript_id = gene_id = None
                for attr in attributes:
                    key, value = attr.strip().split('=')
                    if key == 'ID':
                        transcript_id = value
                    elif key == 'Parent':
                        gene_id = value
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id

print(transcript_to_gene)
