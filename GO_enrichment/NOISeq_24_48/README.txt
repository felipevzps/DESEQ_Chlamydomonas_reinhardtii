# GO Annotations file
# GENE_ID<TAB>GO_ID1,GO_ID2,GO_ID3,...
cat DEGs_NOISeq_Creinhardtii_281_v5.6.annotation_info_24H.txt | tr " " "," | cut -f3,10 > GO_annotations.txt

# Gene List file
# Header (NOISeq list): CONTROL_mean    SALT_mean          M           D      prob      ranking
