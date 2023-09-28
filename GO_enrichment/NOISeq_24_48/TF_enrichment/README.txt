# Gene List
for i in `cat DEGs_NOISeq_CRE_v5.6.hmmsearch_PlnTFDVPfam32.families.24H.txt | cut -f1`; do grep $i ../DEGs_NOISeqTech_SALT24H.mod.mod.csv >> DEGs_NOISeq_CRE_v5.6.hmmsearch_PlnTFDVPfam32.families.24H_scores.txt; done
