# Plano de trabalho

## Encontro 01 (14/10):
- [x] Baixar RNA-SEQ da Chlamydomonas reinhardtii (está armazenado no Google Drive)
  - [x] Descompactar pacotes de dados gerados pelo Google Drive 
  ```bash
  unzip "*.zip"
  ```
  - [x] Transferir os dados brutos de RNA-Seq para o cluster do CENA
  ```bash
  scp -P 1000 *. felipe.peres@bioinfo.cena.usp.br/Storage/data1/felipe.peres/DESEQ_Chlamydomonas_reinhardtii
  ```
- [x] Avaliar a qualidade dos dados brutos com FasQC/MultiQC
  ```bash
  #!/bin/bash
  #$ -cwd
  #$ -V
  #$ -q all.q
  #$ -t 1-32
  #$ -tc 4
  #$ -pe smp 4
  #$ -l h=figsrv
  
  INFILES=`ls -1 *.fastq.gz | head -n $SGE_TASK_ID | tail -n1`
  OUT_FASTQC=${IN/.fastq.gz/.fastqc_report}
  
  module load FastQC/0.11.8
  fastqc -f fastq ${INFILES} -o ${OUT_FASTQC}
  ```
  * OBS: O conteudo GC bate com a literatura (64%). REF: https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=2&id=100524
- [x] Fazer limpeza dos dados com o bbduk2
  ```bash
  #!/bin/bash
  #$ -cwd
  #$ -V
  #$ -q all.q
  #$ -t 1-16
  #$ -tc 4
  #$ -pe smp 4
  #$ -l h=figsrv
  
  OUT=/Storage/data1/felipe.peres/DESEQ_Chlamydomonas_reinhardtii/trimmed_reads
  infileR1=`ls -1 raw_reads/*R1*.fastq.gz |head -n $SGE_TASK_ID|tail -n1`
  infileR2=${infileR1/R1_001.fastq.gz/R2_001.fastq.gz}
  outbase=`basename $infileR1`
  outbase=${outbase/_R1_001.fastq.gz}
  outfileR1=${outbase}_001_R1_trimmed.fastq.gz
  outfileR2=${outbase}_001_R2_trimmed.fastq.gz
  outfile_refstats=${outbase}.trimmed.refstats
  outfile_stats=${outbase}.trimmed.stats
  
  module load bbmap/35.85
  
  bbduk2.sh -Xmx40g threads=$NSLOTS in1=$infileR1 in2=$infileR2 refstats=$OUT/$outfile_refstats stats=$OUT/$outfile_stats
  out1=$OUT/$outfileR1 out2=$OUT/$outfileR2 minlength=75 rref=/Storage/progs/bbmap_35.85/resources/adapters.fa
  fref=/Storage/progs/sortmerna-2.1b/rRNA_databases/rfam-5.8s-database-id98.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/rfam-5s-database-id98.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-bac-23s-id98.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-arc-16s-id95.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-euk-18s-id95.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-arc-23s-id98.fasta,/Storage/progs/sortmerna-2.1b/rRNA_databases/silva-euk-28s-id98.fasta qtrim=w trimq=20 tpe tbo
  ```
- [x] Avaliar a qualidade dos dados limpos com FasQC/MultiQC
  ```bash
  #!/bin/bash
  #$ -cwd
  #$ -V
  #$ -q all.q
  #$ -t 1-32
  #$ -tc 4
  #$ -pe smp 4
  #$ -l h=figsrv
  
  INFILES=`ls -1 trimmed_reads/*.fastq.gz | head -n $SGE_TASK_ID | tail -n1`
  OUT_FASTQC=${IN/.fastq.gz/.fastqc_report}
  
  module load FastQC/0.11.8
  fastqc -f fastq ${INFILES} -o ${OUT_FASTQC}
  ```
- [x] Rodar Kraken2 pra avaliar a contaminação dos dados limpos de RNASeq
  ```bash
  #!/bin/bash
  #$ -q all.q
  #$ -V
  #$ -cwd
  #$ -t 1-16
  #$ -tc 1
  #$ -pe smp 10
  
  DBNAME=krakendb
  R1=`ls -1 trimmed_reads/*R1_trimmed.fastq.gz | head -n $SGE_TASK_ID | tail -n1`
  R2=${R1/R1_trimmed.fastq.gz/R2_trimmed.fastq.gz}
  report_name=`basename $R1`
  kraken_quant=${report_name/R1_trimmed.fastq.gz}kraken_quant
  kraken_output=${report_name/R1_trimmed.fastq.gz}kraken_output
  
  module load Kraken2/2.1.2
  
  /usr/bin/time -v kraken2 --db ./../Kraken2/$DBNAME --report-zero-counts --report $kraken_quant --output $kraken_output --paired $R1 $R2 --threads $NSLOTS
  ```
   * OBS: Dados classificados como "S28" e "S9" tem respectivamente 13% e 26% de reads não classificadas ... O que são?
      * Se forem transcritos importantes, provavelmente vão aparecer na comparação da DESeq Analysis antes/depois da contaminação.

# TODO:
* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
 * Felipe - Montar apresentacao pra explicar os softwares que estamos usando
 * Remover contaminação com ContFree-NGS

## Encontro 02 (19/10)

## Econtro 03 (26/10)
* Apresentacao sobre atividades do primeiro encontro (tudo que foi feito - desde download)
* Comparar dados limpos de RNASeq (antes e depois da remoção de contaminação)
* Quantificaçao ccom Salmon (TPM/CPM)

## Econtro 04 (02/11)
* Felipe - Estudar como usar edgeR ou Biocondutor pra avaliar expressçao diferencial
  * https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html 
* DESEq Analysis (edgeR)

* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
  * Ivan - Procurar microalgas proximas da Chlamydomonas reinhardtii https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=3042 

* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
  * Felipe, Ivan - escrever sobre nosso processo de quantificação e nossa ideia para DESeq Analysis

## Econtro 05 (09/11)
    - O que temos:
        - [x]  Dados crus
        - [x]  Dados limpos
        - [x]  Tabela do Salmon **- Felipe vai fazer**
            
            Sample Name, Sample description, Total Reads, Total Mapped, Ratio of Mapped Reads (%)
            
            ```bash
            #Pegando reads mapeadas
            cat LOM*/logs/salmon_quant.log | grep "%"
            
            #pwd
            /home/felipevzps/Documentos/DESEQ_Chlamydomonas_reinhardtii/DESeq2****************
            ```
            
        - [x]  Foto da Chlamydomonas reinhardtii **- Ivan vai fazer**
        - [x]  Graficos
            - [x]  PCA da distribuicao das amostras
                - [x]  Time
                - [x]  Condition
                - [x]  Time and Condition
                - [ ]  Grafico da distribuicao das medias de counts **- Felipe e Ivan**
                - [ ]  Heatmaps expressao diferencial dos genes da Chlamydomonas reinhardtii -  **Felipe e Ivan**
        - [ ]  Transcriptoma *de novo* da Chlamydomonas reinhardtii
            - [ ]  Comparar nosso transcriptoma com transcriptoma de referencia da Chlamydomonas reinhardtii (Diagrama de Venn - genes BUSCO?) — **Talvez pensar em outra estrategia conversando com Diego**
            - [ ]  Genes BUSCO do transcriptoma - Ivan
            - [ ]  Metricas TRANSRATE do transcriptoma - Ivan
                - [ ]  Rodar Panzer pra pegar os termos GO - **ANTES DO DIA 13**
                    - [ ]  Analisar termos GO enriquecidos - **Exact Fisher’s Test - Perguntar pro Diego qual classe de GO devemos usar (BP, MF, CC) - DIA 17 talvez**
                - [ ]  Por fim, iremos avaliar a conservação dos termos GO enriquecidos nos genes diferencialmente expressos de Chlamydomonas reinhardtii com outras microalgas (a ideia é observar a conservação destes genes).
                - [ ]  **TERMINAR DE ESCREVER ATE QUARTA 16/10**


## Econtro 06 (16/11)
* Ivan vai estar viajando - vamos trabalhar online
* Realizar analise de enriquecimento (termos GO enriquecidos nos transcritos mais expressos)


## Econtro 07 (23/11) - dia antes da apresentacao


# TODO:
- Comparar GO enriquecidos da com outras microalgas 


# References
https://www.nature.com/articles/s41598-021-88954-6







