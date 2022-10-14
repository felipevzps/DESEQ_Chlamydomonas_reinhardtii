# Plano de trabalho

## Encontro 01 (14/10):
* Baixar o dados direto de drive para o computador
 * Descompactar os dados: `unzip "*.zip"`
* Pegar os dados e transferir para o cluster: `scp -P 1000 . felipe.peres@bioinfo.cena.usp.br/`
* Avaliar a qualidade dops dados brutos com FasQC/MultiQC: `fastqc -f fastq ${INPUT} -o ${OUTPUT}`
* Fazer limpeza dos dados com o bbduk2
* Fazer de novo o FastQC/MultiQC para comparar os dados antes e depois do bbduk2
* Rodar Kraken2 pra avaliar a contaminação dos dados limpos de RNASeq
* Remover contaminação com ContFree-NGSrun_fastqc.sh

* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
 * Felipe - Montar apresentacao pra explicar os softwares que estamos usando

## Encontro 02 (19/10):
* Apresentacao sobre atividades do primeiro encontro (tudo que foi feito - desde download)

## Econtro 03 (26/10)
* Comparar dados limpos de RNASeq (antes e depois da remoção de contaminação)
* Quantificaçao ccom Salmon (TPM/CPM)
* DESEq Analysis (edgeR)

* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
  * Ivan - Procurar microalgas proximas da Chlamydomonas reinhardtii
  * Felipe - Estudar como usar edgeR ou Biocondutor pra avaliar expressçao diferencial

## Econtro 04 (02/11)
* Realizar notação de Gene Ontology (GO) com Pannzer
* Realizar analise de enriquecimento (termos GO enriquecidos nos transcritos mais expressos)

* Desenvolver o rascunho do nosso progresso (tópicos para o projeto)
  * Felipe, Ivan - escrever sobre nosso processo de quantificação e nossa ideia para DESeq Analysis

## Econtro 05 (09/10)

## Econtro 06 (16/10)
* Ivan vai estar viajando - vamos trabalhar online

## Econtro 07 (23/10) - dia antes da apresentacao


# TODO:
- Comparar GO enriquecidos da com outras microalgas 


# References
https://www.nature.com/articles/s41598-021-88954-6







