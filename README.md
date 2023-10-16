[DESEQ: Chlamydomonas reinhardtii wiki](https://github.com/felipevzps/DESEQ_Chlamydomonas_reinhardtii/wiki)

- **Analisar TFs possivelmente relacionados ao acúmulo de liídeos**
    - [ ]  Criar uma lista com TFs que “somem” de 24h para 48h
        - **Eles podem estar relacionados ao acúmulo de lipídeos (Envolvidos no switch desse sistema)**
    - [ ]  Procurar funções destes TFs
        - **Eles fazem as células pararem de produzir lipídeos?**
- **Identificar motivos conservados nos promotores dos DEGs (e também nos promotores dos TFs) - MEME**
    - [ ]  Extrair sequência dos promotores de todos os genes
    - [ ]  Analisar quais genes os motivos conservados dos DEGs aparecem
        - Resultando em possíveis genes alvos reguladores
        - Pipeline no artigo de Fábio e Larissa
        - Diego: Seria basicamente buscar por palavras altamente conservadas nos promotores e comparar essas palavras com promotores em todo o genoma.
        - Focar nas condições 24h vs 48h
    - [ ]  Anotar os termos GO dos genes com promotores conservados
    - [ ]  Comparar motivos conservados nos promotores promotores de DEGs **UP** e **DOWN** regulated
        - Comparar os dois grupos além do background, **UP** vs **DOWN**, **DOWN** vs **UP** (resultando em 3 comparações )
- **Comparação de termos GO (estresse salino e deficiência de nitrogênio)**
    - Aparentemente o NaCL ativa o mesmo switch metabolico da deficiência de nitrogênio.
    - Procurar se os mesmos termos GO estão envolvidos nas duas condições.
    - [ ]  Comparar termos GO do paper de nitrogênio e os nossos termos GO. (Comparar em uma tabela)
        - Artigo do Luca com termos do nitrogenio
- **Visualização das rotas metabólicas**
    - O KEGG mostra que a fotossíntese está ativa, mas o que está acontecendo nesta via?
    - [ ]  Usar pathway tools pra visualizar quais genes ou enzimas estão ativas na via dos termos mais significativos (**Score < 0.05**)





