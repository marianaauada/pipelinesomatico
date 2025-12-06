# pipelinesomatico
Pipeline Somático - Do VCF (anotado) até o CGI Classificação

**1. Clonar o git lmabrasil-hg38**

```bash
! git clone https://github.com/renatopuga/lmabrasil-hg38.git
```

Output
```
Cloning into 'lmabrasil-hg38'...
remote: Enumerating objects: 226, done.
remote: Counting objects: 100% (168/168), done.
remote: Compressing objects: 100% (108/108), done.
remote: Total 226 (delta 90), reused 114 (delta 56), pack-reused 58 (from 1)
Receiving objects: 100% (226/226), 8.63 MiB | 20.05 MiB/s, done.
Resolving deltas: 100% (106/106), done.
```

**2. Agora, vá até o github lmabrasil-hg38 na sessão usando o CGI via API Reset no Google Colab**

```bash
%%bash
# Cortar pelas colunas de 1 a 4 e criar um novo arquivo chamado df_WP048-cgi.txt
# 1: CHROM (COnverte CHROM para CHR) formato que o CGI gosta
# 2: POS
# 3: REF
# 4: ALT
cut -f1-4 /content/lmabrasil-hg38/vep_output/liftOver_WP048_hg19ToHg38.vep.filter.tsv | sed -e "s/CHROM/CHR/g"  > df_WP048-cgi.txt

# Listar as 10 primeiras linhas
head df_WP048-cgi.txt
```

Output
```
CHR	POS	REF	ALT
chr1	114716123	C	T
chr9	5073770	G	T
```

**3. Enivar job para Cancer Genome Interpret (CGI) API**
> https://www.cancergenomeinterpreter.org%2Frest_api

Após filtrar apenas as colunas de interesse (CHR, POS, REF e ALT), agora podemos enviar via EST_API as variantes somáticas da amostra WP048

> **Nota:** Altere a variável {SEU_TOKEN} para o Token do CGI criado para sua conta
```python
import requests
headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN}'}
payload = {'cancer_type': 'HEMATO', 'title': 'Somatic MF WP048', 'reference': 'hg38'}
r = requests.post('https://www.cancergenomeinterpreter.org/api/v1',
                headers=headers,
                files={
                        'mutations': open('/content/df_WP048-cgi.txt', 'rb')
                        },
                data=payload)
r.json()
```

Output Job ID:
```
7d09666f743c78387299
```

**4. Status do Job ID (Running, Error ou Done)**

Verifique o ´Status` do seu job id para identificar se a analise terminou ou houve algum erro.

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN}'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers)
r.json()
```

Output:
```
{'status': 'Done',
 'metadata': {'id': '7d09666f743c78387299',
  'user': 'marianabelloauada@gmail.com',
  'title': 'Somatic MF WP048',
  'cancertype': 'HEMATO',
  'reference': 'hg38',
  'dataset': 'input.tsv',
  'date': '2025-12-06 14:05:22'}}
  ```

**5. Log completo do Job**

Aqui podemos verificar o status em cada uma das etapas da análise do CGI

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN}'}
payload={'action':'logs'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
r.json()
```

Output:
```
{'status': 'Done',
 'logs': ['# cgi analyze input.tsv -c HEMATO -g hg38',
  '2025-12-06 15:05:25,296 INFO     Parsing input01.tsv\n',
  '2025-12-06 15:05:29,309 INFO     Running VEP\n',
  '2025-12-06 15:05:30,260 INFO     Check cancer genes and consensus roles\n',
  '2025-12-06 15:05:30,344 INFO     Annotate BoostDM mutations\n',
  '2025-12-06 15:05:30,379 INFO     Annotate OncodriveMUT mutations\n',
  '2025-12-06 15:05:32,672 INFO     Annotate validated oncogenic mutations\n',
  '2025-12-06 15:05:32,828 INFO     Check oncogenic classification\n',
  '2025-12-06 15:05:32,894 INFO     Matching biomarkers\n',
  '2025-12-06 15:05:32,987 INFO     Prescription finished\n',
  '2025-12-06 15:05:32,999 INFO     Aggregate metrics\n',
  '2025-12-06 15:05:35,659 INFO     Compress output files\n',
  '2025-12-06 15:05:35,685 INFO     Analysis done\n']}
```


**6. Download dos Resultados**

Total de 4 arquivos de resultados:
> A definição de cada arquivo pelo CGI (ver no site)

1. alterations.tsv:
2. biomarker.tsv:
3. input01.tsv:
4. summary.txt:

**7. Criar o diretório para resultados para cada amostra**

```bash
%%bash
Criar o diretório com 0 ID da amostra dentro de results
mkdir -p results/WP048
```

**8. Fazer download do arquivo `.zip`**

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN}'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
with open('/content/results/WP048/WP048-cgi.zip', 'wb') as fd:
    fd.write(r._content)
``` 

**9. Descompactar o arquivo '.zip' no diretório de resultados da amostra**

```bash
%%bash
unzip /content/results/WP048/WP048-cgi.zip -d /content/results/WP048/
```

**10. Agora podemos visualizar a tabela `alterations.tsv`e descobrir quais alterações somáticas são `Drivers`, `Passengers` ou `Unclass**

**Visualizar a tabela `alterations.tsv`**

Instalar a lib `pandas`

```bash
! pip install pandas
```

```python
import pandas as pd
pd.read_csv('/content/results/WP048/alterations.tsv',sep='\t',index_col=False, engine= 'python')
```

Output:
|Index|Input ID|CHROMOSOME|POSITION|REF|ALT|CHR|POS|ALT\_TYPE|STRAND|CGI-Sample ID|CGI-Gene|CGI-Protein Change|CGI-Oncogenic Summary|CGI-Oncogenic Prediction|CGI-External oncogenic annotation|CGI-Mutation|CGI-Consequence|CGI-Transcript|CGI-STRAND|CGI-Type|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|input01\_1|1|114716123|C|T|chr1|114716123|snp|+|input01|NRAS|G13D|oncogenic \(predicted and annotated\)|driver \(boostDM: non-tissue-specific model\)|cgi,oncokb,clinvar:13901|chr1:114716123 C\>T|missense\_variant|ENST00000369535|+|SNV|
|1|input01\_2|9|5073770|G|T|chr9|5073770|snp|+|input01|JAK2|V617F|oncogenic \(annotated\)|passenger \(oncodriveMUT\)|cgi,oncokb,clinvar:14662|chr9:5073770 G\>T|missense\_variant|ENST00000381652|+|SNV|


# WP017

|index|Input ID|CHROMOSOME|POSITION|REF|ALT|CHR|POS|ALT\_TYPE|STRAND|CGI-Sample ID|CGI-Gene|CGI-Protein Change|CGI-Oncogenic Summary|CGI-Oncogenic Prediction|CGI-External oncogenic annotation|CGI-Mutation|CGI-Consequence|CGI-Transcript|CGI-STRAND|CGI-Type|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|input01\_1|1|114716127|C|T|chr1|114716127|snp|+|input01|NRAS|G12S|oncogenic \(predicted and annotated\)|driver \(boostDM: non-tissue-specific model\)|cgi,clinvar:177778|chr1:114716127 C\>T|missense\_variant|ENST00000369535|+|SNV|
|1|input01\_2|1|152304661|G|C|chr1|152304661|snp|+|input01|FLG|R3409G|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr1:152304661 G\>C|missense\_variant|ENST00000368799|+|SNV|
|2|input01\_3|11|115209621|G|A|chr11|115209621|snp|+|input01|CADM1|T344I|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr11:115209621 G\>A|missense\_variant|ENST00000331581|+|SNV|
|3|input01\_4|12|132140028|C|T|chr12|132140028|snp|+|input01|DDX51|--|non-protein affecting|non-protein affecting|NaN|chr12:132140028 C\>T|intron\_variant|ENST00000397333|+|SNV|
|4|input01\_5|14|24119817|G|A|chr14|24119817|snp|+|input01|DCAF11|R338H|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr14:24119817 G\>A|missense\_variant|ENST00000446197|+|SNV|
|5|input01\_6|14|59727407|G|A|chr14|59727407|snp|+|input01|RTN1|A426V|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr14:59727407 G\>A|missense\_variant|ENST00000267484|+|SNV|
|6|input01\_7|15|24675943|G|A|chr15|24675943|snp|+|input01|NPAP1|A26T|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr15:24675943 G\>A|missense\_variant|ENST00000329468|+|SNV|
|7|input01\_8|16|67616834|G|A|chr16|67616834|snp|+|input01|CTCF|E348K|oncogenic \(predicted\)|driver \(oncodriveMUT\)|NaN|chr16:67616834 G\>A|missense\_variant|ENST00000646076|+|SNV|
|8|input01\_9|16|71389851|G|A|chr16|71389851|snp|+|input01|CALB2|E268K|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr16:71389851 G\>A|missense\_variant|ENST00000302628|+|SNV|
|9|input01\_10|16|74452195|A|C|chr16|74452195|snp|+|input01|GLG1|--|non-protein affecting|non-protein affecting|NaN|chr16:74452195 A\>C|intron\_variant|ENST00000205061|+|SNV|
|10|input01\_11|19|12943751|GCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAG|-|chr19|12943750|indel|+|input01|CALR|EQRLKEEEEDKKRKEEEE364-381X|oncogenic \(predicted\)|driver \(oncodriveMUT\)|NaN|chr19:12943751-12943751 GCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAG\>-|frameshift\_variant|ENST00000316448|+|DEL|
|11|input01\_12|19|35673516|A|C|chr19|35673516|snp|+|input01|UPK1A|T147P|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr19:35673516 A\>C|missense\_variant|ENST00000616789|+|SNV|
|12|input01\_13|19|45319445|T|G|chr19|45319445|snp|+|input01|CKM|--|non-protein affecting|non-protein affecting|NaN|chr19:45319445 T\>G|intron\_variant|ENST00000221476|+|SNV|
|13|input01\_14|2|137450992|G|A|chr2|137450992|snp|+|input01|THSD7B|R1036Q|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr2:137450992 G\>A|missense\_variant|ENST00000272643|+|SNV|
|14|input01\_15|2|218880905|A|C|chr2|218880905|snp|+|input01|WNT10A|--|non-protein affecting|non-protein affecting|NaN|chr2:218880905 A\>C|5\_prime\_UTR\_variant|ENST00000258411|+|SNV|
|15|input01\_16|20|32434638|-|G|chr20|32434638|indel|+|input01|ASXL1|-642-643X|oncogenic \(predicted\)|driver \(oncodriveMUT\)|NaN|chr20:32434638-32434639 -\>G|frameshift\_variant|ENST00000375687|+|INS|
|16|input01\_17|21|41901750|C|T|chr21|41901750|snp|+|input01|C2CD2|--|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr21:41901750 C\>T|splice\_acceptor\_variant|ENST00000380486|+|SNV|
|17|input01\_18|21|43094667|T|G|chr21|43094667|snp|+|input01|U2AF1|Q157P|oncogenic \(predicted and annotated\)|driver \(boostDM: non-tissue-specific model\)|cgi,oncokb,clinvar:376024|chr21:43094667 T\>G|missense\_variant|ENST00000291552|+|SNV|
|18|input01\_20|3|133380804|G|A|chr3|133380804|snp|+|input01|TMEM108|D365N|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr3:133380804 G\>A|missense\_variant|ENST00000321871|+|SNV|
|19|input01\_22|7|584561|G|A|chr7|584561|snp|+|input01|PRKAR1B|T239M|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr7:584561 G\>A|missense\_variant|ENST00000406797|+|SNV|


# WP019

|Index|Input ID|CHROMOSOME|POSITION|REF|ALT|CHR|POS|ALT\_TYPE|STRAND|CGI-Sample ID|CGI-Gene|CGI-Protein Change|CGI-Oncogenic Summary|CGI-Oncogenic Prediction|CGI-External oncogenic annotation|CGI-Mutation|CGI-Consequence|CGI-Transcript|CGI-STRAND|CGI-Type|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|input01\_1|1|16031913|G|A|chr1|16031913|snp|+|input01|CLCNKA|--|non-protein affecting|non-protein affecting|NaN|chr1:16031913 G\>A|intron\_variant|ENST00000331433|+|SNV|
|1|input01\_2|1|149073703|C|T|chr1|149073703|snp|+|input01|NBPF9|--|non-protein affecting|non-protein affecting|NaN|chr1:149073703 C\>T|intron\_variant|ENST00000615421|+|SNV|
|2|input01\_3|17|76736877|G|T|chr17|76736877|snp|+|input01|SRSF2|P95H|oncogenic \(predicted and annotated\)|driver \(oncodriveMUT\)|cgi,oncokb|chr17:76736877 G\>T|missense\_variant|ENST00000392485|+|SNV|
|3|input01\_6|2|113117932|G|T|chr2|113117932|snp|+|input01|IL1RN|--|non-protein affecting|non-protein affecting|NaN|chr2:113117932 G\>T|5\_prime\_UTR\_variant|ENST00000259206|+|SNV|
|4|input01\_7|4|2341569|G|C|chr4|2341569|snp|+|input01|ZFYVE28|P76R|non-oncogenic|passenger \(oncodriveMUT\)|NaN|chr4:2341569 G\>C|missense\_variant|ENST00000290974|+|SNV|
|5|input01\_11|6|158130554|A|C|chr6|158130554|snp|+|input01|SERAC1|--|non-protein affecting|non-protein affecting|NaN|chr6:158130554 A\>C|intron\_variant|ENST00000647468|+|SNV|
|6|input01\_13|7|111783014|G|A|chr7|111783014|snp|+|input01|DOCK4|--|non-protein affecting|non-protein affecting|NaN|chr7:111783014 G\>A|intron\_variant|ENST00000445943|+|SNV|
|7|input01\_14|9|5073770|G|T|chr9|5073770|snp|+|input01|JAK2|V617F|oncogenic \(annotated\)|passenger \(oncodriveMUT\)|cgi,oncokb,clinvar:14662|chr9:5073770 G\>T|missense\_variant|ENST00000381652|+|SNV|


# WP058

|Index|Input ID|CHROMOSOME|POSITION|REF|ALT|CHR|POS|ALT\_TYPE|STRAND|CGI-Sample ID|CGI-Gene|CGI-Protein Change|CGI-Oncogenic Summary|CGI-Oncogenic Prediction|CGI-External oncogenic annotation|CGI-Mutation|CGI-Consequence|CGI-Transcript|CGI-STRAND|CGI-Type|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|input01\_1|12|57185563|T|G|chr12|57185563|snp|+|input01|LRP1|C2166G|NaN|NaN|NaN|chr12:57185563 T\>G|missense\_variant|ENST00000243077|+|SNV|
|1|input01\_2|15|28272094|A|G|chr15|28272094|snp|+|input01|HERC2|--|non-protein affecting|non-protein affecting|NaN|chr15:28272094 A\>G|intron\_variant|ENST00000261609|+|SNV|
|2|input01\_3|15|43206304|T|G|chr15|43206304|snp|+|input01|AC068724\.4|--|non-protein affecting|non-protein affecting|NaN|chr15:43206304 T\>G|intron\_variant|ENST00000563128|+|SNV|
|3|input01\_4|19|12943751|GCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAG|-|chr19|12943750|indel|+|input01|CALR|EQRLKEEEEDKKRKEEEE364-381X|NaN|NaN|NaN|chr19:12943751-12943751 GCAGAGGCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAG\>-|frameshift\_variant|ENST00000316448|+|DEL|
|4|input01\_5|2|105892889|T|G|chr2|105892889|snp|+|input01|NCK2|--|non-protein affecting|non-protein affecting|NaN|chr2:105892889 T\>G|intron\_variant|ENST00000233154|+|SNV|
|5|input01\_6|20|32434638|-|G|chr20|32434638|indel|+|input01|ASXL1|-642-643X|NaN|NaN|NaN|chr20:32434638-32434639 -\>G|frameshift\_variant|ENST00000375687|+|INS|
|6|input01\_7|7|124892275|T|C|chr7|124892275|snp|+|input01|POT1|K39E|NaN|NaN|NaN|chr7:124892275 T\>C|missense\_variant|ENST00000357628|+|SNV|


# WP068

|Index|Input ID|CHROMOSOME|POSITION|REF|ALT|CHR|POS|ALT\_TYPE|STRAND|CGI-Sample ID|CGI-Gene|CGI-Protein Change|CGI-Oncogenic Summary|CGI-Oncogenic Prediction|CGI-External oncogenic annotation|CGI-Mutation|CGI-Consequence|CGI-Transcript|CGI-STRAND|CGI-Type|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|0|input01\_1|1|114716123|C|T|chr1|114716123|snp|+|input01|NRAS|G13D|oncogenic \(predicted and annotated\)|driver \(boostDM: non-tissue-specific model\)|cgi,oncokb,clinvar:13901|chr1:114716123 C\>T|missense\_variant|ENST00000369535|+|SNV|
|1|input01\_2|9|5073770|G|T|chr9|5073770|snp|+|input01|JAK2|V617F|oncogenic \(annotated\)|passenger \(oncodriveMUT\)|cgi,oncokb,clinvar:14662|chr9:5073770 G\>T|missense\_variant|ENST00000381652|+|SNV|
