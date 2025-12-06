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

**Agora, vá até o github lmabrasil-hg38 na sessão usando o CGI via API Reset no Google Colab**

```bash
%%bash
# Cortar pelas colunas de 1 a 4 e criar um novo arquivo chamado df_WP848-cgi.txt
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

**Enivar job para Cancer Genome Interpret (CGI) API**
> https://www.google.com/url?q=https%3A%2F%2Fwww.cancergenomeinterpreter.org%2Frest_api

Após filtrar apenas as colunas de interesse (CHR, POS, REF e ALT), agora podemos enviar via EST_API as variantes somáticas da amostra WP048

> **Nota:** Altere a variável {SEU_TOKEN} para o Token do CGI criado para sua conta
```python
import requests
headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN'}
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

**Status do Job ID (Running, Error ou Done)**

Verifique o Status` do seu job id para identificar se a analise terminou ou houve algum erro.

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN'}
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

**Log completo do Job**

Aqui podemos verificar o status em cada uma das etapas da análise do CGI

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com be5873853bda53991f05'}
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


**Download dos Resultados**

Total de 4 arquivos de resultados:
> A definição de cada arquivo pelo CGI (ver no site)

1. alterations.tsv:
2. biomarker.tsv:
3. input01.tsv:
4. summary.txt:

**Criar o diretório para resultados para cada amostra**

```bash
%%bash
Criar o diretório com 0 ID da amostra dentro de results
mkdir -p results/WP048
```

**Fazer download do arquivo `.zip`**

```python
import requests
job_id ="7d09666f743c78387299"

headers = {'Authorization': 'marianabelloauada@gmail.com {SEU_TOKEN}'}
payload={'action':'download'}
r = requests.get('https://www.cancergenomeinterpreter.org/api/v1/%s' % job_id, headers=headers, params=payload)
with open('/content/results/WP048/WP048-cgi.zip', 'wb') as fd:
    fd.write(r._content)
``` 

**Descompactar o arquivo '.zip' no diretório de resultados da amostra**

```bash
%%bash
unzip /content/results/WP048/WP048-cgi.zip -d /content/results/WP048/
```

**Agora podemos visualizar a tabela `alterations.tsv`e descobrir quais alterações somáticas são `Drivers`, `Passengers` ou `Unclass**

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
