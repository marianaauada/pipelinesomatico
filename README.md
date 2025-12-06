# pipelinesomatico
Aula Pipeline Somático - Do VCF (anotado) até o CGI Classificação

**1. Clonar o git lmabrasil-hg38**

```bash
! git clone https://github.com/renatopuga/lmabrasil-hg38.git
```

output
```
Cloning into 'lmabrasil-hg38'...
remote: Enumerating objects: 226, done.
remote: Counting objects: 100% (168/168), done.
remote: Compressing objects: 100% (108/108), done.
remote: Total 226 (delta 90), reused 114 (delta 56), pack-reused 58 (from 1)
Receiving objects: 100% (226/226), 8.63 MiB | 20.05 MiB/s, done.
Resolving deltas: 100% (106/106), done.
```

Agora, vá até o gitbut lmabrasil-hg38 na sessão usando o CGI via API Reset no Google Colab

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

output
```
CHR	POS	REF	ALT
chr1	114716123	C	T
chr9	5073770	G	T
```

**Enivar job para Cancer Genome Interpret (CGI) API**
> https://www.google.com/url?q=https%3A%2F%2Fwww.cancergenomeinterpreter.org%2Frest_api

Após filtrar apenas as colunas de interesse (CHR, POS, REF e ALT), agora podemos enviar via EST_API as variantes somáticas da amostra WP048

> Nota: Altere a variável {SEU_TOKEN} para o Token do CGI criado para sua conta
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

output Job ID:
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

output:
´´´
{'status': 'Done',
 'metadata': {'id': '7d09666f743c78387299',
  'user': 'marianabelloauada@gmail.com',
  'title': 'Somatic MF WP048',
  'cancertype': 'HEMATO',
  'reference': 'hg38',
  'dataset': 'input.tsv',
  'date': '2025-12-06 14:05:22'}}
  ´´´
`
**Log completo do Job**


Aqui podemos verificar o status em cada uma das etapas da análise do CGI

** Download dos Resultados**

Total de 4 arquivos de rssultados:

