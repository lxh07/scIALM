# scIALM
A method for sparse scRNA-seq expression matrix imputation using the Inexact Augmented Lagrange Multiplier with low error.

### Dependencies
* Python: Preprocess the expression matrix & Evaluation of results
* Matlab:  Run Algorithm 3
* Dependencies can be installed using [inexact_alm/preprocess/requirements.txt](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/preprocess/requirements.txt)

### Data
* PBMC and Klein have been given.

### 1. Preprocess (Python)
* PBMC
  * Run [inexact_alm/preprocess/data/PBMC/preprocess.py](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/preprocess/data/PBMC/preprocess.py)
  * Output folder: inexact_alm/preprocess/data/PBMC/processed/
* Klein
  * Run [inexact_alm/preprocess/data/Klein/preprocess.py](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/preprocess/data/Klein/preprocess.py)
  * Output folder: inexact_alm/preprocess/data/Klein/processed/
 
### 2. Mask (Python)
* Run [inexact_alm/preprocess/data/mask.py](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/preprocess/data/mask.py)
* Output folder: inexact_alm/preprocess/data/{name}/masked/

### 3. Algorithm 3 (Matlab)
* Run [inexact_alm/test.m](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/test.m)
* Output folder: inexact_alm/preprocess/data/{name}/result/

### 4. Results (Python)
* Run [inexact_alm/preprocess/result.py](https://github.com/lxhfighting/scIALM/blob/main/inexact_alm/preprocess/result.py)

### Set dataset and masking ratio
If you want to use a different dataset and masking ratio, please modify the following code:

mask.py & result.py
```
parser = argparse.ArgumentParser()
parser.add_argument('--masked_prob', default=0.1, type=float)
parser.add_argument('--dataset', default='PBMC', type=str)
```
test.m
```matlab
mask = csvread("preprocess/data/PBMC/masked/PBMC_01.csv",1,1);

writetable(data_table, "preprocess/data/PBMC/result/PBMC_01.csv");
```
