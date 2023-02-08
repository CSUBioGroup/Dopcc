# Dopcc

Detecting overlapping protein complexes via multi-metrics and co-core attachment method

## Datesets

- BioGRID: datasets/BioGRID_network.txt
- IntAct: datasets/IntAct_network.txt

## Reference

- CYC2008: refs/CYC2008_3_236.txt
- Combiner: refs/yeast_combine.txt
- sub reference for overlapping complexes: refs/yeast_overlaps.txt

## embeddings

- BioGRID: emds/BioGRID_network_1024_redu_256.txt
- IntAct: emds/IntAct_network_1024_redu_256.txt

## Usage

run demo:
```
python dopcc.py
```

use specific datasets

```
python dopcc.py --input dataset_path --output result_path --emd embedding_path
```
