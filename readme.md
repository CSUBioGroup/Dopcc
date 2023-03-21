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

1. Run get_embeddings.py to generate embeddings (Hierarchical compression + Node embedding + PCA)
```
python generate_embedding.py --input ./datasets/BioGRID_network.txt --output ./emds/BioGRID_network.txt --representation-size 256 --method deepwalk --hs_num 3
```
2. Run Core Attachment method
```
python dopcc.py --input ./datasets/BioGRID_network.txt --output ./results/BioGRID_network.txt --emd ./emds/BioGRID_network.txt
```
