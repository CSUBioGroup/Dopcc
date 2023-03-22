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

### Run demo:
```
python Dopcc.py
```

### Apply Dopcc on a specific dataset

1. Run get_embeddings.py to generate embeddings (Hierarchical compression + Node embedding + PCA)
```
python generate_embedding.py --input ./datasets/BioGRID_network.txt --output ./emds/BioGRID_network.txt --representation-size 256 --method deepwalk --hs_num 3
```
2. Run Core Attachment method
```
python Dopcc.py --input ./datasets/BioGRID_network.txt --output ./results/BioGRID_network.txt --emd ./emds/BioGRID_network.txt
```

### Parameters
#### generate_embeddings.py

> --input: PPI network
>
> --output: embedding results
>
> --represesntation-size: final embedding size (= each layer embedding size)
>
> --hs_num: the counts of Hierarchical compression

#### Dopcc.py
> --input: PPI network
>
> --output: predicted protein complexes
>
> --emd: embedding results from the generate_embedding.py


## Concat
Please feel free to contact us for any further questions.

- Wenkang Wang wangwk@csu.edu.cn
- Min Li limin@mail.csu.edu.cn
