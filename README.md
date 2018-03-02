# UCLA iGEM 2018
UCLA International Genetically Engineered Machines Team 2018 & Kosuri Lab, UCLA

Promoters are the key drivers of gene expression. We are interested in manipulating operator sequences in the lac
inducible promoter. 

## Deployment

To generate the promoter library:
```
python lac_lib_short.py
```

To filter out any bad promoters:
```
python cleanup_REs.py lac_lib_short.txt RE.fasta output_file output_type
```

To generate the positive controls:
```
python generate_controls.py synthetic_promoter.csv fwd_primers.fasta rev_primers.fasta
```

## Contributors
* Timothy Yu
* Members of UCLA iGEM and the Kosuri Lab
