# C1-C6 reviewer-defense additions (extended panel)

**Date**: 2026-04-20

All in-silico, no wet-lab.


## C1 antiSMASH BGC scan

`results/additional/C1_antismash/antismash_per_MAG.csv` + `_hero_vs_rest.csv`

```
(missing)
```

## C2 DefenseFinder + minced CRISPR

```
              metric  hero_mean  rest_mean  MWU_p
0  n_defense_systems        0.0        0.0    1.0
1    n_crispr_arrays        0.0        0.0    1.0
```

## C3 dN/dS + codon usage

codeml M0 ω per urease gene:

```


```

Codon usage hero vs rest:

```
    metric  hero_mean  rest_mean     MWU_p
0  GC3_pct  51.413327  70.065052  0.043288
1      ENC  49.988579  39.004205  0.002830
```

## C4 ESMFold UreC vs PDB 4CEU

```
(missing)
```

## C5 Pan-MICP environment ANI

```
Ref_file	Query_file	ANI	Align_fraction_ref	Align_fraction_query	Ref_name	Query_name
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_C22.fna	/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_C22.fna	100.00	99.52	99.52	GT14BC2:NODE_2886_length_6948_cov_2.316693	GT14BC2:NODE_2886_length_6948_cov_2.316693
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S23.fna	/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_C22.fna	96.78	46.22	59.76	NODE_193_length_52978_cov_6.563883	GT14BC2:NODE_2886_length_6948_cov_2.316693
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S16.fna	/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_C22.fna	93.54	29.23	37.89	NODE_9_length_193081_cov_8.934825	GT14BC2:NODE_2886_length_6948_cov_2.316693
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S26.fna	/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S26.fna	100.00	99.77	99.77	NODE_16_length_104489_cov_9.279412	NODE_16_length_104489_cov_9.279412
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/Pseudomonas_helleri.fna	/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S26.fna	97.54	83.59	89.36	NZ_JYLD01000010.1 Pseudomonas helleri strain DSM 29165 10_220866_45.2889, whole genome shotgun sequence	NODE_16_length_104489_cov_9.279412
/data/data/Upcycling/research/additional/C5_panMICP_env/refs/HERO_S23.fna	/data/data/Upcycling/research/additional/C5_panMI
```

## C6 In-situ abundance proxy

(SPAdes contig `cov_X` length-weighted; raw FASTQ unavailable.)

```
                metric  hero_mean  rest_mean  hero_median  rest_median     MWU_p
0  length_weighted_cov  26.756106  25.545793     5.751760     9.046683  0.272327
1           median_cov  23.963773  24.645675     5.040609     8.572691  0.333487
2             mean_cov  25.546721  24.735012     5.306006     9.178483  0.314320
```

```
    source  count       mean     median        std
0   cattle     22   5.852739   4.666129   3.624958
1  poultry     42  28.631291   9.748761  47.308926
2    sheep     15  57.427073  33.271047  62.777102
3    swine     32  20.317635  11.993405  24.865754
```
