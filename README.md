# offtargets
Find CRISPRi offtargets from a .fasta file and GenBank files (optionally zipped)



# Quickstart

`python <targets.py> <queries.fasta> <target_organism.gb{.gz}> <number of mismatches allowed>`

`python targets.py ab_key.fasta GCF_009759685.1.gb 1`

```python

| q_name               | q_seq                | q_len | t_locus_tag | t_chromosome | t_seq                | diff_len | diff | coord           | offset | q_dir | t_dir |
|----------------------|----------------------|-------|-------------|--------------|----------------------|----------|------|-----------------|--------|-------|-------|
| AAAAACATGTCCACCAGTAA | AAAAACATGTCCACCAGTAA | 20    | GO593_05510 | CP046654.1   | AAAAACATGTCCACCAGTAc | 1        | c20A | 1167192-1167212 | 51     | F     | R     |
| AAAAACATGTCCACCAGTAC | AAAAACATGTCCACCAGTAC | 20    | GO593_05510 | CP046654.1   | AAAAACATGTCCACCAGTAC | 0        | -    | 1167192-1167212 | 51     | F     | R     |
| AAAAACATGTCCATCAGTAC | AAAAACATGTCCATCAGTAC | 20    | GO593_05510 | CP046654.1   | AAAAACATGTCCAcCAGTAC | 1        | c14T | 1167192-1167212 | 51     | F     | R     |
| AAAAACTTGCTGGCCGCTAC | AAAAACTTGCTGGCCGCTAC | 20    | GO593_09995 | CP046654.1   | AAAAACTTGCTGGCCGCTAC | 0        | -    | 2104826-2104846 | 161    | R     | F     |
| AAAAACTTGCTGGCCGCTCC | AAAAACTTGCTGGCCGCTCC | 20    | GO593_09995 | CP046654.1   | AAAAACTTGCTGGCCGCTaC | 1        | a19C | 2104826-2104846 | 161    | R     | F     |
| AAAAACTTGCTGGCCGCTGC | AAAAACTTGCTGGCCGCTGC | 20    | GO593_09995 | CP046654.1   | AAAAACTTGCTGGCCGCTaC | 1        | a19G | 2104826-2104846 | 161    | R     | F     |
| AAAAAGCGCAAACGATCACC | AAAAAGCGCAAACGATCACC | 20    | GO593_13345 | CP046654.1   | AAAAAtCGCAAACGATCACC | 1        | t6G  | 2774543-2774563 | 17     | F     | R     |
| AAAAATACCTGTCCTATAGT | AAAAATACCTGTCCTATAGT | 20    | GO593_02035 | CP046654.1   | AAAAATACCTGTCCTtTAGT | 1        | t16A | 412669-412689   | 182    | F     | R     |
```
