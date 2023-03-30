# offtargets
Find CRISPRi offtargets from a .txt file and GenBank files (optionally zipped)



# Quickstart

`python find_offtargets.py mismatches.txt ./*gb.gz | tabulate -1`

```python

query_seq             off_target_seq        record_id      position    score  in_gene    gene_name    locus_tag        differences             orientation
--------------------  --------------------  -----------  ----------  -------  ---------  -----------  ---------------  ----------------------  -------------
CGATCAATCAGGTCTTTAAG  CGATCAATCAGGGTTTTCTG  CP023715.1      1878143       16  True       nifU         ZMO1_ZMO1833     T13G, C14T, A18C, A19T  reverse
CGATCAATCAGGGCTTCAAG  CGATCAATCAGGGTTTTCTG  CP023715.1      1878143       16  True       nifU         ZMO1_ZMO1833     C14T, C17T, A18C, A19T  reverse
CGATCAATCAGGGCTTCAAG  CGATCAATTACGGGTTCAAT  CP023715.1      1758185       16  True                    ZMO1_ZMO1705     C9T, G11C, C14G, G20T   reverse
CGATCAATCAGCTCTTCAAG  CATTCAATCAGTTCTTCAAT  CP023715.1      1737107       16  True       purA         ZMO1_ZMO1687     G2A, A3T, C12T, G20T    reverse
CGATCAATCAGCTCTTCAAG  CGATCAAATTGCTCATCAAG  CP023715.1      1634040       16  True                    ZMO1_ZMO1593     T8A, C9T, A10T, T15A    reverse
CGATCAATCAGGTCTGCAAG  CGATAAATCAGGACAGAAAG  CP023715.1      1608057       16  True                    ZMO1_ZMO1573     C5A, T13A, T15A, C17A   reverse
CGATCAATCAGCTCTTCAAG  CGGTCAATCAGCTTTTCACG  CP023715.1       555814       17  True       nusA         ZMO1_ZMO0556     A3G, C14T, A19C         forward
CGATCAATCCGGTCTTCAAG  TGAAGAATCCGGCCTTCAAG  CP023715.1      1402397       16  True       ggt          ZMO1_ZMO1388     C1T, T4A, C5G, T13C     reverse
CGATCAATCAGGGCTTCAAG  AGATCAAACAGAGCTTCTAG  CP023715.1       657015       16  True       purB         ZMO1_ZMO0662     C1A, T8A, G12A, A18T    forward
CGATCAATCAGGTCTTTAAG  CTATCAATCAGGTCATTTTG  CP023715.1      1323884       16  True       fumA         ZMO1_ZMO1307     G2T, T15A, A18T, A19T   reverse
CGATCAATCAGGTCTGCAAG  CGATCAATCAGGTCTTCAAG  CP023715.1      1254050       19  True       recJ         ZMO1_ZMO1231     G16T                    reverse
CGATCAATCAGGTCTTTAAG  CGATCAATCAGGTCTTCAAG  CP023715.1      1254050       19  True       recJ         ZMO1_ZMO1231     T17C                    reverse
CGATCAATCAGGGCTTCAAG  CGATCAATCAGGTCTTCAAG  CP023715.1      1254050       19  True       recJ         ZMO1_ZMO1231     G13T                    reverse
CGATCAATCAGCTCTTCAAG  CGATCAATCAGGTCTTCAAG  CP023715.1      1254050       19  True       recJ         ZMO1_ZMO1231     C12G                    reverse
CGATCAATCCGGTCTTCAAG  CGATCAATCAGGTCTTCAAG  CP023715.1      1254050       19  True       recJ         ZMO1_ZMO1231     C10A                    reverse
CTCCACTAAACCCCGAATTA  CTCCACTAAACCCCGAATTT  CP023715.1      1254022       19  True       recJ         ZMO1_ZMO1231     A20T                    reverse
CTCCACTAAACCCGGAATTT  CTCCACTAAACCCCGAATTT  CP023715.1      1254022       19  True       recJ         ZMO1_ZMO1231     G14C                    reverse
CGATCAATCAGGTCTTTAAG  GGATCAATCTGGTCTTTCCG  CP023715.1       840976       16  True       ddl          ZMO1_ZMO0834     C1G, A10T, A18C, A19C   reverse
CGATCAATCAGCTCTTCAAG  CCTTCAATCCGCTCTTCATG  CP023715.1       817725       16  True       defA         ZMO1_ZMO0813     G2C, A3T, A10C, A19T    reverse
CGATCAATCCGGTCTTCAAG  CCTTCAATCCGCTCTTCATG  CP023715.1       817725       16  True       defA         ZMO1_ZMO0813     G2C, A3T, G12C, A19T    reverse
CGATCAATCAGGTCTTTAAG  GGATCAGTCATGTCTTTAAT  CP023715.1      1241508       16  True                    ZMO1_ZMO1214     C1G, A7G, G11T, G20T    forward
CGATCAATCCGGTCTTCAAG  TTTTCAATCCGGTCTTCCAG  CP023715.1       778090       16  True                    ZMO1_ZMO0780     C1T, G2T, A3T, A18C     reverse
CTCCACTAAACCCGGAATTT  ATCAAATAAACCCGCAATTT  CP023715.1      1379892       16  False                                    C1A, C4A, C6A, G15C     forward
CGATCAATCAGGTCTGCAAG  CGATCAATCAGATCAGCAAC  CP023715.1       654604       17  True       dnaK         ZMO1_ZMO0660     G12A, T15A, G20C        reverse
CGATCAATCAGCTCTTCAAG  CGATCAATCAGATCAGCAAC  CP023715.1       654604       16  True       dnaK         ZMO1_ZMO0660     C12A, T15A, T16G, G20C  reverse
CGATCAATCAGGTCTTTAAG  AGATCACGCAGGTCATTAAG  CP023715.1       624929       16  True       fliC         ZMO1_ZMO0629     C1A, A7C, T8G, T15A     reverse
CGATCAATCAGGTCTGCAAG  CGATCAATAATGTCAGAAAG  CP023715.1       491515       16  True                    ZMO1_ZMO0495     C9A, G11T, T15A, C17A   reverse
CGATCAATCAGGGCTTCAAG  CGATCAACCAAGGCTTGAGG  CP023715.1       438524       16  True                    ZMO1_ZMO0441     T8C, G11A, C17G, A19G   reverse
CGATCAATCAGGGCTTCAAG  AGGTCAACCAGGGTTTCAAG  CP023716.1        18137       16  True                    ZMO1_ZMOp32x017  C1A, A3G, T8C, C14T     forward
CGATCAATCAGGGCTTCAAG  AGGTCAACCAGGGTTTCAAG  CP023717.1        23639       16  True                    ZMO1_ZMOp33x025  C1A, A3G, T8C, C14T     reverse
```
