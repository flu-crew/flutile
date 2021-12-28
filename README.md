![build status](https://github.com/flu-crew/flutile/actions/workflows/python-app.yml/badge.svg)
![PyPI](https://img.shields.io/pypi/v/flutile.svg)

# flutile

## Installation

``` sh
pip install flutile
```

## Commands

### `aadiff`

`aadiff` takes a multiple-sequence alignment as an input and creates a
character difference table. This command is designed for preparing amino acid
difference tables. Below is an example of a comparison of 4 H1 sequences.

`flutile aadiff --subtype=H1 mufile.faa`

| site  | A02479030 | SD0246 | SD0272 | SD0136 |
| ----  | --------- | ------ | -------| ------ |
| -3    | A         |        | T      |        |
| -2    | N         |        | S      |        |
| -1    | A         | T      |        |        |
| 1     | D         |        |        |        |
| 154   | K         |        |        |        |
| 154+1 | -         |        |        | X      |
| 154+2 | -         |        |        | X      |
| 155   | D         | S      |        | X      |
| 156   | D         | G      | N      | X      |

The `--subtype=H1` argument tells `flutile` to align the inputs against an H1
reference (A/United Kingdom/1/1933). The reference is used to determine
relative indices (the `sites` column). The index reference is used only for
indexing and does not appear in the final table. The first three rows (sites
-3, -2, -1) align to the three residues at the end of the signal peptide. Site
1 is the first residue in the mature peptide. Any gaps in the reference
alignment are indexed as `<ref_id>+<offset>`, for example 154+1 and 154+2 are
positions 1 and 2 residues after the reference position 154.

`flutile` uses the references from (Burke 2014):

--- | ------------------------------------------------- | ------------------- | ---------------------- |
H1  | A/United Kingdom/1/1933                           | MKARLLVLLCALAATDA   | DTICIGYHANNS           |
H2  | A/Singapore/1/1957                                | MAIIYLILLFTAVRG     | DQICIGYHANNS           |
H3  | A/Aichi/2/1968                                    | MKTIIALSYIFCLPLG    | QDLPGNDNSTATLCLGHHAVPN |
H4  | A/swine/Ontario/01911–2/1999                      | MLSIAILFLLIAEGSS    | QNYTGNPVICLGHHAVSN     |
H5  | A/Vietnam/1203/2004                               | MEKIVLLFAIVSLVKS    | DQICIGYHANNS           |
H6  | A/chicken/Taiwan/0705/1999                        | MIAIIVIATLAAAGKS    | DKICIGYHANNS           |
H7  | A/Netherlands/219/2003                            | MNTQILVFALVASIPTNA  | DKICLGHHAVSN           |
H8  | A/turkey/Ontario/6118/1968                        | MEKFIAIAMLLASTNA    | YDRICIGYQSNNS          |
H9  | A/swine/Hong Kong/9/1998                          | MEAASLITILLVVTASNA  | DKICIGYQSTNS           |
H10 | A/mallard/bavaria/3/2006                          | MYKIVVIIALLGAVKG    | LDKICLGHHAVAN          |
H11 | A/duck/England/1/1956                             | MEKTLLFAAIFLCVKA    | DEICIGYLSNNS           |
H12 | A/duck/Alberta/60/1976                            | MEKFIILSTVLAASFA    | YDKICIGYQTNNS          |
H13 | A/gull/Maryland/704/1977                          | MALNVIATLTLISVCVHA  | DRICVGYLSTNS           |
H14 | A/mallard/Astrakhan/263/1982                      | MIALILVALALSHTAYS   | QITNGTTGNPIICLGHHAVEN  |
H15 | A/duck/Australia/341/1983                         | MNTQIIVILVLGLSMVRS  | DKICLGHHAVAN           |
H16 | A/black-headed-gull/Turkmenistan/13/1976          | MMIKVLYFLIIVLGRYSKA | DKICIGYLSNNS           |
H17 | A/little-yellow-shouldered bat/Guatemala/060/2010 | MELIILLILLNPYTFVLG  | DRICIGYQANQN           |
H18 | A/flat-faced bat/Peru/033/2010                    | MITILILVLPIVVG      | DQICIGYHSNNS           |

 * The H3 signal peptide appears to actually be `MKTIIALSYIFCLALG`

### `represent`

`represent` takes a multiple-sequence alignment as input and removes entries
that are similar in sequence and time. The function requires that headers have
a date term (with format year/month/day). For example:

```
>A|1990-01-02
GATACA
>B|1990-02-02
CATATA
```

There may be gaps in the alignment. Sequences in the alignment that are
separated by less than or equal to `--max-day-sep` and that are a sequence
identity of greather than or equal to `--min-pident-sep` will be clustered
together. A single representative is sampled from each cluster (the latest one
with ties resolved by order).


### `trim`

Extract the HA1 regions from (currently) H1 and H3 HA proteins.

### References

  1. Burke, Smith (2014) *A Recommended Numbering Scheme for Influenza A HA Subtypes*
