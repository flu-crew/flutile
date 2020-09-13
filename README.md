[![Build Status](https://travis-ci.org/flu-crew/flutile.svg?branch=master)](https://travis-ci.org/flu-crew/flutile)
![PyPI](https://img.shields.io/pypi/v/flutile.svg)

# flutile

## Installation

``` sh
pip install flutile
```

## Commands

### `compare`

`compare` takes a multiple-sequence alignment as an input and creates a
character difference table. This command is designed for preparing amino acid
difference tables (although it will accept DNA as well). Below is an example
output:

| site | A | B | C | D | E | F | G | H |
|------|---|---|---|---|---|---|---|---|
| 2    | L | M | M | M |   |   | M | M |
| 12   | N | I | I | I | T | T | I | T |
| 66   | T |   |   |   |   |   | K |   |
| 110  | S | P | P | P |   |   | P | P |
| 116  | T |   |   |   |   |   | A |   |
| 140  | K | R | R | R |   |   | R | R |
| 146  | K | E | E | E |   |   | E |   |
| 147  | K |   |   |   |   |   |   | M |

Sequences B-H are compared to sequence A. Every difference is reported with the
alternative amino acid, while identities are left blank.

`compare` takes two arguments (`--make-consensus` and
`--use-consensus-as-reference`) which will add a consensus column as the final
column in the table of as the first (in which case it is the reference column).

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

HA trimming is based on:
    Burke, Smith (2014) *A Recommended Numbering Scheme for Influenza A HA Subtypes*


--- | ------------------------------------------------- | ------------------- | ---------------------- |
H1  | A/United Kingdom/1/1933                           | MKARLLVLLCALAATDA   | DTICIGYHANNS           |
H2  | A/Singapore/1/1957                                | MAIIYLILLFTAVRG     | DQICIGYHANNS           |
H3  | A/Aichi/2/1968                                    | MKTIIALSYIFCLPLG    | QDLPGNDNSTATLCLGHHAVPN | * MKTIIALSYIFCL*A*LG
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
