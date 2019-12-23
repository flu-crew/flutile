[![Build Status](https://travis-ci.org/arendsee/flutile.svg?branch=master)](https://travis-ci.org/arendsee/flutile)

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
