Contents:

 * `caton82.txt` - H1 antigenic sites from (Caton 82)

 * `wiley81.txt` - H3 antigenic sites from (Wiley 81)

 * `subtype-refs.faa` and `subtype-refs.fna` - protein and DNA fasta files
   containing the subtype reference strains used in (Burke 2014)

 * `burke2014-index-map.txt` - A mapping between subtween subtype indices from
   supplementary table 2 in (Burke 2014). Columns for strains other than those
   in `subtype-refs.*` files were dropped. The header with the strain names is
   also dropped. The ith column referes to the ith subtype (column 1 is H1).
   Here is a list of the strains:

    A/United/Kingdom/1/1933
    A/Singapore/1/1957
    A/Aichi/2/68
    A/swine/Ontario/01911-1/1999
    A/Vietnam/1203/2004
    A/chicken/Taiwan/0705/1999
    A/Netherlands/219/03
    A/turkey/Ontario/6118/1968
    A/Swine/Hong Kong/9/98
    A/mallard/Bavaria/3/2006
    A/duck/England/1/1956
    A/duck/Alberta/60/1976
    A/gull/Maryland/704/1977
    A/mallard/Astrakhan/263/1982
    A/duck/Australia/341/1983
    A/black-headedgull/Turkmenistan/13/1976
    A/little-yellow-shoulderedbat/Guatemala/060/2010
    A/flat-faced/bat/Peru/033/2010

 * NA subtype refs - chosen based on existence of good crystal structures

  | subtype | pdb  | strain
  +---------+------+----------------
  | H1N1    |      | A/WSN/33
  | H2N2    | 1NN2 | A/Tokyo/3/1967
  | H2N3    | 4HZV | A/swine/Missouri/2124514/2006
  | H10N4   | 2HTV | A/Mink/Sweden/3900/1984
  | H12N5   | 3SAL | A/duck/Alberta/60/1976
  | H3N6    | 4QN4 | A/chicken/Nanchang/7-010/2000
  | H10N7   | 4QN3 | A/mallard/ALB/196/1996
  | H3N8    | 2HT5 | A/duck/Ukraine/1/1963
  | H11N9   | 7NN9 | A/tern/Australia/G70C/1975
  | H17N10  | 4FVK | A/little_yellow-shouldered_bat/Guatemala/153/2009
  | H18N11  | 4K3Y | A/flat-faced_bat/Peru/033/2010

  Here are some alternatives I considered:

    H1N1 3NSS A/California/04/2009
    H5N1 2HTY A/Viet_Nam/1204/2004
    N1   6LXI A/Brevig_Mission/1/1918
    N2   4GZX A/Tanzania/205/2010
    N2   5HUK A/Northern_pintail/Washington/40964/2014
    N2   1IVE A/Tokyo/3/1967
    N8   5HUF A/gyrfalcon/Washington/41088-6/2014
    N9   7NN9 A/tern/Australia/G70C/1975
