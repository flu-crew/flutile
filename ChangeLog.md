0.12.1 [2021-07-15]
-------------------

 * Fix bug in dna2dna conversions and add tests

0.12.0 [2021-05-05]
-------------------

 * Don't spam .tmp temporary alignment files everywhere
 * Add `--count` flag to aggregate counts in aadiff tables

0.11.0 [2021-04-12]
-------------------

 * Add support for for NA subtypes and internal genes 


0.10.0 [2021-03-24]
-------------------

 * Extract motifs by parsing interval expressions such as:
   - "motif=34,56-60"
   - "34"
   - "162,166,167-190"
 * Raise a nice error message when `--subtype` is missing 

0.9.0 [2021-03-19]
------------------

 * Add `trim motif` command for extracting motifs

0.8.0 [2021-03-19]
------------------

 * Fix type in Caton82 annotations
 * Generalize trim ha1 to work for all 18 HA subtypes 

0.7.0 [2020-12-17]
------------------

`annotate` subcommand for mapping between subtypes

0.6.1 [2020-12-03]
------------------

Add wiley81 H3 reference

0.6.0 [2020-12-03]
------------------

allow aadiff annotations 

   - `--caton82` for antigenic sites
   - `--annotation-tables=STR` for list of annotation tables to load
   - `--join-annotations` to join all annotations into one column

0.5.0 [2020-09-14]
------------------

 * Fix bug in h3 trim
 * Add tests for trim
 * Add sequence samples from all h1 and h3 clades
 * Add aadiff subcommand with absolute HA indexing against subtype refs

0.4.0 [2020-09-11]

Add `trim`

0.3.2 [2020-08-13]

Fix bug in --print-groups

0.3.1 [2020-08-12]

Bug fixes

0.3.0 [2020-08-06]
------------------

Changes to represent:
 * Add --print-groups option
 * Add progress bar
 * Only consider dates if --max-day-sep argument is given (no default)
 * If --max-day-sep is given, keep the latest entry from each group, otherwise
   keep first alphabetically
 * Do not parse states unless --same-state is given 

0.2.0 [2020-01-06]
------------------

 * add group by USA state

0.1.0 [2019-12-20]
------------------

 * initial release
