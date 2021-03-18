0.8.0 [2021-03-xx]
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
