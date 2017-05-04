# spector 0.6.6

* use purrr for a speed boost and memory efficiency 
* return `tbl_df` contains `chrom`, `start`, `end`, and `las`  for standard run
* added normalisation of bam file reads to 10^9

# spector 0.6.5

## User facing
* regions can now be passed as `tbl_df` allows for more flexibility
* include silent run option
  * especially useful for quick runs and large cue submissions
* better error messages
* better error reporting if empty vectors passed 

## Package internal

* more robust unit tests
* remove documentation from functions not exported to NAMESPACE

# spector 0.6.0

* added concise documentation
* include vignette
* added path independent unit tests
* included sample data available to user via `spector_sample()`