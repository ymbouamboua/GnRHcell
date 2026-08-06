# Default GnRH developmental marker programs

Returns the curated developmental gene modules used by
[`gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md)
to classify candidate GnRH markers into major developmental programs.

## Usage

``` r
gnrh_stage_modules()
```

## Value

A named list containing four character vectors:

- identity:

  Developmental identity and specification genes.

- migrating:

  Migration and guidance-associated genes.

- mature:

  Markers of differentiated GnRH neurons.

- secreting:

  Genes involved in secretory function.

## Details

The modules represent four key biological stages of the GnRH neuronal
lineage:

- **identity**: specification and developmental identity genes.

- **migrating**: genes involved in neuronal migration and axon guidance.

- **mature**: markers associated with differentiated GnRH neurons.

- **secreting**: genes involved in neuropeptide processing and
  secretion.

All genes are returned in uppercase to facilitate cross-species
comparisons and marker matching.

These modules were manually curated from the GnRH developmental
literature and are intended for exploratory marker discovery rather than
strict cell-state annotation.

Users may modify the returned gene sets and provide a customized module
list to
[`gnrh_marker_programs()`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md).

## References

Wray S. Development of gonadotropin-releasing hormone-1 neurons. Front
Neuroendocrinol. 2010.

Schwanzel-Fukuda M, Pfaff DW. Origin of luteinizing hormone-releasing
hormone neurons. Nature. 1989.

Stevenson EL, Corella KM, Chung WCJ. Ontogenesis of
gonadotropin-releasing hormone neurons. Front Endocrinol. 2013.

Cho HJ, Shan Y, Whittington NC, Wray S. Nasal placode development, GnRH
neuronal migration and Kallmann syndrome. Front Cell Dev Biol. 2019.

Taroc EZM, Prasad A, Lin JM, Forni PE. GnRH-1 neural migration from the
nose to the brain is independent from Slit2, Robo3 and NELL2. Front Cell
Neurosci. 2019.

Li Q et al. Expression of genes for kisspeptin, neurokinin B and
dynorphin in the hypothalamus. Front Endocrinol. 2020.

## See also

[`gnrh_marker_programs`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_marker_programs.md)
[`gnrh_stage_gene_references`](https://ymbouamboua.github.io/GnRHcell/reference/gnrh_stage_gene_references.md)

## Examples

``` r
modules <- gnrh_stage_modules()

names(modules)
#> [1] "identity"  "migrating" "mature"    "secreting"

modules$identity
#>  [1] "FEZF1"  "SOX2"   "HES1"   "OTX2"   "SIX3"   "SIX6"   "ISL1"   "DLX1"  
#>  [9] "DLX2"   "DLX5"   "DLX6"   "ARX"    "FOXG1"  "RAX"    "ISL2"   "MYT1"  
#> [17] "KLF7"   "ZIC2"   "ZIC4"   "ZIC5"   "ZNF483" "EBF3"   "MEIS1" 

length(modules$migrating)
#> [1] 27

custom_modules <- modules
custom_modules$identity <- unique(
  c(custom_modules$identity, "FEZF2")
)
```
