---
title: "gipr_analysis"
author: "Erik Ländström"
date: "7 January 2019"
output: 
  html_document:
    keep_md : TRUE
---




```r
library(tidyverse)
library(knitr)
library(DEP)
```

```
## Warning in fun(libname, pkgname): mzR has been built against a different Rcpp version (0.12.16)
## than is installed on your system (1.0.0). This might lead to errors
## when loading mzR. If you encounter such issues, please send a report,
## including the output of sessionInfo() to the Bioc support forum at 
## https://support.bioconductor.org/. For details see also
## https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.
```



```r
# Read gipr data from maxquant function 
make_tibble_for_DEP <- function(file = "proteinGroups.txt") {
  # Library 
  library(tidyverse)
  
  # Read file
  data_raw <- read_tsv(file)
  
  # Extract LFQ column names
  lfq_col_names <- str_extract(colnames(data_raw), "LFQ.+")
  lfq_col_names <- lfq_col_names[which(!is.na(lfq_col_names))]
  
  data_raw <<- data_raw %>%
    dplyr::select(`Protein IDs`, `Majority protein IDs`,
                  `Fasta headers`, Peptides, 
                  `Razor + unique peptides`, `Unique peptides`, 
                  lfq_col_names,
                  `Only identified by site`, `Reverse`,
                  `Potential contaminant`)
}

# Read gipr data
make_tibble_for_DEP("proteinGroups_qex2.txt")

# Data features
glimpse(data_raw)
```

```
## Observations: 3,477
## Variables: 26
## $ `Protein IDs`             <chr> "CON__ENSEMBL:ENSBTAP00000024146", "...
## $ `Majority protein IDs`    <chr> "CON__ENSEMBL:ENSBTAP00000024146", "...
## $ `Fasta headers`           <chr> NA, ";", NA, "; trypsinogen isoform ...
## $ Peptides                  <int> 4, 8, 6, 6, 4, 5, 12, 17, 1, 1, 5, 1...
## $ `Razor + unique peptides` <int> 1, 3, 1, 6, 3, 1, 5, 8, 1, 1, 2, 1, ...
## $ `Unique peptides`         <int> 1, 0, 1, 6, 3, 1, 4, 3, 1, 1, 1, 1, ...
## $ `LFQ intensity L_1312`    <dbl> 0.0000e+00, 3.3687e+07, 0.0000e+00, ...
## $ `LFQ intensity L_1313`    <dbl> 0, 0, 0, 3977000000, 0, 0, 27841000,...
## $ `LFQ intensity L_1326`    <dbl> 0, 26592000, 0, 3199000000, 0, 0, 40...
## $ `LFQ intensity L_1331`    <dbl> 0, 15657000, 0, 4496600000, 0, 0, 25...
## $ `LFQ intensity L_1488`    <dbl> 0, 0, 0, 3056000000, 0, 0, 17972000,...
## $ `LFQ intensity L_1492`    <dbl> 0, 25533000, 0, 3914300000, 0, 0, 50...
## $ `LFQ intensity L_1503`    <dbl> 0, 33040000, 0, 4879900000, 0, 0, 28...
## $ `LFQ intensity L_1506`    <dbl> 0, 0, 0, 4090700000, 0, 0, 36823000,...
## $ `LFQ intensity P_1311`    <dbl> 0, 0, 0, 2857400000, 0, 1498100, 198...
## $ `LFQ intensity P_1315`    <dbl> 0.0000e+00, 2.7457e+07, 0.0000e+00, ...
## $ `LFQ intensity P_1316`    <dbl> 0.0000e+00, 0.0000e+00, 0.0000e+00, ...
## $ `LFQ intensity P_1317`    <dbl> 0.0000e+00, 2.8221e+07, 0.0000e+00, ...
## $ `LFQ intensity P_1323`    <dbl> 0, 0, 0, 3769700000, 0, 0, 27786000,...
## $ `LFQ intensity P_1497`    <dbl> 0.0000e+00, 3.3001e+07, 0.0000e+00, ...
## $ `LFQ intensity P_1501`    <dbl> 0.0000e+00, 0.0000e+00, 0.0000e+00, ...
## $ `LFQ intensity P_1502`    <dbl> 0.0000e+00, 0.0000e+00, 0.0000e+00, ...
## $ `LFQ intensity P_1505`    <dbl> 2.5194e+07, 4.2765e+07, 1.7672e+07, ...
## $ `Only identified by site` <chr> "+", NA, NA, NA, NA, NA, NA, NA, NA,...
## $ Reverse                   <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, ...
## $ `Potential contaminant`   <chr> "+", "+", "+", "+", "+", "+", "+", "...
```


```r
# Calculate contaminants and reverse hits
n_false_proteins <- data_raw %>%
    summarise(n_reverse = length(Reverse) - sum(is.na(Reverse)),
              n_contaminants = length(`Potential contaminant`) - sum(is.na(`Potential contaminant`)))

kable(n_false_proteins,
  align = "c"
)
```



 n_reverse    n_contaminants 
-----------  ----------------
    36              37       


```r
# Function
remove_contaminants_and_reverse_hits <- function(tb) {
  # Library
  library(tidyverse)
  
  data_raw <<- tb %>%
    filter(is.na(`Potential contaminant`)) %>%
    filter(is.na(Reverse))
}

# Remove contaminants and reverse hits
remove_contaminants_and_reverse_hits(data_raw)

# Check number of proteins
kable(
  length(data_raw$`Protein IDs`),
  col.names = "Number of proteins",
  align = "c"
)
```



| Number of proteins |
|:------------------:|
|        3404        |

Filtered out a total of 73 reverse and potential 
contamination proteins. 


```r
# Remove isoform from refseq identifier
data_raw <- data_raw %>%
  mutate(Protein_ID = str_extract(`Majority protein IDs`, "NP_[0-9]+|XP_[0-9]+"))

convert_NCBIIDs_to_gene_names <- function(vec) {
  # Library
  library(biomaRt)
  library(tidyverse)
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")  
  
  # Get gene names for peptide IDs
  NP <- as_tibble(getBM(attributes = c("refseq_peptide", "external_gene_name"),
                        filters = "refseq_peptide",
                        values = vec,
                        mart = ensembl
  ))
  
  NP <- NP %>%
    rename(Protein_ID = refseq_peptide)
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  XP <- as.tibble(getBM(attributes = c("refseq_peptide_predicted", "external_gene_name"),
                        filters = "refseq_peptide_predicted",
                        values = vec,
                        mart = ensembl
  ))
  
  XP <- XP %>%
    rename(Protein_ID = refseq_peptide_predicted)
  
  # Use full_join to join the two tables
  gene_names <- full_join(NP, XP, by = "Protein_ID")
  
  # Tidy data, gene names in one column and remove NAs
  gene_names <- gene_names %>%
    gather(external_gene_name.x, external_gene_name.y, key = "Source", value = "Gene_name") %>%
    dplyr::select(Protein_ID, Gene_name) %>%
    dplyr::filter(!is.na(Gene_name)) 
  
  # If there are any "" artifacts from the biomaRt search, replace them with NA
  gene_names[gene_names == ""] <- NA
  
  # Save output in a tibble
  gene_names <<- gene_names
}

# Get gene names
convert_NCBIIDs_to_gene_names(data_raw$Protein_ID)

# join gene names with data
data_raw <- full_join(data_raw, gene_names, by = "Protein_ID")
```

Make unique row names.


```r
data_unique <- make_unique(data_raw, "Gene_name", "Protein_ID")
```

## Generate a summarizedExperiment object


```r
lfq_colnums <- grep("LFQ.", colnames(data_unique))

experimental_design <- read_tsv("experimental_design.txt")
```

```
## Parsed with column specification:
## cols(
##   label = col_character(),
##   condition = col_character(),
##   replicate = col_integer()
## )
```

```r
data_se <- make_se(data_unique, lfq_colnums, experimental_design)
```

# Data exploration of

Save nonfiltered, nonimputed output.


```r
non_filtered_data <- get_df_long(data_se)
```

Plot distribution of protein detection.


```r
plot_frequency(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
plot_coverage(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
plot_numbers(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

I realized that this was not such a good idea since I only want to filter on
the treatment variable.

## Generate a summarizedExperiment object, again


```r
# Only include treatment variable in  condition column
experimental_design <- experimental_design %>% 
  mutate(condition = str_extract(condition, "[A-Z]"))
```


```r
lfq_colnums <- grep("LFQ.", colnames(data_unique))

data_se <- make_se(data_unique, lfq_colnums, experimental_design)
```

# Data exploration of  data with focus on treatment group


```r
plot_frequency(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

```r
plot_coverage(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-8-2.png)<!-- -->

```r
plot_numbers(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-8-3.png)<!-- -->

Filter data for proteins detected in at least 5 samples in at least one group.
Since the number of samples per group is not the same, we have to use fractions 
to get what we want. 55% will give at least 5 proteins in both groups.


```r
data_filt <- filter_proteins(data_se, type = "fraction", min = 0.55)

# Check the number of proteins left after filtering
kable(
  dim(data_filt)[1],
  col.names = "Number of proteins",
  align = "c",
  caption = "Proteins left after filtering",
  format = "pandoc"
)
```



Table: Proteins left after filtering

| Number of proteins |
|:------------------:|
|        2509        |


```r
plot_coverage(data_filt)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
plot_numbers(data_filt)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-9-2.png)<!-- -->

This is not what I wanted, after reading the documentation I found out that if
you use "fraction" `filter_proteins` uses this fraction on all samples. That is
not what I wanted. To do what I want, I need to filter before making a 
`SummarizedExperiment` object.


```r
data_filt <- data_unique %>%
  gather(key = "sample", value = "LFQ", `LFQ intensity L_1312`:`LFQ intensity P_1505`) %>% 
  separate(sample, into = c("treatment", "sample"), sep = "_") %>%
  group_by(name, treatment) %>%
  dplyr::filter(sum(LFQ > 0) / length(LFQ) >= 0.55) %>%
  unite(sample, treatment, sample, sep = "_") %>% 
  spread(key = sample, value = LFQ)

# Check the number of proteins left after filtering
kable(
  dim(data_filt)[1],
  col.names = "Number of proteins",
  align = "c",
  caption = "Proteins left after filtering",
  format = "pandoc"
)
```



Table: Proteins left after filtering

| Number of proteins |
|:------------------:|
|        2601        |

With this filtering we got >100 more proteins left after filtering!

## Generate a summarizedExperiment object


```r
lfq_colnums <- grep("LFQ.", colnames(data_filt))

experimental_design <- read_tsv("experimental_design.txt")

data_se <- make_se(data_filt, lfq_colnums, experimental_design)
```

# Data exploration for the third time


```r
plot_frequency(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
plot_coverage(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

```r
plot_numbers(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-12-3.png)<!-- -->

## Normalization with vsn


```r
data_norm <- normalize_vsn(data_se)

plot_normalization(data_se, data_norm)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

```r
meanSdPlot(data_se)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-13-2.png)<!-- -->

```r
meanSdPlot(data_norm)
```

![](gipr_analysis_files/figure-html/unnamed-chunk-13-3.png)<!-- -->

When doing normalization after filtering it seems that more changes than before
filtering. I should look up the function!
