# taxme

[![Build Status](https://travis-ci.org/epruesse/taxme.svg?branch=master)](https://travis-ci.org/epruesse/taxme)

# Install
To get the current development version from github:
```R
# install.packages("devtools")
devtools::install_github("epruesse/taxme")
```

# Usage Examples

```R
> library(taxme)

# Get the taxonomic path for a taxID:
> ncbi.path(7, ranks="standard")
[[1]]
           rank                     name
1:      species Azorhizobium caulinodans
2:        genus             Azorhizobium
3:       family        Xanthobacteraceae
4:        order              Rhizobiales
5:        class      Alphaproteobacteria
6:       phylum           Proteobacteria
7: superkingdom                 Bacteria

# Get the classification for a taxID at a given rank:
> ncbi.group(9606, ranks="superkingdom")
[1] Eukaryota
Levels: Archaea Bacteria Eukaryota root Viroids Viruses

# Given a data.frame with read counts in numReads and NCBI tax IDs in taxID, 
# make a barplot summing the number of reads at genus level, or where the
# taxID indicates human or synthetic reads:
> library(tidyverse)
> df %>% 
    mutate(clade = ncbi.group(taxID, 
                              ranks = c("genus"), 
                              nodes = c("Homo sapiens", "synthetic construct")
           ) %>%
    group_by(clade) %>%
    summarise(numReads = sum(numReads)) %>%
    ggplot(aes(fill = clade, 
               y = numReads)) +
      geom_col()
```
