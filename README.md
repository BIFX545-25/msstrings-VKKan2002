[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/j9pqPnf_)
# MSstrings


Let’s start by pulling a protein amino acid sequence from UniProt.

``` r
library(UniProt.ws)
library(Biostrings)

seqs <- queryUniProt('P51681', 'sequence')$Sequence

AAStringSet(seqs)
```

    AAStringSet object of length 1:
        width seq
    [1]   352 MDYQVSSPIYDINYYTSEPCQKINVKQIAARLLP...FCKCCSIFQQEAPERASSVYTRSTGEQEISVGL

If we are using Trypsin in our Mass Spec protocol, we’ll expect this to
be split at either lysine (K) or arginine (R), except when followed by a
proline (P). So, these are the peptide sequences we might expect to get
in our MS results:

``` r
library(stringr)

# cut points (plus the ends of the sequence)
cp <- c(0, str_locate_all(seqs, '[R|K]')[[1]][,'start'], nchar(seqs))

# possible fragments
frags <- str_sub(seqs, cp[-length(cp)] + 1, cp[-1])

# keep a few that are longer
frags <- frags[nchar(frags) > 15][1:4]
```

Next, we will add some modification metadata and format the same way it
is output from Spectronaut.

``` r
# Glycosylation
frags[1] <- str_replace(frags[1], 'S', 'S[Glycosylation (S)]')

# Oxidation
frags[2] <- str_replace(frags[2], 'M', 'M[Oxidation (M)]')

# Acetylation
frags[3] <- str_replace(frags[3], 'Y', 'Y[Acetylation (Y)]')

# add a second modification (Deamidation) to #3
frags[3] <- str_replace(frags[3], 'Q', 'Q[Deamidation (Q)]')

# add underscores as Spectronaut does
frags <- paste0('_', frags, '_')
```

This is typical of what we would see in Spectronaut output. The
`parse_mods` function in this package will parse these strings and
return a data.frame with the peptide sequences and metadata.

``` r
library(MSstrings)

parse_mods(frags)
```

    22-letter AAString object
    seq: MDYQVSSPIYDINYYTSEPCQK
              |
              Glycosylation

    28-letter AAString object
    seq: LLPPLYSLVFIFGFVGNMLVILILINCK
                          |
                          Oxidation

    64-letter AAString object
    seq: SMTDIYLLNLAISDLFFLLTVPFWAHYAAAQWDFGNTMCQLLTGLYFIGFFSGIFFIILLTIDR
              |                        |
              |                        Deamidation
              Acetylation

    28-letter AAString object
    seq: TVTFGVVTSVITWVVAVFASLPGIIFTR
      No modifications
