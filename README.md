# chipenrichfast

Faster version of the chipenrich (https://github.com/sartorlab/chipenrich) R package for gene set enrichment testing using ChIP-seq data.

Uses vectorised code and pre-filtering of genesets to cut down on run time.

## Usage
```
devtools::install_github("signalbash/chipenrichfast")

##replace input_peaks with a GRanges or data.frame with peak locations

chipenrichfast::chipenrich(peaks=input_peaks, genome='hg38', genesets="GOBP")

```
