# mergeTRA
R functions to read and merge translocations (TRA) called by different structural variations caller.

Inspired by [intansv](http://venyao.github.io/intansv/).

For SV annotation, one can annotate chr1:pos1 and chr2:pos2 speratedly using svAnnotation function from [intansv](http://venyao.github.io/intansv/).



## Usage
Download and souce the scripts.
Remember to load [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) and [plyr](https://cran.r-project.org/web/packages/plyr/index.html) before sourcing.
```R
library(GenomicRanges)
library(plyr)
source("readBreakDancerTRA.R")
source("readDellyTRA.R")
source("methodsMergeTRA.R")
```

### windowsize?
Windowsize should generally be in the range of the read length to 2x read length. e.g. For PE150, windowsize is about 150-300bp.



## Questions
If you have any comments, questions or problems encounterd, please feel free to open an issue.
