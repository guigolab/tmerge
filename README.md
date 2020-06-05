# tmerge
Merge transcriptome read-to-genome alignments into non-redundant transcript models.

tmerge compares transcript structures (or read-to-genome alignments) present in the input and attempts to reduce transcript redundancy, i.e., merge compatible input transcripts into non-redundant transcript models. The program treats spliced and monoexonic reads separately (i.e., those are never merged together).

tmerge is fast and can typically process several millions of aligned long reads in a few minutes.

## Authors
Julien Lagarde, CRG, Barcelona, contact julienlag@gmail.com
Jacob Windsor, contact me@jcbwndsr.com