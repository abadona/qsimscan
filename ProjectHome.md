**QSimScan** (Quick SIMilarity SCANner) is a technology for fast biosequence similarity search. It is a foundation for two software tools:

  * **NSimScan**, for searching similarities in nucleic acid sequences and
  * **PSimScan**, for searching similarities in protein sequences.

These tools perform 5 to 100 times faster than standard NCBI BLAST, depending on chosen parameters. Their sensitivity and selectivity are comparable to NCBI BLAST at the slowest settings, and decrease with the increase of speed, however, they are maintained at the levels reasonable for most tasks.

**QSimScan** has the greatest advantage when used on large collections of query sequences. Comparing entire proteome of E.coli (4,132 proteins) to the NCBI non-redundant protein database of 6,441,864 records takes 1.5 hours on moderately powerful laptop computer (over 76 hours with NCBI BLAST).

[PSimScan article](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0058505) is published in PLOS One in 2012.

For other software developed by SciDM group please see http://scidm.org/