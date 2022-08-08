[![DOI](https://zenodo.org/badge/91003955.svg)](https://zenodo.org/badge/latestdoi/91003955)

METAL is a tool for the meta-analysis of genome-wide association studies.

Documentation is available at:
http://www.sph.umich.edu/csg/abecasis/Metal/

This version makes two changes:

1. It implements the DerSimonian-Laird random effects model. This requires a third pass through the files, calculating tau-square and the random effects parameters, so it takes about 50% longer than the heterogeneity analysis, and three times longer than the fixed effects analysis
2. The results are now all printed to 7 significant figures (previously it was only 4)

To run the random effects model it is necessary to issue the following command:

```
ANALYZE RANDOM
```

See the `examples/GlucoseExample/metal2.txt` for an example.

Note that currently there is only a compiled executable available for Linux machines

If you use this please cite using the DOI above as well as the original METAL software here https://doi.org/10.1093/bioinformatics/btq340
