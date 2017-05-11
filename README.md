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
