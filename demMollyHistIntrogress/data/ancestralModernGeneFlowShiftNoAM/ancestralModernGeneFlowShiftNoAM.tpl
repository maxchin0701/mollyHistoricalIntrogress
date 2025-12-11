//Number of population samples (demes)
2
//Population effective sizes (number of genes)
NPLAT
NPMEX
//Sample sizes
8
4
//Growth rates : negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 PARMIGCUR
PARMIGCUR 0
//Migration matrix 1
0 PARMIGANC
PARMIGANC 0
//Migration matrix 2
0 0
0 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
2 historical event
300000 0 0 0 1 0 1
7800000 1 0 1 NANC 0 2 absoluteResize
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 1.35e-8