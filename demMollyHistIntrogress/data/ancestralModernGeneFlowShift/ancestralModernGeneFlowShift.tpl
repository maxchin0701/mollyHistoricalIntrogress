//Number of population samples (demes)
4
//Population effective sizes (number of genes)
NPLAT
NPFORHLAT
NPFORHMEX
NPMEX
//Sample sizes
8
19
19
4
//Growth rates : negative growth implies population expansion
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
3
//Migration matrix 0
0 0 0 PARMIGCUR
PFORMIGCUR 0 0 0
0 0 0 PFORMIGCUR
PARMIGCUR 0 0 0
//Migration matrix 1
0 0 0 PARMIGANC
0 0 0 0
0 0 0 0
PARMIGANC 0 0 0
//Migration matrix 2
0 0 0 0
0 0 0 0
0 0 0 0
0 0 0 0
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix
5 historical event
300000 1 1 0 1 0 1 absoluteResize
300000 2 2 0 1 0 1 absoluteResize
300000 1 0 1 1 0 1
300000 2 3 1 1 0 1
7800000 3 0 1 NANC 0 2 absoluteResize
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of linkage blocks
1
//per Block: data type, num loci, rec. rate and mut rate + optional parameters
FREQ 1 0 1.35e-8