---
title: headMD

---

# Analysis of Amazon molly historical introgression

This repository contains code related to three analyses published as part of (PAPER HERE). 

* Code associated with demographic inference in fastsimcoal can be found under the `demMollyHistIntrogress/` directory.
* Code assocaited with calling and analysis of shared loss of heterozygosity fragments can be found under the `SNPBasedLOH/` directory.
* Code associated with ancestral genome reconstruction in Cactus can be found under the `AGRHistIntrogress/` directory.

## Conda environments
In the `condaEnvs` directory, you will find `.yml` files associated with conda environments that are used at various points in the analyses. Be sure to build these environments by using:

```unix
conda env create -f ${ENV}.yml
```

## Additional software
In addition to the Conda environments, the following software will need to be installed/loaded as a module and available in your environment.

[bedtools](https://github.com/arq5x/bedtools2)

[GATK](https://github.com/broadinstitute/gatk)

[samtools](https://github.com/samtools/samtools)

[bcftools](https://github.com/samtools/bcftools)

[Cactus](https://github.com/ComparativeGenomicsToolkit/cactus)

[whatshap](https://github.com/whatshap/whatshap)

[fastsimcoal2](https://cmpg.unibe.ch/software/fastsimcoal28/)

## Loading software/exporting paths
At the top of scripts, you'll see software either loaded as modules,

```unix
module load ${MODULE_NAME}
```

or exported to the path,

```unix
#for scripts in SNPBasedLOH or AGRHistIntrogress
PATH=$PATH:/ocean/projects/bio230047p/mchin/software/${SOFTWARE}

#for scripts in demMollyHistIntrogress
export PATH=$PATH::/group/awhitehegrp/max/software/
```
Paths will need to be changed to the appropriate directorys where software is installed.

If you install any of the software above locally but the script loads a module to make it available, be sure to export paths as needed, and vice versa.