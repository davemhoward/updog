# updog
UP and DOwnstream Genetic scoring

Requirements:

The code is currently set to work on a high performance cluster running linux and using a Slurm job scheduler. If you are running Sun Grid Engine or something similar for job scheduling please reach out (David.Howard@kcl.ac.uk) and I can help set up UPDOG for your system.

UPDOG requires PLINK version ≥ v1.07. The PLINK version needs to be able to read in the format of your data (bfile, pfile, or vcf). PLINK is required to be set up as an environmental variable which is called using $module load

UPDOG requires R version ≥ 3.0.2. R is required to be set up as an environmental variable which is called using $module load
