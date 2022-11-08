# updog - UP and DOwnstream Genetic scoring

## Requirements

The code is currently written for a high performance cluster running linux and using a Slurm job scheduler. If you are running Sun Grid Engine or something similar for job scheduling please reach out (David.Howard@kcl.ac.uk) and I can help set up UPDOG for your system.

updog requires PLINK version ≥ v1.07. The PLINK version needs to be able to read in the format of your data (bfile, pfile, or vcf). PLINK is required to be set up as an environmental variable which is called using $module load

updog requires R version ≥ 3.0.2. R is required to be set up as an environmental variable which is called using $module load

updog requires 4 sources of data:
1. Test data. This is the genetic dataset you want to make prediction in to. This data is required to be split by chromosome and in bfile, pfile or vcf format.
2. Linkage Disequilibrium (ld) reference data. This has to match the ancestry of the summary stats (typically the 1000 genomes or the HRC reference data panels. This data is required to be split by chromosome and in bfile, pfile or vcf format.
3. Genome-wide summary statistics
4. Genome-wide genetic scores for making predictions (typically created from passing the summary statistics file to packages such as PRS-CS or GCTB, or by using *P*-value thresholding and clumping)


## Getting Started

Clone this repository using the following git command:
```
git clone https://github.com/davemhoward/updog.git
```

To run updog type 
```
./updog
```
and supply the following flags and arguments

-t [location of test data]

Test data should be split by chromosome. Chromosome number should be the last part of the filename prefix and omitted along with the file type. So if the full filename was UK_biobank_chr1.vcf then type `-t UK_biobank_chr`

-u [filetype of test data] `-u bfile, -u pfile, or -u vcf`

-l [location of ld reference data for summary statistics]

LD reference data should be split by chromosome, chromosome number should be the last part of the filename prefix and omitted along with the file type. So if the full file name was 1000g_eur_chr1.vcf then type `-t 1000g_eur_chr`

-m [filetype of ld reference data for summary statistics] `-m bfile, -m pfile, or -m vcf`


  -s <location of summary statistics> {optional flag -n}
     File should be genome-wide, space separated, with a header row. The first 4
     columns must contain SNP Name, A1 allele, A2 allele, Effect Size. The -n
     flag indicates that effect sizes are odds ratios else beta coefficients are
     assumed. The header row is ignored along with any additional columns.
  -g <location of genetic scores> {optional flag -i}
     File should be genome-wide, space separated, with a header row. The first 6
     columns must contain chromosome number, basepair position, SNP Name, A1
     allele, A2 allele, Genetic Score. The -i flag indicates that genetic scores
     are odds ratios else beta coefficients are assumed. The header row is
     ignored along with any additional columns
  -p <location of plink environment module>
     i.e. what you enter after entering module load on the command line
  -r <location of R environment module>
     i.e. what you enter after entering module load on the command line
     Use R version ≥ 3.0.2
  -o <name for output>


