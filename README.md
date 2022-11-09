# updog - UP and DOwnstream Genetic scoring

## Requirements

The code is currently written for a high performance cluster running linux and using a Slurm job scheduler. If you are running Sun Grid Engine or something similar for job scheduling please reach out (David.Howard@kcl.ac.uk) and I can help set up UPDOG for your system.

updog requires PLINK version ≥ v1.07. The PLINK version needs to be able to read in the format of your data (bfile, pfile, or vcf). PLINK is required to be set up as an environmental variable which is called using $module load

updog requires R version ≥ 3.0.2. R is required to be set up as an environmental variable which is called using $module load

updog requires 4 sources of data:
1. Test data. This is the genetic dataset you want to make prediction in to. This data is required to be split by chromosome and in bfile, pfile or vcf format.
2. Linkage Disequilibrium (ld) reference data. This has to match the ancestry of the summary stats (typically the 1000 genomes or the HRC reference data panels. This data is required to be split by chromosome and in bfile, pfile or vcf format.
3. Genome-wide summary statistics. The summary statistics file should be genome-wide, space separated, with a header row. The first 4 columns must contain SNP Name, A1 allele, A2 allele, Effect Size. The header row is ignored along with any additional columns. The summary statistics file must use either of the following formats:

```
SNP A1 A2 beta
rs1000000 A G 0.0179
rs1000002 T C -0.0460
rs10000037 A G -0.0049
...
```

or

```
SNP A1 A2 OR
rs1000000 A G 1.0181
rs1000002 T C 0.9550
rs10000037 A G 0.9951
...
```

4. Genome-wide genetic scores for making predictions (typically created from passing the summary statistics file to packages such as PRS-CS or GCTB, or by using *P*-value thresholding and clumping). The genetic scores file should be genome-wide, space separated, with a header row. The first 6 columns must contain chromosome number, basepair position, SNP Name, A1 allele, A2 allele, Genetic Score. The header row is ignored along with any additional columns. The genetic scores file must use either of the following formats:

```
Chr BP SNP A1 A2 beta
5 124252658 rs10041789 A G 0.0262
8 61718717 rs10103011 T C 0.0298
14 104017953 rs10149470 A G -0.0280
...
```

or

```
Chr BP SNP A1 A2 OR
5 124252658 rs10041789 A G 1.0265
8 61718717 rs10103011 T C 1.0302
14 104017953 rs10149470 A G 0.9724
...
```


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

Chromosome number should be the last part of the filename prefix and omitted along with the file type. So if the full filename was UK_biobank_chr1.vcf then use `-t UK_biobank_chr`

-u [filetype of test data] `-u bfile, -u pfile, or -u vcf`

-l [location of ld reference data for summary statistics]

Chromosome number should be the last part of the filename prefix and omitted along with the file type. So if the full file name was 1000g_eur_chr1.vcf then use `-t 1000g_eur_chr`

-m [filetype of ld reference data for summary statistics] `-m bfile, -m pfile, or -m vcf`

-s [location of summary statistics] {optional flag -n}

The optional -n flag indicates that effect sizes are odds ratios else beta coefficients are assumed.

-g [location of genetic scores] {optional flag -i}

The -i flag indicates that genetic scores are odds ratios else beta coefficients are assumed. 

-p [location of plink environment module]

So if you'd normally enter `module load plink/1.9` on the command line use `-p plink/1.9`

-r [location of R environment module]

So if you'd normally enter `module load r/4.1.1` on the command line use `-r r/4.1.1`

-a "[any additional scheduler requirements]"

This flag is optional and can be used to provide the scheduler with any additional information such as for placing jobs on a particular partition or for accounting purposes. updog already requests runtime and memory and so these do not need to be requested.
 
-o [name for output]

## Example code to run updog

```
./updog \
-t ~/scratch/UK_biobank_chr \
-u bfile \
-l ~/scatch/1000g_eur_chr \
-m vcf \
-s DepressionSummaryStatistics.txt \
-g LeadVariantsExtractedfromSummaryStatistics.txt \
-p plink/1.9 \
-r r/4.1.1 \
-a "--partition=cpu" \
-o PredictionOfDepressionInUKBiobank
```

## Output

Once updog has run there should be a results file called genomewidescores_[name for output].txt in the directory. This will be a four column file with the family id, individual id, and phenotype from the test data in the first three columns. The fourth column contains the updog calculated risk score for each individual.

updog cuts the genome in to chunks and if any chunks fail to complete there will be a file called resubmitjobs_{outname} in the directory. If this is produced, enter `./resubmitjobs_[name for outname]` on the command line.

updog creates a log directory and will store additional output in there

## Support

Please send any questions about updog to: [David.Howard@kcl.ac.uk](mailto:David.Howard@kcl.ac.uk)

## Citation

Citation will appear here in due course
