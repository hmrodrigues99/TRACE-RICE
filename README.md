# Whole-genome polymorphisms and relatedness of rice varieties circulating in the Mediterranean market
Authors: Hugo M. Rodrigues, M. Beatriz Vieira, Pedro M. Barros, M. Margarida Oliveira
(ITQB NOVA, Oeiras Portugal)

This repository holds all scripts created and used for the elaboration of the **Whole-genome polymorphisms and relatedness of rice varieties circulating in the Mediterranean market** research article.
*Currently under revision*

Required Tool(s)

* bwa-mem (v0.7.17)
* gatk (v4.1.3.0) OR docker (20.10.21)
* samtools (v1.7)

Download Data
```
mkdir tracerice
./getVarieties.sh   # this might take a while
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz   # Oryza sativa reference genome
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-59/variation/vcf/oryza_sativa/oryza_sativa.vcf.gz   # known Oryza sativa variants
```

Mapping and processing
```
bwa index Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
./map_bwamem.sh
./sam_to_bam.sh sample_list.txt
```

Additional steps before GATK BaseRecalibrator
```
# Add sample read groups
samtools addreplacerg -r ID:Bomba -r LB:BombaLB -r SM:Bomba -r PL:ILLUMINA -o align_Bomba/Bomba_RGdedup.bam -O BAM align_Bomba/Bomba_dedup.bam
samtools addreplacerg -r ID:Puntal -r LB:PuntalLB -r SM:Puntal -r PL:ILLUMINA -o Puntal_RGdedup.bam -O BAM Puntal_dedup.bam
# Connect with GATK using Docker
docker pull broadinstitute/gatk:latest
docker run -v ~/tracerice:/gatk/data -it broadinstitute/gatk:latest
# Prepare required files
samtools faidx Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
gatk CreateSequenceDictionary -R Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
gzip -d oryza_sativa.vcf.gz
gatk IndexFeatureFile -I oryza_sativa.vcf
```

GATK BaseRecalibrator
```
./bqsr.sh sample_list.txt
```

GATK HaplotypeCaller
```
./haplotype_caller.sh
```

GATK CombineGVCFs & GenotypeGVCFs
```
./get_cohort_vcf.sh
```

Filter variants using GATK SelectVariants & VariantFiltration
```
./filter_vcf.sh
```

Required Tool(s)

* bcftools (v1.7)
* SnpEff (v5.1d)

Annotating variants with SnpEff and obtain a HIGH impact SNP table
```
java -jar snpEff.jar -c /path_to/snpEff/snpEff.config -v Oryza_sativa PASS_cohort.vcf.gz > annotated_PASS_cohort.vcf
java -jar SnpSift.jar filter "((ANN[*].IMPACT = 'HIGH'))" annotated_PASS_cohort.vcf > HIGH_PASS_cohort.vcf
bcftools view -H HIGH_PASS_cohort.vcf > HIGH_PASS_cohort.tab
```

Generating SNP density heatmaps
```
# Install necessary packages and open scripts to change working directory as needed
heatmap.R
```

QTL annotation and HIGH impact gene enrichment analysis
```
# eatingqualityQTL.tab and seedQTL.tab files retrieved from SnpEff database
# Install necessary packages and open scripts to change working directory and/or sample list as needed
snp_in_qtl.R
variant_enrichment.R
```



