# TRACE-RICE
Collection of scripts used for the elaboration of the "TRACE-RICE" article.

Required Tool(s)

* bcftools (v1.7)
* bwa-mem (v0.7.17)
* docker (20.10.21) OR gatk (v4.1.3.0)
* samtools (v1.7)
* SnpEff (v5.1d)

To run the script pipeline, do the following steps:

Download Data
```
mkdir trace_testing
# Example varieties
wget Bomba
wget Puntal
# Oryza sativa genome and known variants
wget Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
wget oryza_sativa.vcf.gz
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
docker run -v ~/trace_testing:/gatk/data -it broadinstitute/gatk:latest
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

Annotating variants with SnpEff and obtain a HIGH impact SNP table
```
java -jar snpEff.jar -c /path_to/snpEff/snpEff.config -v Oryza_sativa PASS_cohort.vcf.gz > annotated_PASS_cohort.vcf
java -jar SnpSift.jar filter "((ANN[*].IMPACT = 'HIGH'))" annotated_PASS_cohort.vcf > HIGH_PASS_cohort.vcf
bcftools view -H HIGH_PASS_cohort.vcf > HIGH_PASS_cohort.tab
```





