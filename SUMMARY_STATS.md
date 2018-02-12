# Descriptions of HGSC BCFtools Plugins Summary Stats

## Variant Level

```
bcftools +HGSC_variant_summary input.vcf > variant_summary.tsv
```

### Individual Variant (Row)
* **chr**:  contig
* **pos**:  position
* **pass_homref**:  number of genotypes with "No_variant" or "PASS" in FT format field and "0/0" in GT format field
* **pass_hetvar**:  number of genotypes with "No_variant" or "PASS" in FT format field and e.g. "0/1", "1/2" in GT format field
* **pass_homvar**:  number of genotypes with "No_variant" or "PASS" in FT format field and e.g. "1/1", "2/2" in GT format field
* **fail_homref**:  number of genotypes without "No_variant" or "PASS" in FT format field and "0/0" in GT format field
* **fail_hetvar**:  number of genotypes without "No_variant" or "PASS" in FT format field and e.g. "0/1", "1/2" in GT format field
* **fail_homvar**:  number of genotypes without "No_variant" or "PASS" in FT format field and e.g. "1/1", "2/2" in GT format field
* **missing**:  all other genotypes, e.g. "./.", "./0"
* **minor_allele_freq**:  number of "1" alleles observed divided by all alleles observed ("."s are not included)
* **is_monomorphic**: True if all non-"." alleles are "1"

### Totals
* **num_samples**:  number of sample columns parsed from header
* **num_variant_sites**:    number of non-header rows parsed from VCF
* **pass_homref**:  number of genotypes with "No_variant" or "PASS" in FT format field and "0/0" in GT format field
* **pass_hetvar**:  number of genotypes with "No_variant" or "PASS" in FT format field and e.g. "0/1", "1/2" in GT format field
* **pass_homvar**:  number of genotypes with "No_variant" or "PASS" in FT format field and e.g. "1/1", "2/2" in GT format field
* **fail_homref**:  number of genotypes without "No_variant" or "PASS" in FT format field and "0/0" in GT format field
* **fail_hetvar**:  number of genotypes without "No_variant" or "PASS" in FT format field and e.g. "0/1", "1/2" in GT format field
* **fail_homvar**:  number of genotypes without "No_variant" or "PASS" in FT format field and e.g. "1/1", "2/2" in GT format field
* **het_hom_ratio**:  (pass_hetvar + fail_hetvar) / (pass_homvar + fail_homvar)
* **pass_het_hom_ratio**:  pass_hetvar / pass_homvar
* **fail_het_hom_ratio**:  fail_hetvar / fail_homvar
* **monomorphic_sites**:  number of rows where all non-"." alleles are "1"


## Sample Level

```
bcftools +HGSC_sample_summary input.vcf > sample_summary.tsv
bcftools +HGSC_sample_summary input.vcf -- --fail > sample_summary.tsv
bcftools +HGSC_sample_summary input.vcf -- --both > sample_summary.tsv
```

By default, only uses passing genotypes ("No_var" or "PASS" in FT format field) for *ALL* metrics, including average_coverage. With --fail option, only looks at failing genotypes. With --both, includes all genotypes.

* **sample**: sample name
* **variant_count**:  number of variant genotypes observed (sum of hetvar and homvar)
* **passing_variant_count**:  number of variant genotypes observed where "No_var" or "PASS" is in FT format field
* **ti_tv_ratio**:  number of transition variant alleles over transversion variant alleles
* **homref**: number of "0/0" genotypes
* **hetvar**: number of e.g. "0/1", "1/2" genotypes
* **homvar**: number of "1/1", "2/2", "3/3" genotypes
* **missing**:  number of all other genotypes (e.g. "./.", "./0")
* **het_hom_ratio**: hetvar / homvar
* **missing_rate**: missing / (missing + hetvar + homvar + homref)
* **average_coverage**: coverage_numerator / coverage_denominator
* **coverage_numerator**: sum of all non-"." values in DP format field (sum of all coverage)
* **coverage_denominator**: count of all non-"." values in DP format field (number of non-blank DP fields)

