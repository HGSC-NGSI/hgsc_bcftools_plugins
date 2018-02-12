## Filtering and Formatting

This document describes some typical post-processing we might perform on VCF files. BCFtools is fast, requires strict and consistent formatting, and gives us a nice interface to write our own small C plugins for interacting with VCFs, so we make frequent use of it.

Requires BCFtools 1.6.

Make sure plugins environment variable is set:

```
export BCFTOOLS_PLUGINS=/path/to/bcftools/plugins
```

### 1. Filtering

```
bcftools +HGSC_filt_w_dotdots input.vcf $MIN_DEPTH > filtered.vcf
```

This will change the "FT" subfield of each sample"s FORMAT field in the
following ways:

* If the "FT" subfield is not empty, do nothing
* Else write "PASS" if the "GT" subfield has a variant allele
* Else if the "GT" subfield is "./0", write "."
* Else write "No_var" if the "GT" subfield is "0/0" and "DP" is >=
$MIN_DEPTH
* Else write "No_data" if the "GT" subfield is "0/0" and "DP" is <
$MIN_DEPTH

This step will also change the "GT" subfield to "./." if the "FT"
subfield was marked as "No_data."

Just let MIN_DEPTH=1 for this.

### 2. Append Full Set Summary to INFO Field

```
bcftools +HGSC_append_gtcounts filtered.vcf > filtered_summarized.vcf
```

This appends a number of INFO subfields to the header, e.g.:

```
##INFO=<ID=PASS_33_15621,Number=.,Type=Integer,Description="Count of samples w/filter PASS and genotype 3/3 in 15621-sample set">
##INFO=<ID=FAIL_.._15621,Number=.,Type=Integer,Description="Count of samples w/filter FAIL and genotype ./. in 15621-sample set">
```

These subfields take the form of FILTER_GENOTYPE_NUMSAMPLES.

Possbile filter values:
* PASS
* NDAT: "No_data"
* NVAR: "No_var"
* NFLT: "."
* FAIL: Any other value

The "NUMSAMPLES" part allows us to distinguish between runs of this
plugin on the full set versus runs on the subset.

If the count for any of these combinations of filter and genotype is
nonzero, the INFO field will be updated to include the counts for that
combination.

#### **BEWARE!** 
The HGSC_append_gtcounts plugin assumes SNPs (i.e., it assumes only 4 possible alleles)! It will die if there is a row in the input VCF with 5 or more alleles.

### 3. Sample Renaming

```
bcftools reheader -h $PATH_TO_HEADER filtered_summarized.vcf > filtered_summarized_reheadered.vcf
```

We frequently need to rename samples at some point to match other components of a dataset.

Just take the old header and replace the sample names with what they
should be (in the order of the original sample names).

We also add descriptions of filters to the header in this step if they aren't already there:

```
##FILTER=<ID=low_snpqual,Description="SNP posterior probability is less than 0.50">
##FILTER=<ID=low_VariantReads,Description="Variant read depth is less than 2">
##FILTER=<ID=low_VariantRatio,Description="Variant read ratio is less than 0.25">
##FILTER=<ID=low_coverage,Description="Total coverage is less than 6">
##FILTER=<ID=high_coverage,Description="Total coverage is greater than 8000">
##FILTER=<ID=single_strand,Description="All variant reads are in a single strand direction">
##FILTER=<ID=No_data,Description="No valid reads on this site">
##FILTER=<ID=No_var,Description="No valid variants reads on this site">
```

### 4. Biallelic/Multiallelic Split

Specify minimum/maximum alleles to bcftools view:

```
bcftools view -m2 -M2 filtered_summarized_reheadered.vcf > biallelic.vcf
bcftools view -m3 filtered_summarized_reheadered.vcf > multiallelic.vcf
```

### 5. Subsetting

With "-s", use comma-delimited list of sample names:

```
bcftools view filtered_summarized_reheadered.vcf -s sample_name1,sample_name2,sample_name3 > subset.vcf
```

or use a newline-delimited list in a file with "-S":

```
bcftools view input.vcf -S $PATH_TO_LIST_OF_SAMPLES
```

### 6. Append Subset Summary to INFO Field

```
bcftools +HGSC_append_gtcounts subset.vcf > subset_summarized.vcf
```

### 7. Biallelic/Multiallelic Split

```
bcftools view -m2 -M2 subset_summarized.vcf > subset_biallelic.vcf
bcftools view -m3 subset_summarized.vcf > subset_multiallelic.vcf
```

### 8. Mappability Annotation

```
bcftools annotate input.vcf -a mappability.bed.gz -c 'CHROM,FROM,TO,-,MS' -h mappability.hdr > annotated.vcf
```

### 9. Split Multiallelic into Separate Lines

```
bcftools norm -m- multiallelic.vcf > split_multiallelic.vcf
```
