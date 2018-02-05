#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#define MAX_ALLELES 256 // I can't fathom there being more than this

int nsamples;
int total_sites;
int pass_het;  // Number of genotypes, use this + hom to calculate # of variants GTs
int pass_hom;
int pass_ref;
int missing;
int fail_het;  // Number of genotypes, use this + hom to calculate # of variants GTs
int fail_hom;
int fail_ref;
int monomorphic;
int *var_allele_counts;
bcf_hdr_t *header;

/*
 *     This short description is used to generate the output of `bcftools plugin -l`.
 *     */
const char *about(void)
{
    printf("Prints per-variant summary stats for each variant, e.g. count of passing genotypes. Totals at bottom.");
}

/*
 *     Called once at startup, allows to initialize local variables.
 *         Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
 *         */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    total_sites = 0;
    pass_het = 0;
    pass_hom = 0;
    pass_ref = 0;
    fail_het = 0;
    fail_hom = 0;
    fail_ref = 0;
    missing = 0;
    monomorphic = 0;
    header = in;
    var_allele_counts = calloc(MAX_ALLELES, sizeof(int));
    
    printf("chr,pos,pass_homref,pass_hetvar,pass_homvar,fail_homref,fail_hetvar,fail_homvar,missing,minor_allele_freq,is_monomorphic\n");
 
    return 1;
}


/*
 *     Called for each VCF record. Return rec to output the line or NULL
 *         to suppress output.
 *         */
bcf1_t *process(bcf1_t *rec)
{
    total_sites++;

    int bcf_get_success_check;

    char **filter_data = NULL;
    int32_t num_filter_data = 0;
    bcf_get_success_check = bcf_get_format_string(header, rec, "FT", &filter_data, &num_filter_data);
 
    int32_t *gt_data = NULL;
    int32_t num_gt_data = 0;
    bcf_get_success_check = bcf_get_genotypes(header, rec, &gt_data, &num_gt_data);

    // printf("numgtdata:%d, numfiltdata:%d\n", num_gt_data, num_filter_data);

    if (bcf_get_success_check <= 0)
    {
        printf("Error getting gt");
        return NULL;
    }

    int var_het_pass = 0;
    int var_hom_pass = 0;
    int var_ref_pass = 0;
    int var_het_fail = 0;
    int var_hom_fail = 0;
    int var_ref_fail = 0;
    int var_miss = 0;
    bool is_monomorphic = false;  // all genotypes 1/1

    size_t i;
    for (i = 0; i < MAX_ALLELES; i++)
        var_allele_counts[i] = 0;

    for (i = 0; i < nsamples; i++)
    {
        int all1 = bcf_gt_allele(gt_data[2*i + 0]);
        int all2 = bcf_gt_allele(gt_data[2*i + 1]);
        bool is_pass = strcmp(filter_data[i], "PASS") == 0;
        is_pass = is_pass || (strcmp(filter_data[i], "No_var") == 0);

        if (all1 >= 0 && all1 < MAX_ALLELES)
            var_allele_counts[all1]++;
        if (all2 >= 0 && all1 < MAX_ALLELES)
            var_allele_counts[all2]++;

        if (all1 == 0 && all2 == 0)
        {
            if (is_pass)
            {
                pass_ref++;
                var_ref_pass++;
            }
            else
            {   
                fail_ref++;
                var_ref_fail++;
            }
        }
        else if (all1 != all2 && all1 >= 0 && all2 >= 0 && all1 < MAX_ALLELES && all2 < MAX_ALLELES)
        {
            if (is_pass)
            {
                pass_het++;
                var_het_pass++;
            }
            else
            {   
                fail_het++;
                var_het_fail++;
            }
        }
        else if (all1 == all2 && all1 >= 0 && all2 >= 0 && all1 < MAX_ALLELES && all2 < MAX_ALLELES)
        {
            if (is_pass)
            {
                pass_hom++;
                var_hom_pass++;
            }
            else
            {   
                fail_hom++;
                var_hom_fail++;
            }
        }
        else
        {
            missing++;
            var_miss++;
        }
    }


    int non_1_allele_sum = 0;
    for (i = 2; i < MAX_ALLELES; i++)
        non_1_allele_sum += var_allele_counts[i];

    if (var_allele_counts[0] == 0 && var_allele_counts[1] > 0 && non_1_allele_sum == 0)
    {
        is_monomorphic = true;
        monomorphic++;
    }

    // printf("chr,pos,pass_homref,pass_hetvar,pass_homvar,fail_homref,fail_hetvar,fail_homvar,missing,minor_allele_freq,is_monomorphic\n");
    printf("%s,", bcf_hdr_id2name(header, rec->rid));
    printf("%d,", rec->pos + 1);
    printf("%d,", var_ref_pass);
    printf("%d,", var_het_pass);
    printf("%d,", var_hom_pass);
    printf("%d,", var_ref_fail);
    printf("%d,", var_het_fail);
    printf("%d,", var_hom_fail);
    printf("%d,", var_miss);

    int total_alleles_observed = 0;
    for(i = 0; i < MAX_ALLELES; i++)
        total_alleles_observed += var_allele_counts[i];
    if (total_alleles_observed != 0)
        printf("%lf,", var_allele_counts[1] / (double) total_alleles_observed);
    else
        printf("0,");

    if (is_monomorphic)
        printf("True\n");
    else
        printf("False\n");

    free(gt_data);
    free(filter_data[0]);
    free(filter_data);

    return NULL;
}


/*
 *     Clean up.
 *     */
void destroy(void)
{
    printf("TOTALS:\n");
    printf("num_samples,num_variant_sites,pass_homref,pass_hetvar,pass_homvar,fail_homref,fail_hetvar,fail_homvar,"
           "het_hom_ratio,pass_het_hom_ratio,fail_het_hom_ratio,monomorphic_sites\n");
    printf("%d,%d,%d,%d,%d,%d,%d,%d,%lf,%lf,%lf,%d\n", nsamples, total_sites, pass_ref, pass_het, pass_hom, fail_ref, fail_het, fail_hom,
           (pass_het + fail_het) / (double) (pass_hom + fail_hom), pass_het / (double) pass_hom, 
           fail_het / (double) fail_hom,  monomorphic);
}
