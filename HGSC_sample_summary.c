#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

#define MAX_COVERAGE 1000   //coverage values above this won't be used for calculating the average
#define MIN_COVERAGE 0      //coverage values below this won't be used for calculating the average

typedef struct
{
    int genotypes_with_depth;
    int passing_variants;
    int ref;  // 0/0
    int het;  // 0/1, 0/2, etc.
    int var;  // 1/1, 1/2, 2/2, etc.
    int missing;  // ./., ./0, 1/., etc.
    int transitions;
    int transversions;
    long total_coverage;  // divide this by genotypes_with_depth sites to get average coverage
} bucket_t;

typedef struct
{
    bool use_pass;
    bool use_fail;
    bool is_indel_file;
} args_t;

int nsamples;
int num_sites;
bcf_hdr_t *header;
bucket_t **samp_buckets;
args_t *args;


/*
 *     This short description is used to generate the output of `bcftools plugin -l`.
 *     */
const char *about(void)
{
    return "Calculate per-sample summary metrics.\n";
}

const char *usage(void)
{
    return "Usage: bcftools +plugin_name GENERAL_OPTIONS INPUT.bcf -- PLUGIN_OPTIONS\n"
           "About:\n"
           "\tCalculate per-sample summary metrics.\n"
           "\tDefault is to only include data where FT is 'PASS' or 'No_var'. --fail to use only failed data or --both to use both.\n"
           "\tDefault also assumes file is SNP only. You can give --indel to turn off ti/tv counts.\n" 
           "bcftools +plugin_name GENERAL_OPTIONS INPUT.bcf -- PLUGINS_OPTIONS\n"
           "bcftools +HGSC_sample_summary INPUT.bcf\n"
           "bcftools +HGSC_sample_summary INPUT.bcf -- --fail\n"
           "bcftools +HGSC_sample_summary INPUT.bcf -- --both --indel\n";
}

/*
 *     Called once at startup, allows to initialize local variables.
 *         Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
 *         */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    header = in;
    args = calloc(1, sizeof(args_t));
    args->use_pass = true;
    args->use_fail = false;
    args->is_indel_file = false;

    static struct option long_options[] =
    {
        {"help", no_argument, NULL, 'h'},
        {"fail", no_argument, NULL, 'f'},
        {"both", no_argument, NULL, 'b'},
        {"indel", no_argument, NULL, 'i'}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "hfbi", long_options, NULL)) >= 0)
    {
        switch(opt)
        {
            case 'f': args->use_pass = false; args->use_fail = true; break;
            case 'b': args->use_fail = true; break;
            case 'i': args->is_indel_file = true; break;
            default: error("%s", usage()); break;
        }
    } 

    num_sites = 0;
    
    samp_buckets = calloc(nsamples, sizeof(bucket_t *));
    int i;
    for (i = 0; i < nsamples; i++)
    {
        samp_buckets[i] = calloc(nsamples, sizeof(bucket_t));
    }

    return 1;
}


/*
 *     Called for each VCF record. Return rec to output the line or NULL
 *         to suppress output.
 *         */
bcf1_t *process(bcf1_t *rec)
{
    int bcf_get_success_check;

    num_sites++;

    char **filter_data = NULL;
    int32_t num_filter_data = 0;
    bcf_get_success_check = bcf_get_format_string(header, rec, "FT", &filter_data, &num_filter_data);
    
    int32_t *gt_data = NULL;
    int32_t num_gt_data = 0;
    bcf_get_success_check = bcf_get_genotypes(header, rec, &gt_data, &num_gt_data);

    int32_t *depth_data = NULL;
    int32_t num_depth_data = 0;
    bcf_get_success_check = bcf_get_format_int32(header, rec, "DP", &depth_data, &num_depth_data);

    int ref_base_num = bcf_acgt2int(*rec->d.allele[0]);

    int i;
    for (i = 0; i < nsamples; i++)
    {
        int all1 = bcf_gt_allele(gt_data[2*i + 0]);
        int all2 = bcf_gt_allele(gt_data[2*i + 1]);
        bool is_pass = (strcmp(filter_data[i], "PASS") == 0) || (strcmp(filter_data[i], "No_var") == 0);

        if (is_pass && !args->use_pass)
            continue;

        if (!is_pass && !args->use_fail)
            continue;

        if (depth_data[i] >= MIN_COVERAGE && depth_data[i] <= MAX_COVERAGE) //errors are a big negative number, so skip
        {
            samp_buckets[i]->total_coverage += depth_data[i];
            samp_buckets[i]->genotypes_with_depth++;
        } 

        if (all1 == 0 && all2 == 0)
        {
            samp_buckets[i]->ref++;
        }
        else if (all1 != all2 && all1 >= 0 && all2 >= 0)
        {
            samp_buckets[i]->het++;
            
            if (is_pass)
                samp_buckets[i]->passing_variants++;           
 
            // stored as 0, 1, 2, 3 for A, C, G, T, respectively, so we can do this small madness
            if (!args->is_indel_file)
            {
                int all2_base_num = bcf_acgt2int(*rec->d.allele[all2]);
                if (abs(ref_base_num - all2_base_num) == 2)                  
                    samp_buckets[i]->transitions++;
                else
                    samp_buckets[i]->transversions++;
            }
        }
        else if (all1 == all2 && all1 >= 0 && all2 >= 0)
        {
            samp_buckets[i]->var++;
    
            if (is_pass)
                samp_buckets[i]->passing_variants++;

            if (!args->is_indel_file)
            {        
                int all1_base_num = bcf_acgt2int(*rec->d.allele[all1]);
                if (abs(ref_base_num - all1_base_num) == 2)                  
                    samp_buckets[i]->transitions++;
                else
                    samp_buckets[i]->transversions++;
                
                int all2_base_num = bcf_acgt2int(*rec->d.allele[all2]);
                if (abs(ref_base_num - all2_base_num) == 2)                  
                    samp_buckets[i]->transitions++;
                else
                    samp_buckets[i]->transversions++;
            } 
        }
        else
        {
            samp_buckets[i]->missing++;
        }
    }

    free(gt_data);
    free(depth_data);
    free(filter_data[0]);
    free(filter_data);

    return NULL;
}


/*
 *     Clean up.
 *     */
void destroy(void)
{
    printf("sample,variant_count,passing_variant_count,ti_tv_ratio,homref,hetvar,homvar,missing,het_hom_ratio,"
           "missing_rate,average_coverage,coverage_numerator,coverage_denominator\n");

    int i;
    for (i = 0; i < nsamples; i++)
    {
        printf("%s,", header->samples[i]);
        printf("%d,", samp_buckets[i]->het + samp_buckets[i]->var);
        printf("%d,", samp_buckets[i]->passing_variants);
        printf("%lf,", samp_buckets[i]->transitions / (double) samp_buckets[i]->transversions);
        printf("%d,", samp_buckets[i]->ref);
        printf("%d,", samp_buckets[i]->het);
        printf("%d,", samp_buckets[i]->var);
        printf("%d,", samp_buckets[i]->missing);
        printf("%lf,", samp_buckets[i]->het / (double) samp_buckets[i]->var);
        printf("%lf,", samp_buckets[i]->missing / (double) num_sites);
        printf("%lf,", samp_buckets[i]->total_coverage / (double) samp_buckets[i]->genotypes_with_depth);
        printf("%ld,", samp_buckets[i]->total_coverage);
        printf("%d,", samp_buckets[i]->genotypes_with_depth);
        printf("\n");
    }

    for (i = 0; i < nsamples; i++)
        free(samp_buckets[i]);
    free(samp_buckets);

    free(args);
}

