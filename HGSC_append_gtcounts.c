#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>

//len(PASS, FAIL, NVAR, NDAT, NFLT) = 5
#define NUM_FILTER_STRINGS 5
//similarly for genotypes
#define NUM_GT_STRINGS 5

#define INDEX_PASS 0
#define INDEX_FAIL 1
#define INDEX_NVAR 2
#define INDEX_NDAT 3
#define INDEX_NFLT 4

#define INDEX_MISS 0
#define INDEX_0 1
#define INDEX_1 2
#define INDEX_2 3
#define INDEX_3 4

int nsamples;
bcf_hdr_t *header;
char *filter_strings[] = {"PASS", "FAIL", "NVAR", "NDAT", "NFLT"};
char *genotype_strings[] = {".", "0", "1", "2", "3"};


/*
 *     This short description is used to generate the output of `bcftools plugin -l`.
 *     */
const char *about(void)
{
    return "Add INFO subfields with counts of various FT and GT combinations.\n"
           "Be warned that this assumes SNPs! This plugin won't work if there are more than 4 alleles in one row.\n";
}

/*
 *     Called once at startup, allows to initialize local variables.
 *         Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
 *         */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    header = out;
  
    int ret;
    int filt_counter;
    int gt1_counter;
    int gt2_counter;
    
    for (filt_counter = 0; filt_counter < NUM_FILTER_STRINGS; filt_counter++)
    {
        for (gt1_counter = 0; gt1_counter < NUM_GT_STRINGS; gt1_counter++)
        {
            for (gt2_counter = 0; gt2_counter < NUM_GT_STRINGS; gt2_counter++)
            {
                char *header_string;
                ret = asprintf(&header_string, "##INFO=<ID=%s_%s%s_%d,Number=.,Type=Integer,Description=\"Count of "
                                               "samples w/filter %s and genotype %s/%s in %d-sample set\">",
                               filter_strings[filt_counter], genotype_strings[gt1_counter], 
                               genotype_strings[gt2_counter], nsamples,
                               filter_strings[filt_counter], genotype_strings[gt1_counter], 
                               genotype_strings[gt2_counter], nsamples);
                if (ret < 0)
                {
                    fprintf(stderr, "Error updating header.\n");
                    return -1;
                }
 
                ret = bcf_hdr_append(header, header_string);
                if (ret != 0)
                {
                    fprintf(stderr, "Error updating header.\n");
                    return -1;
                }
                
                free(header_string);
            }
        }
    }

    return 0;
}


/*
 *     Called for each VCF record. Return rec to output the line or NULL
 *         to suppress output.
 *         */
bcf1_t *process(bcf1_t *rec)
{
    int num_returned;

    char **filter_data = NULL;
    int32_t num_filter_data = 0;
    num_returned = bcf_get_format_string(header, rec, "FT", &filter_data, &num_filter_data);

    int32_t *gt_data = NULL;
    int32_t num_gt_data = 0;
    num_returned = bcf_get_genotypes(header, rec, &gt_data, &num_gt_data);

    //char *filter_string;  //one of PASS, FAIL, NVAR, NDAT
    //char* gt_string;    //one of .., .0, .1, 00, 01, or 11 for ./., ./0, ./1, 0/0, etc 

    int buckets[NUM_FILTER_STRINGS][NUM_GT_STRINGS][NUM_GT_STRINGS];
    int i, j, k;
    for (i = 0; i < NUM_FILTER_STRINGS; i++)
        for (j = 0; j < NUM_GT_STRINGS; j++)
            for (k = 0; k < NUM_GT_STRINGS; k++)
                buckets[i][j][k] = 0;

    int filt_index;
    int gt1_index;
    int gt2_index;

    for (i = 0; i < nsamples; i++)
    {
        int all1 = bcf_gt_allele(gt_data[2*i + 0]);
        int all2 = bcf_gt_allele(gt_data[2*i + 1]);

        //CHECK FILTER
        if(strcmp(filter_data[i], "PASS") == 0)
            filt_index = INDEX_PASS;
        else if(strcmp(filter_data[i], "No_var") == 0)
            filt_index = INDEX_NVAR;
        else if(strcmp(filter_data[i], "No_data") == 0)
            filt_index = INDEX_NDAT;
        else if(strcmp(filter_data[i], ".") == 0)
            filt_index = INDEX_NFLT;
        else
            filt_index = INDEX_FAIL;

        //CHECK GT
        if (all1 < 0)
            gt1_index = INDEX_MISS;
        if (all1 == 0)
            gt1_index = INDEX_0;
        if (all1 == 1)
            gt1_index = INDEX_1;
        if (all1 == 2)
            gt1_index = INDEX_2;
        if (all1 == 3)
            gt1_index = INDEX_3;
        
        if (all2 < 0)
            gt2_index = INDEX_MISS;
        if (all2 == 0)
            gt2_index = INDEX_0;
        if (all2 == 1)
            gt2_index = INDEX_1;
        if (all2 == 2)
            gt2_index = INDEX_2;
        if (all2 == 3)
            gt2_index = INDEX_3;

        buckets[filt_index][gt1_index][gt2_index]++;
    }

    for (i = 0; i < NUM_FILTER_STRINGS; i++)
        for (j = 0; j < NUM_GT_STRINGS; j++)
            for (k = 0; k < NUM_GT_STRINGS; k++)
            {
                if (buckets[i][j][k] == 0)
                    continue;        
    
                char *id_str = NULL;
                int ret = asprintf(&id_str, "%s_%s%s_%d", 
                                   filter_strings[i], genotype_strings[j], genotype_strings[k], nsamples);

                if (ret < 0 || bcf_update_info_int32(header, rec, id_str, &buckets[i][j][k], 1) < 0) 
                {
                    fprintf(stderr, "Error adding %s:-/\n", id_str); 
                    exit(1); 
                }

                free(id_str);
            } 

    free(filter_data[0]);
    free(filter_data);
    free(gt_data);

    return rec;
}


/*
 *     Clean up.
 *     */
void destroy(void)
{
}

