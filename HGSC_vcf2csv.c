#include <stdio.h>
#include <stdlib.h>
//#include "../htslib-1.6/htslib/vcf.h"
#include <htslib/vcf.h>

int nsamples;


bcf_hdr_t *header;
/*
 *     This short description is used to generate the output of `bcftools plugin -l`.
 *     */
const char *about(void)
{
    return "";
}

/*
 *     Called once at startup, allows to initialize local variables.
 *         Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
 *         */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    header = in;
   
    if (argc != 1)
    {
        printf("Prints out a CSV of genotypes. Use at your peril!");
        exit(1);
    }
  
    //printf("var_id,var_type,");
 
    int i;
    for (i = 0; i < nsamples; i++)
    {
        printf("%s", header->samples[i]);

        if (i != nsamples - 1)
            printf(",");
        else
            printf("\n");

    }
 
    return 1;
    //return 0;
}


/*
 *     Called for each VCF record. Return rec to output the line or NULL
 *         to suppress output.
 *         */
bcf1_t *process(bcf1_t *rec)
{
    int bcf_get_success_check;

    char* var_id = rec->d.id;
    char* sv_type = rec->d.allele[1];

    int32_t *gt_data = NULL;
    int32_t num_gt_data = 0;
    bcf_get_success_check = bcf_get_genotypes(header, rec, &gt_data, &num_gt_data);

    char **filter_data = NULL;
    int32_t num_filter_data = 0;
    bcf_get_success_check = bcf_get_format_string(header, rec, "FT", &filter_data, &num_filter_data);

    //printf("%s,%s,", var_id, sv_type);

    int i;
    for (i = 0; i < nsamples; i++)
    {
        int all1 = bcf_gt_allele(gt_data[2*i + 0]);
        int all2 = bcf_gt_allele(gt_data[2*i + 1]);

        if ((strcmp(filter_data[i], "PASS") != 0) && (strcmp(filter_data[i], "No_var") != 0))
        {
            printf("-10");
        }
        else if (all1 < 0 || all2 < 0)
        {
            printf("-10");
        }
        else if (all1 == 0 && all2 == 0)
        {
            printf("0");
        }
        else if (all1 != all2)
        {
            printf("0.5");
        }
        else if (all1 == all2)
        {
            printf("1");
        }
        else  // Something is terribly wrong, probably?
        {
            printf("-100");
        }

        if (i != nsamples - 1)
            printf(",");
        else
            printf("\n");
    }

    free(filter_data[0]);
    free(filter_data);
    free(gt_data);

    return NULL;
}


/*
 *     Clean up.
 *     */
void destroy(void)
{
}

