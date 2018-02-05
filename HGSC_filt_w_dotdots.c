#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <htslib/vcf.h>

int nsamples, nsnps, nindels, nmnps, nothers, nsites;

int min_depth;

bcf_hdr_t *header;
/*
 *     This short description is used to generate the output of `bcftools plugin -l`.
 *     */
const char *about(void)
{
    return "First arg should be filename, second arg should be minimum depth.\n"
           "This will change the \"FT\" subfield of each sample's FORMAT field in the\n"
           "following ways\n:"
           "* If the \"FT\" subfield is not empty, do nothing\n"
           "* Else write \"PASS\" if the \"GT\" subfield has a variant allele\n"
           "* Else if the \"GT\" subfield is \"./0\", write \".\"\n"
           "* Else write \"No_var\" if the \"GT\" subfield is \"0/0\" and \"DP\" is >= $MIN_DEPTH\n"
           "* Else write \"No_data\" if the \"GT\" subfield is \"0/0\" and \"DP\" is < $MIN_DEPTH\n"
           "This step will also change the \"GT\" subfield to \"./.\" if the \"FT\"\n"
           "subfield was marked as \"No_data.\"\n";
}

/*
 *     Called once at startup, allows to initialize local variables.
 *         Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
 *         */
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    nsamples = bcf_hdr_nsamples(in);
    nsnps = nindels = nmnps = nothers = nsites = 0;
    header = in;
   
    if (argc < 2)
    {
        printf("%s", about());
        exit(1);
    }

    char *endptr;
    errno = 0;
    min_depth = (int) strtol(argv[1], &endptr, 10); 
    if (errno != 0 || argv[1] == endptr)
    {
        printf("Invalid input.\n");
        exit(1);
    }
 
    //return 1;
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

    int32_t *depth_data = NULL;
    int32_t num_depth_data = 0;
    num_returned = bcf_get_format_int32(header, rec, "DP", &depth_data, &num_depth_data);
 
    const char **new_filter_data = calloc(nsamples, sizeof(char*));

    int i;
    for (i = 0; i < nsamples; i++)
    {
        int all1 = bcf_gt_allele(gt_data[2*i + 0]);
        int all2 = bcf_gt_allele(gt_data[2*i + 1]);

        if(strcmp(filter_data[i], ".") != 0)
            new_filter_data[i] = filter_data[i];
        else if (all1 == 0 && all2 == 0) 
        {
            if (depth_data[i] < min_depth)
            {
                gt_data[2*i + 0] = bcf_gt_missing;
                gt_data[2*i + 1] = bcf_gt_missing;

                new_filter_data[i] = "No_data";
            }
            else
                new_filter_data[i] = "No_var"; 
        }
        else if (all1 < 0 && all2 >= 0)
            new_filter_data[i] = ".";
        else if (all1 < 0 && all2 < 0)
            new_filter_data[i] = "No_data";
        else
            new_filter_data[i] = "PASS";
    }
 
    bcf_update_format_string(header, rec, "FT", new_filter_data, nsamples);
    free(filter_data[0]);
    free(filter_data);

    bcf_update_genotypes(header, rec, gt_data, num_gt_data); 

    free(gt_data);
    free(depth_data);
    free(new_filter_data);

    return rec;
}


/*
 *     Clean up.
 *     */
void destroy(void)
{
}

