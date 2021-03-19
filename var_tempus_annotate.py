"""
Tempus coding excerise
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

import requests
import vcf

# function containing all the VCF parsing steps

#List of consquences
#https://gemini.readthedocs.io/en/latest/content/database_schema.html,
# used here as sets for faster computation of whether an element is present

SEV_HIGH = ("exon_loss_variant", "frameshift_variant",
            "splice_acceptor_variant", "splice_donor_variant", "start_lost",
            "stop_gained", "stop_lost", "initiator_codon_variant",
            "rare_amino_acid_variant", "chromosomal_deletion")
SEV_MED = ("missense_variant", "inframe_insertion", "inframe_deletion",
           "coding_sequence_variant", "disruptive_inframe_deletion",
           "disruptive_inframe_insertion",
           "5_prime_UTR_truncation + exon_loss_variant",
           "3_prime_UTR_truncation + exon_loss_variant",
           "splice_region_variant", "mature_miRNA_variant",
           "regulatory_region_variant", "TF_binding_site_variant",
           "regulatory_region_ablation", "regulatory_region_amplification",
           "TFBS_ablation", "TFBS_amplification")
SEV_LOW = ("stop_retained_variant", "synonymous_variant",
           "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant",
           "coding_sequence_variant", "upstream_gene_variant",
           "downstream_gene_variant", "intergenic_variant",
           "intragenic_variant", "gene_variant", "transcript_variant",
           "exon_variant", "5_prime_UTR_premature_start_codon_gain_variant",
           "start_retained_variant", "conserved_intron_variant",
           "nc_transcript_variant", "NMD_transcript_variant",
           "incomplete_terminal_codon_variant",
           "non_coding_transcript_exon_variant", "transcript_ablation",
           "transcript_amplification", "feature_elongation",
           "feature_truncation", "unsure")


def bulk_api_anntation(csv_annotate, bulk_string, find_pattern_array):
    """
      This function  calls the hms harvard api and annotates the dict with it's effects in bulk
      input: 
	csv_annotate(dict): a dictionary of fields for the variant file of the vcf
        bulk_string(string):  a string to push through the api to get the bulk callback
        find_pattern_array(list): Important for keeping order since the api returns an unorderd bulk
      output:
	csv_annotate(dict): a dictionary of fields added with the effects of the variant added
    """
    url = 'http://exac.hms.harvard.edu/rest/bulk/variant'
    bulk_lines = bulk_string
    bulk_lines += "]'"
    r = requests.post(url, eval(bulk_lines))
    exac_annotations = json.loads(r.text)
    variant_consequences = ""

    #go through the orded list so the data matches the dict lists in csv_annotate
    for variant in find_pattern_array:

        # Get the allele frequencies from EXaC
        if 'allele_freq' in exac_annotations[variant]['variant']:
            allele_freq = exac_annotations[variant]['variant']['allele_freq']
        else:
            allele_freq = 'NA'

        csv_annotate['AlleFreq'].append(allele_freq)

        # Get the variant consequence from EXaC
        api_consquence = exac_annotations[variant]['consequence']
        if api_consquence is None or api_consquence == {}:
            csv_annotate['Consequence'].append('Unknown')
            csv_annotate['Genes'].append('NA')
        else:
            # Loop through consequences for the variants that have more than one consequence
            if len(api_consquence) > 1:
                #if a high severity is returned, we keept that and move on, if not we continue to search until we run out of consquences or find one.
                level = 0
                for cons in api_consquence.keys():
                    if cons in SEV_HIGH:
                        variant_consequences = cons
                        break
                    #will only update if the medium severity is higher than the current consquence
                    elif cons in SEV_MED and level < 1:
                        variant_consequences = cons
                        level = 1
                    elif level < 1:
                        variant_consequences = cons

            else:
                # Get the consequence for variants that only have one, also get the genes
                variant_consequences = list(api_consquence.keys())[0]
            csv_annotate['Genes'].append(
                list(api_consquence[variant_consequences].keys())[0])
            csv_annotate['Consequence'].append(variant_consequences)

    return csv_annotate


#function for parsing
def vcf_annotate(var_file, out_file):
    """
      this function parses the var_file and adds annotations
      input: 
       var_file(string): name of file to parse
       out_file(string) : name of file to save to 
 
      output:
         saved_file(string): a file of annoations from vcf file
    """
    # open and parse vcf file
    vcf_reader = vcf.Reader(open(var_file, 'r'))
    num_lines = 0
    #count the amount of lines in the VCF so the program knows when to stop
    for line in vcf_reader:
        num_lines += 1

    # parse through vcf file
    vcf_reader = vcf.Reader(open(var_file, 'r'))

    csv_annotate = defaultdict(list)
    find_pattern_array = []
    for counter, line in enumerate(vcf_reader):

        # We will send out the bulk api every 350 entries ( or less if its at the end), so we restart it here
        if counter % 351 == 0:
            bulk_lines = "'[\\"
        else:
            bulk_lines += ",\\"
        # collect chrom, position, ref, alt for the line , and the orded array for the csv/tsv annotations.
        csv_annotate['Chr'].append(str(line.CHROM))
        csv_annotate['Position'].append(str(line.POS))
        csv_annotate['Ref'].append(str(line.REF))
        csv_annotate['Alt'].append(str(line.ALT[0]))

        plug_in = csv_annotate['Chr'][-1] + "-" + csv_annotate['Position'][-1] + "-"
        plug_in += csv_annotate['Ref'][-1] + "-" + csv_annotate['Alt'][-1]
        bulk_lines += f'"{plug_in}\\"'
        find_pattern_array.append(plug_in)
        # obtain the Annotation from the type
        try:
            csv_annotate['Annotation'].append(str(line.INFO['TYPE'][0]))
        except Exception as error:
            print(" no Annotation, exception: " + str(error))
            csv_annotate['Annotation'].append('NA')
        # calculate the depth of the read
        try:
            csv_annotate['DP'].append(str(line.INFO['DP']))
        except Exception as error:
            print(" no depth, exception: " + str(error))
            csv_annotate['DP'].append('NA')
        # calculate number of variant supporting reads

        try:
            return_box = []
            for i in line.INFO['AO']:
                return_box.append(i)
            csv_annotate['VSReads'].append(return_box)
        except Exception as error:
            print(" no variant supporting reads, exception: " + str(error))
            csv_annotate['VSReads'].append('NA')

        # calculate number of reference allele observations
        try:
            csv_annotate['RO'].append(str(line.INFO['RO']))
        except Exception as error:
            print(" no refrence allele, exception: " + str(error))
            csv_annotate['RO'].append('NA')

        # calculate total read depth at variant location seq depth
        try:
            return_box = []
            for i in csv_annotate['VSReads'][-1]:
                return_box.append(
                    str(float(csv_annotate['RO'][-1]) + float(i)))
            csv_annotate['SeqDepth'].append(return_box)
        except Exception as error:
            print(" no sequence total read, exception: " + str(error))
            csv_annotate['SeqDepth'].append('NA')

        # calculate percent of variant supporting reads, to 4 floats
        try:
            return_box = []
            for i in range(len(csv_annotate['SeqDepth'][-1])):

                return_box.append(
                    str(
                        round(
                            float(csv_annotate['VSReads'][-1][i]) / float(
                                csv_annotate['SeqDepth'][-1][i]) * 100, 4)))
            csv_annotate['PVSReads'].append(return_box)
        except Exception as error:
            print(" no percent reads, exception: " + str(error))
            csv_annotate['PVSReads'].append('NA')

        # call API every 350 lines
        if counter % 351 == 350:
            csv_annotate = bulk_api_anntation(csv_annotate, bulk_lines,
                                              find_pattern_array)
            find_pattern_array = []

        if counter > num_lines:
            csv_annotate = bulk_api_anntation(csv_annotate, bulk_lines,
                                              find_pattern_array)
            break

    #needed to hardcode this to make sure it keeps order on the csv/tsv

    fields = [
        'Chr', 'Position', 'Ref', 'Alt', 'Annotation', 'DP', 'VSReads', 'RO',
        'SeqDepth', 'PVSReads', 'Consequence', 'AlleFreq', 'Genes'
    ]
    delim = ','
    if out_file[out_file.index(".") + 1:] == 'tsv':
        delim = '\t'
    with open(out_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=delim)
        writer.writerow(fields)
        writer.writerows(zip(*[csv_annotate[key] for key in fields]))


if __name__ == '__main__':

    # Allow users to specify the location of VCF file
    parser = argparse.ArgumentParser()
    parser.add_argument("-v",
        'vcf', help='VCF file location, put the input file here')
    parser.add_argument("--out",
        'outputfile',
        help='output file name, put the name of the output file you would like'
    )

    #shows help if theres nothing to do
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    if not os.path.isfile(args.vcf):
        print('Unable to locate VCF file: {0}'.format(args.vcf))
        sys.exit(1)
    vcf_annotate(args.vcf, args.outputfile)
