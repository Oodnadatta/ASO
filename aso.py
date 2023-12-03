#! /usr/bin/env python3

# Anne-Sophie Denommé-Pichon
# asdeno@hotmail.fr

import csv
import pprint
import sys

def get_green_from_panelapp_file(panelapp_file_path):
    panelapp_dict = {}
    with open(panelapp_file_path) as panelapp_file:
        for row in csv.DictReader(panelapp_file, delimiter='\t'):
            if 'Expert Review Green' in row['Sources(; separated)'].split(';'): 
                if row['Gene Symbol'] != '':
                    panelapp_dict[row['Gene Symbol']] = row
    return panelapp_dict
    # {'AARS': {'Entity Name': 'AARS', 'Entity type': 'gene', 'Gene Symbol': 'AARS', 'Sources(; separated)': 'Wessex and West Midlands GLH;NHS GMS;Victorian C\linical Genetics Services;Expert Review Green;Literature', 'Level4': 'Early onset or syndromic epilepsy', 'Level3': '', 'Level2': '', 'Model_Of_Inherita\nce': 'BIALLELIC, autosomal or pseudoautosomal', 'Phenotypes': 'Developmental and epileptic encephalopathy 29, OMIM:616339;Developmental and epileptic e\ncephalopathy, 29, MONDO:0014593', 'Omim': '601065', 'Orphanet': '', 'HPO': '', 'Publications': '25817015;28493438', 'Description': '', 'Flagged': '', '\GEL_Status': '3', 'UserRatings_Green_amber_red': '', 'version': '4.0', 'ready': '', 'Mode of pathogenicity': '', 'EnsemblId(GRch37)': 'ENSG00000090861',\ 'EnsemblId(GRch38)': 'ENSG00000090861', 'HGNC': 'HGNC:20', 'Position Chromosome': '', 'Position GRCh37 Start': '', 'Position GRCh37 End': '', 'Position\ GRCh38 Start': '', 'Position GRCh38 End': '', 'STR Repeated Sequence': '', 'STR Normal Repeats': '', 'STR Pathogenic Repeats': '', 'Region Haploinsuffi\ciency Score': '', 'Region Triplosensitivity Score': '', 'Region Required Overlap Percentage': '', 'Region Variant Type': '', 'Region Verbose Name': ''}

def enumerate_clinvar_data(clinvar_file_path):
    with open(clinvar_file_path) as clinvar_file:
        for line in clinvar_file:
            if not line.startswith('#'):
                info = line.rstrip().split('\t')[7]
                row_dict = {}
                for key_value in info.split(';'):
                    key, value = key_value.split('=')
                    if key == 'GENEINFO':
                        geneinfo = []
                        for gene in value.split('|'):
                            geneinfo.append(gene.split(':'))
                        value = geneinfo
                    elif key == 'MC':
                        mc = []
                        for molecular_consequence in value.split(','):
                            mc.append(molecular_consequence.split('|'))
                        value = mc
                    row_dict[key] = value
                yield row_dict

def add_clinvar_info(clinvar_file_path, green_genes_dict):
    for green_gene in green_genes_dict.values():
        green_gene['P/LP_missense_count'] = 0 # Add 'missense_count' in green_genes_dict
        green_gene['P/LP_premature_stop_codon_count'] = 0
    for clinvar_row in enumerate_clinvar_data(clinvar_file_path):
        if 'GENEINFO' in clinvar_row:
            for gene, omim in clinvar_row['GENEINFO']:
                if gene in green_genes_dict.keys():
                    if 'CLNSIG' in clinvar_row:
                        if clinvar_row['CLNSIG'] in ['Likely_pathogenic',
                                                     'Pathogenic',
                                                     'Pathogenic/Likely_pathogenic',
                                                     'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance',
                                                     'Pathogenic/Likely_risk_allele',
                                                     'Pathogenic|drug_response|other']:
                            variant_types = []
                            if 'MC' in clinvar_row:
                                is_missense = False
                                is_premature_stop_codon = False
                                for molecular_consequence in clinvar_row['MC']:
                                    # clinvar_row['MC'] : value from 'MC' = [['SO:0001624', '3_prime_UTR_variant']]
                                    variant_type = molecular_consequence[1]
                                    #variant_types.append(variant_type)
                                # print(clinvar_row['CLNSIG'], gene, variant_types) = 5 Uncertain_significance GNB1 ['3_prime_UTR_variant']
                                    if variant_type in ['inframe_deletion',
                                                        'inframe_indel',
                                                        'inframe_insertion',
                                                        'missense_variant']:
                                        is_missense = True
                                    elif variant_type in ['frameshift_variant',
                                                          'nonsense',
                                                          'splice_acceptor_variant',
                                                          'splice_donor_variant']:
                                        is_premature_stop_codon = True
                                    elif variant_type in ['3_prime_UTR_variant',
                                                          '5_prime_UTR_variant',
                                                          'genic_downstream_transcript_variant',
                                                          'genic_upstream_transcript_variant',
                                                          'initiator_codon_variant',
                                                          'intron_variant',
                                                          'no_sequence_alteration',
                                                          'non-coding_transcript_variant',
                                                          'stop_lost',
                                                          'synonymous_variant']:
                                        pass
                                    else:
                                        print(f'Unknown variant_type: {variant_type}', file=sys.stderr)
                                        sys.exit(1)
                                if is_missense:
                                    green_genes_dict[gene]['P/LP_missense_count'] += 1
                                if is_premature_stop_codon:
                                    green_genes_dict[gene]['P/LP_premature_stop_codon_count'] += 1
    
                
def get_nm_mane_from_gnomad_constraint_file(gnomad_file_path):
    gnomad_dict = {}
    with open(gnomad_file_path) as gnomad_file:
        for row in csv.DictReader(gnomad_file, delimiter='\t'):
            if row['mane_select'] == 'true':
                if row['transcript'].startswith('NM_'):
                    gnomad_dict[row['gene']] = row
    return gnomad_dict
                    
def add_gnomad_info(gnomad_file_path, green_genes_dict):
    gnomad_dict = get_nm_mane_from_gnomad_constraint_file(gnomad_file_path)
    for gene_dict in green_genes_dict.values():
        if gene_dict['Gene Symbol'] in gnomad_dict: #FIXME with HGNC
            gene_dict['lof.pLI'] = gnomad_dict[gene_dict['Gene Symbol']]['lof.pLI']
            gene_dict['lof.oe'] = gnomad_dict[gene_dict['Gene Symbol']]['lof.oe']
            gene_dict['lof.oe_ci.upper'] = gnomad_dict[gene_dict['Gene Symbol']]['lof.oe_ci.upper']
        else: #FIXME with HGNC
            gene_dict['lof.pLI'] = 'missing'
            gene_dict['lof.oe'] = 'missing'
            gene_dict['lof.oe_ci.upper'] = 'missing'
     
def display_genes_dict(genes_dict):
    print('Gene_Symbol\tPhenotype\tP/LP_missense_count\tP/LP_PSC_count\tlof.pLI\tlof.oe\tLOEUF')
    for gene_dict in genes_dict.values():
        print('\t'.join([
            gene_dict['Gene Symbol'],
            gene_dict['Phenotypes'],
            str(gene_dict['P/LP_missense_count']),
            str(gene_dict['P/LP_premature_stop_codon_count']),
            gene_dict['lof.pLI'],
            gene_dict['lof.oe'],
            gene_dict['lof.oe_ci.upper']
        ]))
    
if __name__ == "__main__":
    green_genes_dict = get_green_from_panelapp_file(sys.argv[1])
    add_clinvar_info(sys.argv[2], green_genes_dict)
    add_gnomad_info(sys.argv[3], green_genes_dict)
    display_genes_dict(green_genes_dict)

    
