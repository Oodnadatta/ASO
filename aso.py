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

def get_gene_phenotype_from_omim_file(omim_genemap_file_path):
    omim_dict = {}
    with open(omim_genemap_file_path) as omim_file:
        for row in csv.DictReader(omim_file, delimiter='\t'):
            gene_symbol_list = []
            for gene_symbol in row['Gene/Locus And Other Related Symbols'].split(','):
                gene_symbol_list.append(gene_symbol.strip())
            phenotype_list = []
            for phenotype in row['Phenotypes'].split(';'):
                phenotype_list.append(phenotype.strip())
                # TODO FIXME: parse phenotype_list:
                # 1) phenotype = phenotype.strip().split(',')[0]
                # 2) mim_number = phenotype.strip().split(',')[1]
                # 3) transmission_mode = phenotype.strip().split(',')[2]
            gene_phenotype = {
                'Gene/Locus And Other Related Symbols': gene_symbol_list,
                'Approved Gene Symbol': row['Approved Gene Symbol'],
                'Phenotypes': phenotype_list,
            }
            for gene in gene_symbol_list:
                omim_dict[gene] = gene_phenotype
    return omim_dict

def add_omim_info(omim_genemap_file_path, green_genes_dict):
    omim_dict = get_gene_phenotype_from_omim_file(omim_genemap_file_path)
    for gene_dict in green_genes_dict.values():
        if gene_dict['Gene Symbol'] in omim_dict:
            gene_dict['Approved Gene Symbol'] = omim_dict[gene_dict['Gene Symbol']]['Approved Gene Symbol']
            gene_dict['OMIM Phenotypes'] = omim_dict[gene_dict['Gene Symbol']]['Phenotypes']
            gene_dict['Transmission mode'] = 'not implemented'
        else:
            gene_dict['Approved Gene Symbol'] = 'missing'
            gene_dict['OMIM Phenotypes'] = 'missing'
            gene_dict['Transmission mode'] = 'missing'
    
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
    
                
def get_nm_mane_from_gnomad4_constraint_file(gnomad4_file_path):
    gnomad4_dict = {}
    with open(gnomad4_file_path) as gnomad4_file:
        for row in csv.DictReader(gnomad4_file, delimiter='\t'):
            if row['mane_select'] == 'true':
                if row['transcript'].startswith('NM_'):
                    gnomad4_dict[row['gene']] = row
    return gnomad4_dict
                    
def add_gnomad4_info(gnomad4_file_path, green_genes_dict):
    gnomad4_dict = get_nm_mane_from_gnomad4_constraint_file(gnomad4_file_path)
    for gene_dict in green_genes_dict.values():
        if gene_dict['Approved Gene Symbol'] in gnomad4_dict:
            gene_dict['lof.pLI'] = gnomad4_dict[gene_dict['Approved Gene Symbol']]['lof.pLI']
            gene_dict['lof.oe'] = gnomad4_dict[gene_dict['Approved Gene Symbol']]['lof.oe']
            gene_dict['lof.oe_ci.upper'] = gnomad4_dict[gene_dict['Approved Gene Symbol']]['lof.oe_ci.upper']
        else:
            gene_dict['lof.pLI'] = 'missing'
            gene_dict['lof.oe'] = 'missing'
            gene_dict['lof.oe_ci.upper'] = 'missing'

def get_constraint_from_gnomad2_constraint_file(gnomad2_file_path):
    gnomad2_dict = {}
    with open(gnomad2_file_path) as gnomad2_file:
        for row in csv.DictReader(gnomad2_file, delimiter='\t'):
            gnomad2_dict[row['gene']] = row
    return gnomad2_dict
                    
def add_gnomad2_info(gnomad2_file_path, green_genes_dict):
    gnomad2_dict = get_constraint_from_gnomad2_constraint_file(gnomad2_file_path)
    for gene_dict in green_genes_dict.values():
        if gene_dict['Gene Symbol'] in gnomad2_dict:
            gene_dict['pLI'] = gnomad2_dict[gene_dict['Gene Symbol']]['pLI']
            gene_dict['oe_lof'] = gnomad2_dict[gene_dict['Gene Symbol']]['oe_lof']
            gene_dict['oe_lof_upper'] = gnomad2_dict[gene_dict['Gene Symbol']]['oe_lof_upper']
        else:
            gene_dict['pLI'] = 'missing'
            gene_dict['oe_lof'] = 'missing'
            gene_dict['oe_lof_upper'] = 'missing'
            
def display_genes_dict(genes_dict):
    print('Gene_Symbol\tApproved_gene_symbol\tOMIM_Phenotype\tTransmission_mode\tPanelapp_Phenotype\tP/LP_missense_count\tP/LP_PSC_count\tlof.pLI.v4\tlof.oe.v4\tLOEUF.v4\tlof.pLI.v2\tlof.oe.v2\tLOEUF.v2')
    for gene_dict in genes_dict.values():
        print('\t'.join([
            gene_dict['Gene Symbol'],
            gene_dict['Approved Gene Symbol'],
            str(gene_dict['OMIM Phenotypes']),
            gene_dict['Transmission mode'],
            '"' + gene_dict['Phenotypes'] + '"',
            str(gene_dict['P/LP_missense_count']),
            str(gene_dict['P/LP_premature_stop_codon_count']),
            gene_dict['lof.pLI'],
            gene_dict['lof.oe'],
            gene_dict['lof.oe_ci.upper'],
            gene_dict['pLI'],
            gene_dict['oe_lof'],
            gene_dict['oe_lof_upper']
        ]))
    
if __name__ == "__main__":
    green_genes_dict = get_green_from_panelapp_file(sys.argv[1])
    add_omim_info(sys.argv[2], green_genes_dict)
    add_clinvar_info(sys.argv[3], green_genes_dict)
    add_gnomad4_info(sys.argv[4], green_genes_dict)
    add_gnomad2_info(sys.argv[5], green_genes_dict)
    display_genes_dict(green_genes_dict)
