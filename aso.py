#! /usr/bin/env python3

# Anne-Sophie Denommé-Pichon
# asdeno@hotmail.fr

import csv
import sys

def open_panelapp_file(file_path):
    with open(file_path) as panelapp_file:
        for row in csv.DictReader(panelapp_file, delimiter='\t'):
            if 'Expert Review Green' in row['Sources(; separated)'].split(';'): 
                print(row['Gene Symbol'] + '\t' + row['Phenotypes'])

def open_clinvar_file(file_path):
    with open(file_path) as clinvar_file:
        for line in clinvar_file:
            if not line.startswith('#'):
                info = line.split('\t')[7]
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
                print(row_dict)
#                CLNSIG
#                GENEINFO
#                MC
                break

            
if __name__ == "__main__":
    print("hello world")
    open_panelapp_file(sys.argv[1])
    open_clinvar_file(sys.argv[2])
    
