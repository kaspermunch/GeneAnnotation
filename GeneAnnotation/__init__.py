

import sys
import pandas as pd

#########################################################################
# Entreaz api
#########################################################################

from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "kaspermunch@birc.au.dk"

def retrieve_annotation(id_list):
    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was:", e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)
    return annotations

def get_uid(name):
    handle = Entrez.esearch(db="gene", retmax=10, 
        term="Homo sapiens[ORGN] AND current only[Filter] AND {}[GENE]".format(name))
    record = Entrez.read(handle)
#    print(record)
    handle.close()
    assert len(record['IdList']), "name {} maps to more than one uid".format(name)
    return record['IdList'][0]



def get_coordinates(gene_names):

    uid_list = [get_uid(name) for name in gene_names]
    annot_data = retrieve_annotation(uid_list)
    #print(annot_data)

    coord = dict()
    for record in annot_data['DocumentSummarySet']['DocumentSummary']:
        # for d in (record['LocationHist']):
        #     print(d)
        coord_info = record['GenomicInfo'][0]
        coord[record['NomenclatureSymbol']] = \
            (coord_info['ChrLoc'], coord_info['ChrStart'], coord_info['ChrStop'])
    return coord

#########################################################################
# Mygene API
#########################################################################

import mygene

def get_gene_coords(gene_list):
    mg = mygene.MyGeneInfo()

    records = list()
    for gene in mg.querymany(gene_list, scopes='symbol,accession,entrezgene,ensemblgene', 
        species='human', 
        fields='symbol,genomic_pos_hg19,genomic_pos,type_of_gene,summary'):
        #fields='all'):

        if 'genomic_pos_hg19' not in gene:
            print('skipping', gene['query'])
            continue

        coord = gene['genomic_pos_hg19']
        if not isinstance(coord, list):
            coord = [coord]
        for c in coord:
            records.append((gene['query'], gene['symbol'], c['chr'], c['start'], c['end'], gene['type_of_gene'], gene['summary']))

    return pd.DataFrame.from_records(records, columns=['query', 'symbol', 'chrom', 'start', 'end', 'type_of_gene', 'summary'])



# with open('../Downloads/allVIPs_withfactors_file') as f:
# 	gene_list = f.read().split()
# 	df = get_gene_coords(gene_list)
# print(df)
# assert 0

print(get_gene_coords(['FOXP2', 'BRCA2']))

# #print(mg.getgenes(['FOXP2', 'ENSG00000123374']))
# print(mg.getgenes([1017, '1018','ENSG00000148795'], scopes='symbol,reporter,accession', species='human', as_dataframe=True))

# print(mg.getgenes(['FOXP2', 'BRCA2']))
# out = mg.querymany(['FOXP2', 'BRCA2', 'ENSG00000123374'], scopes='symbol,reporter,accession', fields='entrezgene', species='human')
#out = mg.querymany(['FOXP2', 'BRCA2', 'ENSG00000123374'],  fields='entrezgene', species='human')
# print(out)


#print(get_coordinates(['FOXP2', 'BRCA2']))
#Entrez.read(Entrez.esearch(db="gene", retmax=10, term="Homo sapiens[ORGN] AND FOXP2[GENE]"))


