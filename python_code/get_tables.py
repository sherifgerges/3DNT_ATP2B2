import pandas as pd
import glob

gene_symbol_to_ensg = {
    'SETD1A': 'ENSG00000099381',
    'CUL1': 'ENSG00000055130',
    'XPO7': 'ENSG00000130227',
    'TRIO': 'ENSG00000038382',
    'CACNA1G': 'ENSG00000006283', 
    'SP4': 'ENSG00000105866',
    'GRIA3': 'ENSG00000125675',
    'GRIN2A': 'ENSG00000183454',
    'HERC1': 'ENSG00000103657',
    'RB1CC1': 'ENSG00000023287',
    'ATP2B2': 'ENSG00000157087',
    }

schema10 = [
    'SETD1A',
    'CUL1',
    'XPO7',
    'TRIO',
    'CACNA1G',
    'SP4',
    'GRIA3',
    'GRIN2A',
    'HERC1',
    'RB1CC1',
]

aa_mapping = {
    'Ala': 'A', 'Cys': 'C', 'Asp': 'D', 'Glu': 'E', 'Phe': 'F', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Lys': 'K', 'Leu': 'L', 'Met': 'M', 'Asn': 'N', 'Pro': 'P', 'Gln': 'Q', 'Arg': 'R', 'Ser': 'S',
    'Thr': 'T', 'Val': 'V', 'Trp': 'W', 'Tyr': 'Y', 'del': 'l', 'dup': 'p',
}

def standardize_mutation(hgvs):
    parts = hgvs.split('.')
    
    mutation = parts[1]  # Get the mutation part e.g., Ala4854Thr
    ref_aa = aa_mapping[mutation[:3]]  # Convert Ala -> A
    pos = ''.join([char for char in mutation if char.isdigit()])  # Extract the position e.g., 4854
    alt_aa = aa_mapping[mutation[-3:]]  # Convert Thr -> T
    return f"{ref_aa}{pos}{alt_aa}"  # Return in format A4854T

def get_schema_results(gene_symbol_list):
    # has to be in gene_symbol_to_ensg
    to_concat = []
    for gene_name in gene_symbol_list:
        ensg = gene_symbol_to_ensg[gene_name]
        schema_file = glob.glob(f'data/meta_{ensg}_variants_*.csv')[0]
        df_gene = pd.read_csv(schema_file)
        df_gene['gene_symbol'] = gene_name
        to_concat.append(df_gene)
    df_schema = pd.concat(to_concat)

    df_schema['locus'] = [':'.join(x.split('-')[0:2]) for x in df_schema['Variant ID']]
    df_schema = df_schema[df_schema['AC Case'] + df_schema['AC Control'] < 10]
    df_schema = df_schema[df_schema['HGVSp/c'].str.match('^p\.[A-Za-z]{3}\d+[A-Za-z]{3}$')]
    df_schema['Mutation'] = df_schema['HGVSp/c'].apply(standardize_mutation)
    df_schema['aa_pos'] = [int(x[1:-1]) for x in df_schema['Mutation']]
    return df_schema

def get_vsms(vsm_files):
    df_vep = pd.concat([pd.read_csv(f, sep='\t') for f in vsm_files], axis=0)   
    df_vep['Variant ID'] = list(['-'.join([a,b,c,d]) 
        for (a,b,c,d) in zip(
            [x.split(':')[0] for x in df_vep.locus],
            [x.split(':')[1] for x in df_vep.locus],
            [x.split('"')[1] for x in df_vep.alleles],
            [x.split('"')[3] for x in df_vep.alleles],
        )
    ])
    df_vep = df_vep.drop(['locus', 'alleles'], axis=1)
    return df_vep

def explode_table(dfh):
    df_case = dfh.loc[dfh.index.repeat(dfh['AC Case'])].reset_index(drop=True)
    df_case['is_case'] = 1
    df_control = dfh.loc[dfh.index.repeat(dfh['AC Control'])].reset_index(drop=True)
    df_control['is_case'] = 0
    dfe = pd.concat([df_case, df_control], axis=0)
    return dfe

def get_schema_and_vsms():
    df_schema = get_schema_results(schema10)
    df_vep = get_vsms(['data/vep_scores_schema10.tsv'])
    dfh = pd.merge(
        df_schema,
        df_vep,
        how='left',
        on=['gene_symbol','Variant ID'],
    )
    df = explode_table(dfh)
    return df

def get_schema_asd_from_mark():
    dfcase = pd.read_csv('data/atp2b2_case_variants_schema_asd.tsv')
    dfcase['is_case'] = 1
    dfcontrol = pd.read_csv('data/atp2b2_control_variants_schema_asd.tsv')
    dfcontrol['is_case'] = 0
    df = pd.concat([dfcase, dfcontrol])
    df['Mutation'] = df['HGVSp/c'].apply(standardize_mutation)
    df['aa_pos'] = [int(x[1:-1]) for x in df.Mutation]
    df_stability = pd.read_csv('data/ATP2B2_ThermoMPNN.csv')
    df_stability = df_stability.rename(columns={'ddG (kcal/mol)':'ddG'})
    df = df.merge(df_stability[['Mutation', 'ddG']])
    return(df)