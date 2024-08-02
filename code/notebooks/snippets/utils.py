import plotly.express as px
from sklearn.decomposition import PCA
from umap import UMAP
import pandas as pd

import requests
from time import sleep 
import json


# Oligo and LCM colors:
clrs = px.colors.qualitative.G10
syncom_colors = {'YL32': '#149AB3',
 'KB18': '#616161',
 'I48': '#C26215',
 'YL27': '#F59C46',
 'YL45': '#E30C4B',
 'I46': '#2C5D52',
 'YL58': '#163A1A',
 'YL31': '#099334',
 'YL2': '#282E68',
 'I49': '#3DB077',
 'YL44': '#AC7FB6',
 'KB1': '#42AB34', 
 'Salmonella': clrs[0],
 'ASF519': clrs[1],
 'Clostridium indolis Y18184 (?)': clrs[8],
 'unclassified Lachnospiraceae':clrs[4],
 'Turicibacter': clrs[2],
 'Staphylococcus': clrs[6],
 'contig_21': clrs[4],
 'test':'white',
}



# Genome Map

genome_map = {'CP015399.2': 'YL32',
                  'CP015400.2': 'KB18',
                  'CP015401.2': 'I48',
                  'CP015402.2': 'YL27',
                  'CP015403.2': 'YL45',
                  'CP015404.2': 'I46',
                  'CP015405.2': 'YL58',
                  'CP015406.2': 'YL31',
                  'CP015407.2': 'YL2',
                  'CP015408.2': 'I49',
                  'CP015409.2': 'YL44',
                  'CP015410.2': 'KB1',
                  'GCF_000364265': 'ASF519',
                  'FQ312003.1': 'SL1344',
                  'FQ312003.1;FQ312003.1': 'SL1344',
                  'HE654725.1': 'SL1344',
                  'HE654726.1': 'SL1344',
                  'HE654724.1': 'SL1344',
                  'contig_15': 'contig_15',
                  'contig_21': 'contig_21',
                  'contig_26': 'contig_26',
                  'contig_46': 'contig_46',
                  'AQFU02000001.1': 'ASF 502',
                  'AQFU02000002.1': 'ASF 502',
                  'AQFU02000003.1': 'ASF 502',
                  'CP097573.1': 'ASF500',
                  'NZ_CP097810.1': 'ASF356',
                  'NZ_AQFR02000001.1': 'ASF360',
                  'NZ_AQFR02000002.1': 'ASF360',
                  'NZ_AQFR02000003.1': 'ASF360',
                  'NZ_CP097561.1': 'ASF361',
                  'NZ_AQFT02000001.1': 'ASF492',
                  'NZ_AQFT02000002.1': 'ASF492',
                  'NZ_AQFT02000003.1': 'ASF492',
                  'NZ_AQFT02000004.1': 'ASF492',
                  'NZ_AQFT02000005.1': 'ASF492',
                  'NZ_AQFT02000006.1': 'ASF492',
                  'NZ_AQFV02000001.1': 'ASF519',
                  'NZ_AQFV02000002.1': 'ASF519',
                  'NZ_AQFV02000003.1': 'ASF519',
                  'NZ_AQFV02000004.1': 'ASF519',
                  'NZ_AQFV02000005.1': 'ASF519',
                  'NZ_AQFV02000006.1': 'ASF519',
                  'NZ_CP097562.1': 'ASF457'
                  }

# mOTUs (v3) map

motus_to_strains = {'Bacteroides caecimuris [ref_mOTU_v31_03476]': 'I48', 
                    'Muribaculum intestinale [ref_mOTU_v31_10099]': 'YL27',
                    'Blautia sp. YL58 [ref_mOTU_v31_02153]': 'YL58',
                    'Akkermansia muciniphila [ref_mOTU_v31_03591]':'YL44',
                    'Hungateiclostridiaceae bacterium KB18 [ref_mOTU_v31_10098]':'KB18',
                    'Flavonifractor plautii [ref_mOTU_v31_05238]': 'YL31',
                    '[Clostridium] clostridioforme/bolteae [ref_mOTU_v31_03442]': 'YL32',
                    'Enterococcus faecalis [ref_mOTU_v31_00318]': 'KB1',
                    'Salmonella enterica [ref_mOTU_v31_00099]': 'SL1344',
                    'Parabacteroides goldsteinii [ref_mOTU_v31_01679]': 'ASF519',
                    'Lactobacillus reuteri [ref_mOTU_v31_04085]': 'I49',
                    }



# NCBI map

ncbi_taxid_map= {'I48': '1796613', 'YL32': '1834196', 'YL58':'1796616', 'YL27': '1796646'}

# PCA functions 

def find_pcs(df, num_pcs=2, num_genes=500, choose_by='variance'):
    """
    :param numPCs:
    :param numGenes:
    :return:
    """
    if num_genes:
        # calculate var for each, pick numGenes top var across samples -> df
        if choose_by == 'variance':
            genes = df.var(axis=1).sort_values(ascending=False).head(num_genes).index
            fdf = df.loc[genes].T
        else:
            pass
            # todo implement log2fc selection
    else:
        fdf = df.T
    pca = PCA(n_components=num_pcs)
    principal_components = pca.fit_transform(fdf)
    pcs = [f'PC{i}' for i in range(1, num_pcs + 1)]
    pc_df = (pd.DataFrame(data=principal_components, columns=pcs).set_index(fdf.index))
    pc_var = {pcs[i]: round(pca.explained_variance_ratio_[i] * 100, 2) for i in range(0, num_pcs)}
    #pc_df = pc_df.merge(_self.sd, left_index=True, right_index=True)
    return pc_df, pc_var

def get_strain_pca(df, strain, sd, color_by='Treatment', symbol_by='saturation'):
    strain_df = df[df.genome == strain]
    strain_df = strain_df[['gene_id', 'sample_id', 'num_reads']].pivot(index='gene_id', columns='sample_id')
    strain_df.columns = [f[1] for f in strain_df.columns]
    strain_df = np.log2(strain_df/strain_df.sum()*1000000 + 0.5)
    pc_df, pc_var = find_pcs(strain_df, 2, 1000)
    pc_df = pc_df.reset_index().rename({'index':'sample_id'}, axis=1)
    pc_df = pc_df.merge(sd[sd.genome == strain], on='sample_id', how='left')
    fig = px.scatter(pc_df, x='PC1', y='PC2', color=color_by, symbol=symbol_by,
                    template='plotly_white', width=800, height=800, 
                     hover_data=[c for c in pc_df.columns if 'PC' not in c], title=strain)
    fig.update_traces(marker=dict(size=12,
                              line=dict(width=2,
                                        color='DarkSlateGrey')),
                  selector=dict(mode='markers'), )
    return fig, pc_df



def find_umap_df(df):
    umap_2d = UMAP(n_components=2, init='random', random_state=0)
    proj_2d = umap_2d.fit_transform(df.T)
    proj_2d = (pd.DataFrame(data=proj_2d, columns=["UMAP1", "UMAP2"]).set_index(df.T.index))
    return proj_2d


# STRING functions

def get_term_df(func_analysis, description):
    df = func_analysis[func_analysis.description.str.contains(description)]
    fdf = pd.concat([df[['inputGenes', 'description']].explode('inputGenes'), df[['preferredNames']].explode('preferredNames')], axis=1).drop_duplicates()
    fdf = fdf.rename(columns={'inputGenes': 'locus_tag'})
    return fdf

def process_strain(res, strain, taxid_map):
    sres = res[res.genome == strain].copy()
    sres['locus_tag'] = sres['ID'].str.replace("gene-", '')
    taxid = taxid_map[strain]
    sup = sres.query("log2FoldChange > 0.6 & padj < 0.05")
    sdown = sres.query("log2FoldChange < -0.6 & padj < 0.05")
    link_up = link_to_string(sup.locus_tag.values, taxid)
    link_down = link_to_string(sdown.locus_tag.values, taxid)
    func_up = pd.DataFrame(string_function(sup.locus_tag.values, taxid))
    func_down = pd.DataFrame(string_function(sdown.locus_tag.values, taxid))
    return sres, link_up, link_down, func_up, func_down

def highlight_term(sres, func_up, func_down, up_terms, down_terms):
    ups = pd.concat([get_term_df(func_up, term) for term in up_terms])
    downs = pd.concat(get_term_df(func_down, term) for term in down_terms)
    highlights = pd.concat([ups, downs])
    sres = sres.merge(highlights, on='locus_tag', how='left')
    sres['description'] = sres['description'].fillna('Other')
    return sres



def link_to_string(gene_names, species):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = 'tsv-no-header'
    method = 'get_link'
    request_url = "/".join([string_api_url, output_format, method])
    if len(gene_names) < 550:
        params = {
            "identifiers": "\r".join(gene_names),  # your protein
            "species": species,  # species NCBI identifier
            "network_flavor": "confidence",  # show confidence links
            "caller_identity": "explodata"  # your app name
        }
        network = requests.post(request_url, data=params)
        network_url = network.text.strip()
        sleep(0.5)
    else:
        return 'Oh No'
    return network_url


def string_function(gene_names, species):
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "json"
    method = "enrichment"


    ##
    ## Construct the request
    ##

    request_url = "/".join([string_api_url, output_format, method])

    ##
    ## Set parameters
    ##


    params = {

        "identifiers" : "\r".join(gene_names),  # your protein
        "species" : species, # species NCBI identifier 
        "caller_identity" : "test_api" # your app name

    }

    ##
    ## Call STRING
    ##

    response = requests.post(request_url, data=params)

    ##
    ## Read and parse the results
    ##

    data = json.loads(response.text)
    return data

def func_graph(df, c='Blues_r'):
    df = df[(df.term.str.startswith('GO')) | (df.term.str.startswith('KW') |(df.category == 'KEGG'))].copy()
    df['description'] = df['description'] +" | " + df['term']
    fig = px.bar(df.sort_values('number_of_genes'), y='description', x='number_of_genes', color='fdr', orientation='h',
    width =900, height=500, template= 'plotly_white', color_continuous_scale=c, 
    labels={'number_of_genes':'Number of genes', 'description': ''})
    return fig