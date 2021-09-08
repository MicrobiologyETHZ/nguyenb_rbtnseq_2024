import requests
from requests.exceptions import HTTPError
import datetime as dt
import pandas as pd
import os
import re
import sys
import aiohttp
import asyncio

class Gene(object):

    def __init__(self, gene_id=""):
        self.gene_id = gene_id  # This could be either locus tag or gene name
        self.entry = ""
        self.name = ""
        self.definition = ""
        self.pathways = []
        self.modules = []
        self.aa_seq = ""
        self.nt_seq = ""
        self.pathway_string = ""
        self.module_string = ""

    def get_full_info(self):
        return "\t".join([self.entry, self.name, self.definition,
                          self.pathway_string]).rstrip('\t')

    def get_short_info(self):
        return "\t".join([self.entry, self.name, self.definition])

    def get_aa_seq(self):
        return ">{}|{}\n{}\n".format(self.gene_id, self.name, self.aa_seq)

    def get_nt_seq(self):
        return ">{}|{}\n{}\n".format(self.gene_id, self.name, self.nt_seq)

    def __str__(self):
        return "GENE:\n{}\t{}\t{}\nPATHWAYS:\n{}".format(self.entry, self.name,
                                                         self.definition,
                                                         "\n".join(self.pathways))


    def kegg_get(self):
        try:
            url = "http://rest.kegg.jp/get/{}".format(self.gene_id)

            entry = requests.get(url).text.strip()
            lines = entry.split("\n")
            for i, line in enumerate(lines):
                field = line.split()[0]
                if field == 'ENTRY':
                    self.entry = lines[i].split()[1]
                elif field == "NAME":
                    self.name = lines[i].split("NAME        ")[1]
                elif field == "DEFINITION":
                    p = lines[i].split("DEFINITION  ")[1]
                    self.definition = re.sub(r"^\(.*?\)", "", p).strip()
                elif field == "PATHWAY":
                    p = lines[i].split("PATHWAY     ")[1]
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.pathways = p.split("            ")
                    self.pathway_string = ", ".join(self.pathways).rstrip(", ")
                elif field == "MODULE":
                    p = lines[i].split("MODULE     ")[1]
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.modules = p.split("            ")
                    self.module_string = ", ".join(self.pathways).rstrip(", ")
                elif field == "AASEQ":
                    p = ''
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.aa_seq = p.replace("            ", '')
                elif field == "NTSEQ":
                    p = ''
                    count = 1
                    while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                        p += lines[i + count]
                        count += 1
                    self.nt_seq = p.replace("            ", '')
        except HTTPError:
            print("Bad gene identifier")


class GeneSet(object):

    """
    Given either a set of genes or a genome (ex. sey, eco) retrieve and store KEGG info
    Gene Set is consists of Gene objects
    kegg ids
    gene names

    """
    def __init__(self, gene_id_list="", genome="", out_dir="."):
        if genome:
            self.genome = genome
            if gene_id_list:
                self.kegg_ids = [f"{genome}:{gene}" for gene in gene_id_list]
            else:
                self.kegg_ids = []
        elif gene_id_list:
            self.genome = ""
            self.kegg_ids = gene_id_list
        else:
            print('Need to specify either genome or gene list')
            sys.exit()
        self.name = ""
        self.genes = []
        self.pathways = {}
        self.out_dir = out_dir
        self.gene_set = []


    def get_kegg_genes(self):
        try:
            if self.genome:
                url = f"http://rest.kegg.jp/list/{self.genome}"
                entry = requests.get(url).text.strip()
                self.kegg_ids = [gene.split('\t')[0] for gene in entry.split("\n")]
        except HTTPError:
            print("Bad genome identifier")

    def write_gmt(self, file_path, p=True):
        for gene in self.genes:
            name = gene.entry
            if p:
                group = gene.pathways
            else:
                group = gene.modules

            # Go through the list that is saved in the dict:
            for pathway in group:
                # Check if in the inverted dict the key exists
                if pathway not in self.pathways:
                    # If not create a new list

                    self.pathways[pathway] = [name]
                else:
                    self.pathways[pathway].append(name)
        with open(file_path, 'w') as fo:
            for key in self.pathways.keys():
                print(key.split())
                pathway_id = key.split()[0]
                pathway_desc = " ".join(key.split()[1:])
                fo.write(f"{pathway_id.strip()}\t{pathway_desc.strip()}\t" + "\t".join(set(self.pathways[key])) + "\n")
            # massage and write to gmt.
        return self.pathways

    def get_info_for_each_gene(self):

        for gene_id in self.kegg_ids:
            print(gene_id)
            fg = Gene(gene_id)
            fg.kegg_get()
            self.gene_set.append(fg)

    def write_nt_seq(self):
        today = dt.datetime.today().strftime("%Y_%m_%d")
        filename = os.path.join(self.out_dir,  "{}_gene_set_nt_seq.fasta".format(today))
        with open(filename, "w") as fo:
            for gene in self.gene_set:
                fo.write(gene.get_nt_seq())
        return filename

    def write_aa_seq(self):
        today = dt.datetime.today().strftime("%Y_%m_%d")
        filename = os.path.join(self.out_dir,  "{}_gene_set_nt_seq.fasta".format(today))
        with open(filename, "w") as fo:
            for gene in self.gene_set:
                fo.write(gene.get_aa_seq())
        return filename

    def get_info_df(self):
        genes = []
        labels = ["Entry", "Name", "Function", "Pathways"]
        for gene in self.gene_set:
            gene_info = (gene.entry, gene.name, gene.definition, gene.pathway_string)
            genes.append(gene_info)
        return pd.DataFrame.from_records(genes, index="Entry", columns=labels)

    def __str__(self):
        s = ""
        for g in self.gene_set:
            s += g.get_short_info() + "\n"
        return s.strip()


if __name__ == "__main__":
    #gene_list = ["sey:fumA", "sey:SL1344_1398"]
    #gene = Gene()
    #gene.get_gene_info(['sey:fumA'])
    gs = GeneSet(genome="sey",)
    #print(gs.kegg_ids)
    print('Getting genes')
    gs.get_kegg_genes()
    print(gs.kegg_ids[0:10])
    print("getting gene info")
    gs.get_info_for_each_gene()
    print('processing pathway info')
    gs.write_gmt('./sey_pathways.gmt', True)
    #print(gs.pathway_set)


    #gs.write_nt_seq()


# def get_pathways_for_genome(genome="sey"):
#     response = requests.get(f"http://rest.kegg.jp/list/module/{genome}")
#     return response
#
#
#
# def get_pathway_gene_set(pathway="path:sey00010"):
#     response = requests.get(f"http://rest.kegg.jp/get/{pathway}")
#     return response
#
# if __name__ == "__main__":
#     print(get_pathways_for_genome("sey").text)