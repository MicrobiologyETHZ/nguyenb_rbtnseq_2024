import aiohttp
import asyncio
import time
from kegg_api import Gene
import re
import requests
from requests.exceptions import HTTPError
import sys


async def get_gene(session, url):
    async with session.get(url) as resp:
        gene_text = await resp.text()
        gene = Gene()
        lines = gene_text.strip().split("\n")
        for i, line in enumerate(lines):
            if len(line.strip()) < 1:
                continue
            field = line.strip().split()[0]
            print(field)
            if field == 'ENTRY':
                gene.entry = lines[i].split()[1]
            elif field == "NAME":
                gene.name = lines[i].split("NAME        ")[1]
            elif field == "DEFINITION":
                p = lines[i].split("DEFINITION  ")[1]
                gene.definition = re.sub(r"^\(.*?\)", "", p).strip()
            elif field == "PATHWAY":
                p = lines[i].split("PATHWAY     ")[1]
                count = 1
                while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                    p += lines[i + count]
                    count += 1
                gene.pathways = p.split("            ")
                gene.pathway_string = ", ".join(gene.pathways).rstrip(", ")
            elif field == "MODULE":
                p = lines[i].split("MODULE     ")[1]
                count = 1
                while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                    p += lines[i + count]
                    count += 1
                gene.modules = p.split("            ")
                gene.module_string = ", ".join(gene.pathways).rstrip(", ")
            elif field == "AASEQ":
                p = ''
                count = 1
                while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                    p += lines[i + count]
                    count += 1
                gene.aa_seq = p.replace("            ", '')
            elif field == "NTSEQ":
                p = ''
                count = 1
                while lines[i + count].startswith("           ") and len(lines[i + count]) > 0:
                    p += lines[i + count]
                    count += 1
                gene.nt_seq = p.replace("            ", '')
        return gene


async def main(gene_list):
    async with aiohttp.ClientSession() as session:
        tasks = []
        for gene_id in gene_list:
            url = f'http://rest.kegg.jp/get/{gene_id}'
            tasks.append(asyncio.ensure_future(get_gene(session, url)))
        original_gene = await asyncio.gather(*tasks)
        for gene in original_gene:
            print(gene.name)
    return original_gene


def write_gmt(gene_set, file_path, set_by="pathway"):
    pathways = {}

    for gene in gene_set:
        print(gene.name)
        if set_by == "pathway":
            group = gene.pathways
        else:
            group = gene.modules
        # Go through the list that is saved in the dict:
        for pathway in group:
            # Check if in the inverted dict the key exists
            if pathway not in pathways:
                # If not create a new list
                pathways[pathway] = [gene.entry]
            else:
                pathways[pathway].append(gene.entry)
    with open(file_path, 'w') as fo:
        for key in pathways.keys():
            pathway_id = key.split()[0]
            pathway_desc = " ".join(key.split()[1:])
            fo.write(f"{pathway_id.strip()}\t{pathway_desc.strip()}\t" + "\t".join(set(pathways[key])) + "\n")
        # massage and write to gmt.
    return pathways


if __name__ == "__main__":
    # Get gene list given genome
    genome = sys.argv[1]
    out_file = sys.argv[2]
    set_by = sys.argv[3]

    try:
        print(f"Retrieving genes for {genome}")
        start_time = time.time()
        url = f"http://rest.kegg.jp/list/{genome}"
        entry = requests.get(url).text.strip()
        gene_list = [gene.split('\t')[0] for gene in entry.split("\n")]
        print(f"Found {len(gene_list)} genes")
        print("--- %s seconds ---" % (time.time() - start_time))
        #for gene_list_chunk in gene_list[0:100]
        print("Retrieving pathway annotations for each gene")
        start_time = time.time()
        gene_set = asyncio.run(main(gene_list=gene_list))
        print(gene_set[0])
        # print("--- %s seconds ---" % (time.time() - start_time))
        # start_time = time.time()
        # print(f"Processing pathway information and writing gmt file to {out_file}")
        # write_gmt(gene_set, out_file, set_by)
        # print("--- %s seconds ---" % (time.time() - start_time))

    except HTTPError:
        print("Bad genome identifier")
