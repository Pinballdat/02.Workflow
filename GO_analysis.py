# Import the OBO parser from GOATools
# from goatools import obo_parser
# import wget
# import os
import pandas as pd
from collections import defaultdict
# from collections import OrderedDict
from itertools import islice
import argparse

GO_basic_db = "/data18tb/datnguyen/[WGS]_MinhIBT/Results/TT2_GO/go-basic.obo"
with open(GO_basic_db, "r") as hd:
    content = hd.read().split("\n\n[Term]\n")[1:]
# print(content)
bio_pro = {}
mol_fun = {}
cel_com = {}
Name_sp_db = {}
for c in content:
    if len(c.split("\n")) >= 3:
        Name_sp_db[c.split("\n")[0].split(": ")[1]] = c.split("\n")[2].split(": ")[1]
        if c.split("\n")[2].split(": ")[1] == "biological_process":
            bio_pro[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]
        elif c.split("\n")[2].split(": ")[1] == "molecular_function":
            mol_fun[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]
        elif c.split("\n")[2].split(": ")[1] == "cellular_component":
            cel_com[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]

# parser = argparse.ArgumentParser(prog="GO_analysis.py", description="Program for analysing Gene Ontology function in Eggnog_mapper")
# parser.add_argument("-i", "--input", type=str, help="Input table of %(prog)s program", required=True)
# parser.add_argument("-o", "--output", type=str, help="Output table of %(prog)s program", required=True)

# args = parser.parse_args()

# GO_analysis_input = args.input
GO_analysis_input = "/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/14.GO/KC2_6_GO_id.txt"
# GO_results = "/data18tb/datnguyen/[WGS]_MinhIBT/Results/TT2_GO/TT2_GO_emapper.tsv"
GO_df = pd.read_csv(GO_analysis_input, sep="\t")
GO_dict = {}
for i in GO_df.index:
    GO_dict[GO_df.loc[i,"#query"]] = GO_df.loc[i,"GOs"].split(",")

new_list = []
for value in GO_dict.values():
    for v in value:
        if v != "-":
            new_list.append(v)
            
BioPro_dict = defaultdict(int)
MolFun_dict = defaultdict(int)
CelCom_dict = defaultdict(int)
for element in new_list:
    if element in bio_pro.keys():
        BioPro_dict[bio_pro[element]] += 1
    if element in mol_fun.keys():
        MolFun_dict[mol_fun[element]] += 1
    if element in cel_com.keys():
        CelCom_dict[cel_com[element]] += 1
sorted_BioPro_dict = dict(sorted(BioPro_dict.items(), key=lambda item: item[1]))
sorted_MolFun_dict = dict(sorted(MolFun_dict.items(), key=lambda item: item[1]))
sorted_CelCom_dict = dict(sorted(CelCom_dict.items(), key=lambda item: item[1]))
# print(sorted_BioPro_dict)
merged_dict = {**sorted_BioPro_dict, **sorted_MolFun_dict, **sorted_CelCom_dict}
sorted_merged_dict = dict(sorted(merged_dict.items(), key=lambda item: item[1]))
# print(merged_dict)
# print(len(sorted_CelCom_dict))
# final_dict = dict(islice(reversed(sorted_merged_dict.items())))
final_df = pd.DataFrame({"Number of genes": sorted_merged_dict})
for i in final_df.index:
    if i in sorted_BioPro_dict.keys():
        final_df.loc[i, "Category"] = "Biological Process"
    if i in sorted_MolFun_dict.keys():
        final_df.loc[i, "Category"] = "Molecular Function"
    if i in sorted_CelCom_dict.keys():
        final_df.loc[i, "Category"] = "Cellular Component"
#         print(final_df.loc[i,"Number of genes"])
sorted_final_df = final_df.sort_values(["Category","Number of genes"], ascending=[True, True])
print(sorted_final_df)
# GO_analysis_output = args.output
sorted_final_df.to_csv("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/14.GO/KC2_6_GO_analysis.csv")

