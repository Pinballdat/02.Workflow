import pandas as pd

file = "/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.txt"
with open(file, "r") as hd:
    content = hd.readlines()
index_1 = []
index_2 = []
index_3 = []
for index in range(len(content)):
    if content[index][0].isalpha() == True and content[index + 1][0].isalpha() == True:
        index_1.append(index)
    if content[index][0].isalpha() == True and content[index + 1][0].isdigit() == True:
        index_2.append(index)
    if index not in index_1 and index not in index_2:
        index_3.append(index)

dict_1 = {}
dict_2 = {}
valid_gr_2 = []

for i in range(len(index_2)):
    dict_2[content[index_2[i]].strip("\n")] = []
    for j in range(len(index_3)):
        if i < len(index_2)-1:
            if index_3[j] > index_2[i] and index_3[j] < index_2[i+1]:
                dict_2[content[index_2[i]].strip("\n")].append(int(str(content[index_3[j]].strip("\n").split(" ")[-1].strip("()"))))
                valid_gr_2.append(index_3[j])

remaining_gr_2 = [j for j in index_3 if j not in valid_gr_2]
for gr in remaining_gr_2:
    dict_2[content[index_2[-1]].strip("\n")].append(int(str(content[gr].strip("\n")).split(" ")[-1].strip("()")))

Level_2 = []
Genes_lv2 = []

for k in dict_2:
    Level_2.append(k)
    Genes_lv2.append(sum(dict_2[k]))

df_2 = pd.DataFrame({"Pathway":Level_2, "Number of genes":Genes_lv2}).set_index("Pathway")

valid_gr_1 = []
for i in range(len(index_1)):
    dict_1[content[index_1[i]].strip("\n")] = []
    for j in range(len(index_2)):
        if i < len(index_1)-1:
            if index_2[j] > index_1[i] and index_2[j] < index_1[i+1]:
                dict_1[content[index_1[i]].strip("\n")].append(sum(dict_2[str(content[index_2[j]].strip("\n"))]))
                valid_gr_1.append(index_2[j])
                
remaining_gr_1 = [j for j in index_2 if j not in valid_gr_1]
for gr in remaining_gr_1:
    dict_1[content[index_1[-1]].strip("\n")].append(sum(dict_2[str(content[gr].strip("\n"))]))

Level_1 = []
Genes_lv1 = []

for k in dict_1:
    Level_1.append(k)
    Genes_lv1.append(sum(dict_1[k]))

df_1 = pd.DataFrame({"Pathway":Level_1, "Number of genes":Genes_lv1}).set_index("Pathway")



### Ouput: XLSX file
df_1.to_excel("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv1.xlsx")
df_2.to_excel("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv2.xlsx")