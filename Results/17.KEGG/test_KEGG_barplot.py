# Draw bar chart for KEGG analysis

import pandas as pd
import matplotlib.pyplot as plt
# import KEGG_analysis

KEGG_lv1 = "/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv1.xlsx"


fig = plt.figure(figsize = (25,10))

plt.rcParams.update({'font.size': 25})
plt.rcParams.update({'font.family':"serif"})

colors = ["#EE4E4E", "#3572EF", "#8576FF", "#41B06E", "#FFD23F", "#B67352"]
df_lv1 = pd.read_excel(KEGG_lv1)
barhs = plt.barh(df_lv1["Pathway"], df_lv1["Number of genes"], color= colors)

for barh, value in zip(barhs, df_lv1["Number of genes"]):
    plt.text(barh.get_width() + 50, barh.get_y() + barh.get_height()/5, value, 
             ha="center", va="bottom", fontsize=20)
# plt.xticks(rotation=45)
plt.title("KEGG Pathway Enrichment", fontsize = 25, fontweight="bold")
plt.xlabel("Number of genes", fontsize = 25, fontweight="bold")
plt.ylabel("Pathway Name", fontsize = 25, fontweight="bold")

plt.tight_layout()
plt.savefig("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv1.png",
            format="png",
            dpi=400)


## Metabolism
KEGG_lv2 = "/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv2.xlsx"
plt.rcParams.update({'font.family':"serif"})
df_lv2 = pd.read_excel(KEGG_lv2)
# print(df_lv2)a

# print(KEGG_analysis.dict_1)
fig = plt.figure(figsize = (20,12))
# print(df_lv2)
colors = ["#EE4E4E", "#3572EF", "#8576FF", "#41B06E", "#FFD23F", "#B67352"]

Metabolism = df_lv2.loc[0:11]
GIP = df_lv2.loc[12:15]
EIP = df_lv2.loc[16:18]
CP = df_lv2.loc[19:22]
OS = df_lv2.loc[23:29]
HD = df_lv2.loc[30:]
# print(GIP)
# # print(OS)
categories = ["Metabolism", "GIP", "EIP", "CP", "OS", "HD"]

meta_barhs = plt.barh(Metabolism["Pathway"], Metabolism["Number of genes"], color=colors[0])

for barh, value in zip(meta_barhs, Metabolism["Number of genes"]):
    plt.text(barh.get_width() + 25, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

# plt.title("Biểu đồ phân nhóm genes theo Human Diseases", fontsize = 16)
# plt.xlabel("Số lượng genes", fontsize = 16)
# plt.ylabel("Nhóm genes", fontsize = 16)

# plt.tight_layout()
# plt.savefig("C:/Users/Admin/Desktop/Python_code/lobi_analysis/Monthon_HD.png", format="png")


GIP_barhs = plt.barh(GIP["Pathway"], GIP["Number of genes"], color=colors[1])

for barh, value in zip(GIP_barhs, GIP["Number of genes"]):
    plt.text(barh.get_width() + 15, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

# plt.title("Biểu đồ phân nhóm genes theo Genetic Information Processing", fontsize = 16)
# plt.xlabel("Số lượng genes", fontsize = 16)
# plt.ylabel("Nhóm genes", fontsize = 16)

# plt.tight_layout()
# plt.savefig("C:/Users/Admin/Desktop/Python_code/lobi_analysis/Ri3_GIP.png", format="png")

EIP_barhs = plt.barh(EIP["Pathway"], EIP["Number of genes"], color=colors[2])

for barh, value in zip(EIP_barhs, EIP["Number of genes"]):
    plt.text(barh.get_width() + 20, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

# plt.title("Biểu đồ phân nhóm genes theo Environmental Information Processing", fontsize = 16)
# plt.xlabel("Số lượng genes", fontsize = 16)
# plt.ylabel("Nhóm genes", fontsize = 16)

# plt.tight_layout()
# plt.savefig("C:/Users/Admin/Desktop/Python_code/lobi_analysis/Ri3_EIP.png", format="png")

CP_barhs = plt.barh(CP["Pathway"], CP["Number of genes"], color=colors[3])

for barh, value in zip(CP_barhs, CP["Number of genes"]):
    plt.text(barh.get_width() + 15, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

# plt.title("Biểu đồ phân nhóm genes theo Cellular Processes", fontsize = 16)
# plt.xlabel("Số lượng genes", fontsize = 16)
# plt.ylabel("Nhóm genes", fontsize = 16)

# plt.tight_layout()
# plt.savefig("C:/Users/Admin/Desktop/Python_code/lobi_analysis/Ri3_CP.png", format="png")

OS_barhs = plt.barh(OS["Pathway"], OS["Number of genes"], color=colors[4])

for barh, value in zip(OS_barhs, OS["Number of genes"]):
    plt.text(barh.get_width() + 15, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

# plt.title("Biểu đồ phân nhóm genes theo Organismal Systems", fontsize = 16)
# plt.xlabel("Số lượng genes", fontsize = 16)
# plt.ylabel("Nhóm genes", fontsize = 16)

# plt.tight_layout()
# plt.savefig("C:/Users/Admin/Desktop/Python_code/lobi_analysis/Ri3_OS.png", format="png")

HD_barhs = plt.barh(HD["Pathway"], HD["Number of genes"], color=colors[5])

for barh, value in zip(HD_barhs, HD["Number of genes"]):
    plt.text(barh.get_width() + 15, barh.get_y() + barh.get_height()/10, value, ha="center", va="bottom", fontsize=15)

plt.title("KEGG Pathway Enrichment", fontsize=20, fontweight="bold")
plt.xlabel("Number of genes", fontsize=16, fontweight="bold")
plt.ylabel("Pathway Name", fontsize=16, fontweight="bold")
# plt.ylabel("Nhóm genes", fontsize = 16)
plt.legend(["Metabolism", 
            "Genetic Information Processing", 
            "Environment Information Processing", 
            "Cellular Processes", 
            "Organismal Systmes",
            "Human Diseases"], title="KEGG Pathway Groups", fontsize=13, title_fontsize=13)
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)

plt.gca().set_axisbelow(True)
# plt.grid(axis="both", alpha=0.3)
plt.tight_layout()
plt.savefig("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/17.KEGG/K2_6_KEGG.pathway.lv2.png", format="png")