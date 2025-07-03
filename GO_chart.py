import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.figure(figsize=(30,10))

# parser = argparse.ArgumentParser(prog="GO_chart.py", description="Program for drawing Gene Ontology level 1 graph")
# parser.add_argument("-i", "--input", type=str, help="Input table of the %(prog)s program", required=True)
# parser.add_argument("-o", "--output", type=str, help="Output graph of the %(prog)s program", required=True)
# args = parser.parse_args()

GO_analysis_output = "/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/14.GO/KC2_6_GO_analysis.csv"
# GO_analysis_output = args.input

content = pd.read_csv(GO_analysis_output)
df = pd.DataFrame(content.groupby(["Category"])["Number of genes"].sum())
print(df)
fig = sns.barplot(data=df, x="Number of genes", y=df.index, orient = 'h', palette=["#3274a1", "#e1812c", "#3a923a"])

fig.bar_label(fig.containers[0], size=25)
fig.bar_label(fig.containers[1], size=25)
fig.bar_label(fig.containers[2], size=25)

fig.set_xlabel("Number of genes", fontsize=25, fontweight="bold")
fig.set_ylabel("Function", fontsize=25, fontweight="bold")

fig.tick_params(labelsize=20)

plt.title("GO Enrichment Analysis", fontsize=30, fontweight="bold")
plt.tight_layout()
# plt.savefig("test.png")

# GO_lv1_graph = args.output
plt.savefig("/data18tb/datnguyen/WGS_MinhIBT_KC2_6/Results/14.GO/KC2_6_GO_graph.lv1.png", dpi=400) 
