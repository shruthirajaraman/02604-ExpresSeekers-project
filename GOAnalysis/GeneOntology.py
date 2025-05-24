# Functional Enrichment Analysis
# Gene Ontology (GO) and KEGG
# Shruthi Rajaraman

# Import Statements
import pandas as pd
import gseapy as gp
from gseapy import dotplot
import matplotlib.pyplot as plt
import os
import mygene

# ReadFile
# INPUT: a string file directory to a file
# OUTPUT: a list containing the contents of filename
def ReadFile(filename):
    f = open(filename)
    genes = f.read()
    f.close()
    return genes

# MakeNewFolder
# INPUT: a string representing a folder name
# OUTPUT: a new directory with folder_name as the name of the folder
def MakeNewFolder(folder_name): 
    try: # Folder created if it doesn't already exist
        os.mkdir(folder_name)
        print("Folder %s created!" % folder_name)
    except FileExistsError: # Error raised if folder already exists
        print("Folder %s already exists" % folder_name)

# Convert EnsemblIDs into gene symbols
# GetGeneSymbol
# INPUT: a list of EnsemblIDs representing a set of genes
# OUTPUT: a dictionary mapping each EnsemblID to its corresponding gene symbol
def GetGeneSymbol(ensemblIDs):
    mg = mygene.MyGeneInfo()
    newEnsemblIDs = []
    geneSymbol_dict = {}

    # Strip the version number from each EnsemblID
    for id in ensemblIDs:
        newID = id.split(".")[0]
        newEnsemblIDs.append(newID)

    # Obtain the entry associated with each EnsemblID
    results = mg.querymany(newEnsemblIDs, scopes='ensembl.gene', fields='symbol', species='human')

    # Create a dictionary mapping each EnsemblID to its corresponding gene symbol
    for entry in results:
        geneSymbol_dict[entry['query']] = entry.get('symbol', None)
    
    return geneSymbol_dict

# ObtainDictValues
# INPUT: a dictionary dict
# OUTPUT: a list containing all the values in dict
def ObtainDictValues(dict):
    values = []
    for key in dict:
        values.append(dict[key])
    return values

# # Obtain all available gene sets
# names = gp.get_library_name()
# for name in names:
#     print(name, "\n")

# Perform Gene Ontology on a specified set of genes
# Read file and initialize list
print("\nReading file and obtaining gene symbols...")
geneFile = ReadFile("intersection_sig_genes.txt").split("\n")
geneDict = GetGeneSymbol(geneFile)
geneSymbol_list = ObtainDictValues(geneDict)

cleaned_geneSymbol_list = []
for symbol in geneSymbol_list:
    if symbol != None:
        cleaned_geneSymbol_list.append(symbol)

# LUAD Upregulated Genes
geneFile_luad = ReadFile("LUAD_upregulated_genes.txt").split("\n")
geneDict_luad = GetGeneSymbol(geneFile_luad)
geneSymbol_list_luad = ObtainDictValues(geneDict_luad)

cleaned_geneSymbol_list_luad = []
for symbol in geneSymbol_list_luad:
    if symbol != None:
        cleaned_geneSymbol_list_luad.append(symbol)

# LUSC Upregulated Genes
geneFile_lusc = ReadFile("LUSC_upregulated_genes.txt").split("\n")
geneDict_lusc = GetGeneSymbol(geneFile_lusc)
geneSymbol_list_lusc = ObtainDictValues(geneDict_lusc)

cleaned_geneSymbol_list_lusc = []
for symbol in geneSymbol_list_lusc:
    if symbol != None:
        cleaned_geneSymbol_list_lusc.append(symbol)

# Use the Enrichr module to perform gene ontology analysis on the gene symbols
print("\nPerforming gene ontology analysis...")
# Molecular Function
enr1 = gp.enrichr(gene_list=cleaned_geneSymbol_list,
                 gene_sets=['GO_Molecular_Function_2025'],
                 organism='human',
                )
enr1_luad = gp.enrichr(gene_list=cleaned_geneSymbol_list_luad,
                 gene_sets=['GO_Molecular_Function_2025'],
                 organism='human',
                )
enr1_lusc = gp.enrichr(gene_list=cleaned_geneSymbol_list_lusc,
                 gene_sets=['GO_Molecular_Function_2025'],
                 organism='human',
                )

# Biological Processes
enr2 = gp.enrichr(gene_list=cleaned_geneSymbol_list,
                 gene_sets=['GO_Biological_Process_2025'],
                 organism='human',
                )
enr2_luad = gp.enrichr(gene_list=cleaned_geneSymbol_list_luad,
                 gene_sets=['GO_Biological_Process_2025'],
                 organism='human',
                )
enr2_lusc = gp.enrichr(gene_list=cleaned_geneSymbol_list_lusc,
                 gene_sets=['GO_Biological_Process_2025'],
                 organism='human',
                )

# KEGG Pathways
enr3 = gp.enrichr(gene_list=cleaned_geneSymbol_list,
                 gene_sets=['KEGG_2021_Human'],
                 organism='human',
                )
enr3_luad = gp.enrichr(gene_list=cleaned_geneSymbol_list_luad,
                 gene_sets=['KEGG_2021_Human'],
                 organism='human',
                )
enr3_lusc = gp.enrichr(gene_list=cleaned_geneSymbol_list_lusc,
                 gene_sets=['KEGG_2021_Human'],
                 organism='human',
                )

# Print gene symbols and GO analysis results to the terminal
print("\nMolecular Function:")
print(enr1.results.head())
print("Molecular Function (LUAD):")
print(enr1_luad.results.head())
print("Molecular Function (LUSC):")
print(enr1_lusc.results.head())

print("\nBiological Processes:")
print(enr2.results.head())
print("Biological Processes (LUAD):")
print(enr2_luad.results.head())
print("Biological Processes (LUSC):")
print(enr2_lusc.results.head())

print("\nKEGG Pathways:")
print(enr3.results.head())
print("KEGG Pathways (LUAD):")
print(enr3_luad.results.head())
print("KEGG Pathways (LUSC):")
print(enr3_lusc.results.head())

# Make folder to store GO results
print("\nWriting GO results to CSV files...")
MakeNewFolder("GO_Results")

# Create CSV files for each GO term

filePath = os.path.join("GO_Results", "molecular_function.csv")
enr1.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "molecular_function_luad.csv")
enr1_luad.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "molecular_function_lusc.csv")
enr1_lusc.results.to_csv(filePath)

filePath = os.path.join("GO_Results", "biological_process.csv")
enr2.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "biological_process_luad.csv")
enr2_luad.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "biological_process_lusc.csv")
enr2_lusc.results.to_csv(filePath)

filePath = os.path.join("GO_Results", "kegg.csv")
enr3.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "kegg_luad.csv")
enr3_luad.results.to_csv(filePath)
filePath = os.path.join("GO_Results", "kegg_lusc.csv")
enr3_lusc.results.to_csv(filePath)

# Plotting Bar Graphs
print("\nCreating bar graphs...")
MakeNewFolder("GO_Plots")

# Molecular Function Bar Graph
filePath = os.path.join("GO_Plots", "molecular_function_bar.png")
enr1.results = enr1.results[:10].sort_values(by="Combined Score")
geneInfo1 = enr1.results["Term"]
combinedScores1 = enr1.results["Combined Score"]
plt.barh(geneInfo1, combinedScores1, color="teal")
plt.title("Molecular Function")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "molecular_function_luad_bar.png")
enr1_luad.results = enr1_luad.results[:10].sort_values(by="Combined Score")
geneInfo1_luad = enr1_luad.results["Term"]
combinedScores1_luad = enr1_luad.results["Combined Score"]
plt.clf()
plt.barh(geneInfo1_luad, combinedScores1_luad, color="teal")
plt.title("Molecular Function (Upregulated Genes in LUAD)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "molecular_function_lusc_bar.png")
enr1_lusc.results = enr1_lusc.results[:10].sort_values(by="Combined Score")
geneInfo1_lusc = enr1_lusc.results["Term"]
combinedScores1_lusc = enr1_lusc.results["Combined Score"]
plt.clf()
plt.barh(geneInfo1_lusc, combinedScores1_lusc, color="teal")
plt.title("Molecular Function (Upregulated Genes in LUSC)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

# Biological Process Bar Graph
filePath = os.path.join("GO_Plots", "biological_process_bar.png")
enr2.results = enr2.results[:10].sort_values(by="Combined Score")
geneInfo2 = enr2.results["Term"]
combinedScores2 = enr2.results["Combined Score"]
plt.clf()
plt.barh(geneInfo2, combinedScores2, color="orchid")
plt.title("Biological Process")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "biological_process_luad_bar.png")
enr2_luad.results = enr2_luad.results[:10].sort_values(by="Combined Score")
geneInfo2_luad = enr2_luad.results["Term"]
combinedScores2_luad = enr2_luad.results["Combined Score"]
plt.clf()
plt.barh(geneInfo2_luad, combinedScores2_luad, color="orchid")
plt.title("Biological Process (Upregulated Genes in LUAD)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "biological_process_lusc_bar.png")
enr2_lusc.results = enr2_lusc.results[:10].sort_values(by="Combined Score")
geneInfo2_lusc = enr2_lusc.results["Term"]
combinedScores2_lusc = enr2_lusc.results["Combined Score"]
plt.clf()
plt.barh(geneInfo2_lusc, combinedScores2_lusc, color="orchid")
plt.title("Biological Process (Upregulated Genes in LUSC)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

# KEGG Bar Graph
filePath = os.path.join("GO_Plots", "KEGG_bar.png")
enr3.results = enr3.results[:10].sort_values(by="Combined Score")
geneInfo3 = enr3.results["Term"]
combinedScores3 = enr3.results["Combined Score"]
plt.clf()
plt.barh(geneInfo3, combinedScores3, color="coral")
plt.title("KEGG")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "KEGG_luad_bar.png")
enr3_luad.results = enr3_luad.results[:10].sort_values(by="Combined Score")
geneInfo3_luad = enr3_luad.results["Term"]
combinedScores3_luad = enr3_luad.results["Combined Score"]
plt.clf()
plt.barh(geneInfo3_luad, combinedScores3_luad, color="coral")
plt.title("KEGG (Upregulated Genes in LUAD)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

filePath = os.path.join("GO_Plots", "KEGG_lusc_bar.png")
enr3_lusc.results = enr3_lusc.results[:10].sort_values(by="Combined Score")
geneInfo3_lusc = enr3_lusc.results["Term"]
combinedScores3_lusc = enr3_lusc.results["Combined Score"]
plt.clf()
plt.barh(geneInfo3_lusc, combinedScores3_lusc, color="coral")
plt.title("KEGG (Upregulated Genes in LUSC)")
plt.xlabel("Combined Score")
plt.ylabel("Term")
plt.savefig(filePath, bbox_inches='tight')

# Plotting Colormaps
print("\nCreating colormaps...")

# Molecular Function Colormap
filePath = os.path.join("GO_Plots", "molecular_function_cmap.png")
ax1 = dotplot(enr1.res2d, title='GO Molecular Function',cmap='viridis_r', size=5, figsize=(3,5))
fig1 = ax1.get_figure()
fig1.savefig(filePath, bbox_inches='tight')

# Biological Process Colormap
filePath = os.path.join("GO_Plots", "biological_process_cmap.png")
ax2 = dotplot(enr2.res2d, title='GO Biological Process',cmap='viridis_r', size=5, figsize=(3,5))
fig2 = ax2.get_figure()
fig2.savefig(filePath, bbox_inches='tight')

# KEGG Colormap
filePath = os.path.join("GO_Plots", "KEGG_cmap.png")
ax3 = dotplot(enr3.res2d, title='KEGG',cmap='viridis_r', size=5, figsize=(3,5))
fig3 = ax3.get_figure()
fig3.savefig(filePath, bbox_inches='tight')

print("\nGO analysis has successfully been completed!")