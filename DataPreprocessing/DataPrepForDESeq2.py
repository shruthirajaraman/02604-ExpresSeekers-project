import os
import pandas as pd

# Define the main folder containing all sample subfolders
main_folder = "../Datasets/BreastCancer"
output_file = "../Datasets/DESeq2Input/breastCancerCounts.csv"

# Dictionary to store expression data
expression_dfs = []  # Store each sample as a DataFrame with gene names

# Iterate over each sample folder
for sample_folder in os.listdir(main_folder):
    sample_path = os.path.join(main_folder, sample_folder)
    
    # Ensure it's a directory
    if os.path.isdir(sample_path):
        
        # Find the RNA-seq data file (excluding the logs folder)
        expression_file = None
        for file in os.listdir(sample_path):
            if file.endswith(".tsv") and "logs" not in file.lower():
                expression_file = os.path.join(sample_path, file)
                break  # We are expecting 1 tsv file

        if expression_file:
            # Read the file (1st row skipped since it contains metadata)
            df = pd.read_csv(expression_file, sep='\t', skiprows=1)

            # Ensure required columns exist
            if "gene_id" not in df.columns or "unstranded" not in df.columns:
                raise ValueError(f"Missing expected columns in {expression_file}.")

            # Rename the expression column to the sample ID (using the folder name)
            df = df[["gene_id", "unstranded"]].rename(columns={"unstranded": sample_folder})

            expression_dfs.append(df)

# Merge all samples using an outer join to keep all genes. We will remove the null elements later
merged_df = expression_dfs[0]
for df in expression_dfs[1:]:
    merged_df = merged_df.merge(df, on="gene_id", how="outer")

merged_df = merged_df[:-4] # Remove the summary statistic rows
merged_df.index = merged_df['gene_id'] # Set index as gene_id
merged_df= merged_df.drop('gene_id', axis = 1) # Remove the gene_id which is present as a separate column
merged_df.dropna() # Remove null unstranded values
merged_df = merged_df.loc[merged_df.iloc[:, 1:].sum(axis=1) > 0] # Remove those genes that have 0 expression across all samples
print("Dataset Head:")
print(merged_df.head())

# Save merged dataset
merged_df.to_csv(output_file, index=True)
# Display summary
print(f"Processed {len(expression_dfs)} samples with {merged_df.shape[0]} unique genes.")