import pandas as pd
import random

# Function to extract the first 12 characters of the sample id to get the patient id
def extract_patient_id(patient_id):
    return patient_id[:12]

def process_networks(df, output_path):
    '''
    Transforms networks in the input format for PARADYS

    Parameters
    ----------
    df : DataFrame
        Original file with all networks
    output_path : string
        Path where the processed networks files will be stored

    Returns
    -------
    None.

    '''
    
    # Extract the first 12 characters from the "patient id" column
    df["patient"] = df["patient id"].apply(extract_patient_id)
    
    # Create an empty list to store the transformed data
    transformed_data = []
    
    # Iterate over the rows of the original dataframe
    for _, row in df.iterrows():
        patient = row["patient"]
        
        # Iterate over the columns representing gene pairs (edges)
        for col_name, edge in row.items():
            # Skip the "patient id" and "patient" columns
            if col_name in ["patient id", "patient"]:
                continue
            
            # Extract gene1 and gene2 from the column name (edge)
            gene1, gene2 = eval(col_name)  # Convert the column name string to a tuple
            
            # If the edge is present (value is 1), add it to the transformed_data
            if edge > 0:
                transformed_data.append({
                    "patient": patient,
                    "edge": col_name,
                    "tf": gene1,
                    "gene": gene2
                })
    
    # Create a new dataframe from the transformed data
    transformed_dataframe = pd.DataFrame(transformed_data)
    
    
    # Save the dataframe to a CSV file
    output_file = output_path+"networks.csv"
    transformed_dataframe.to_csv(output_file, index=False)
    
    print(f"CSV file '{output_file}' has been created.")

def process_mutations(cnvs,muts,output_path):
    '''
    Transforms cnvs and mutations in the input format for PARADYS

    Parameters
    ----------
    cnvs : DataFrame
        
    muts : DataFrame
   
    output_path : String
        
    '''
    
    # Create an empty list to store the transformed data
    transformed_data = []
    
    # Iterate over the rows of the original dataframe
    for _, row in cnvs.iterrows():
        gene_id = row["Gene Symbol"]
        
        # Iterate over the columns representing patient ids
        for col_name, value in row.items():
            # Skip the "gene Symbol" column
            if col_name == "gene Symbol":
                continue
            
            # Extract the patient id and value for the gene
            patient_id = col_name
            gene_value = value
            
            # Check if the value is 1 or -1, then add to the transformed_data
            if gene_value == 1 or gene_value == -1:
                transformed_data.append({
                    "patient": extract_patient_id(patient_id),
                    "gene": gene_id
                })
    
    # Create a new dataframe from the transformed data
    transformed_dataframe = pd.DataFrame(transformed_data)    
    
    # Filter the matrix to include only non-synonymous_variant effects
    transformed_data = muts[muts["effect"] != "synonymous_variant"]
    
    # Extract the first 12 characters from the "Sample ID" column
    transformed_data["patient"] = transformed_data["Sample_ID"].apply(extract_patient_id)
    
    # Create the final dataframe with selected columns
    transformed_dataframe2 = transformed_data[["patient", "gene"]]
    
    # Concatenate the two dataframes
    combined_dataframe = pd.concat([transformed_dataframe, transformed_dataframe2])
    
    # Drop duplicates if any (to keep only one instance of each patient-gene pair)
    combined_dataframe.drop_duplicates(subset=['patient', 'gene'], inplace=True)
    
    combined_dataframe.to_csv(output_path+'mutations.csv', index=False)

# Function to generate random mutations  
def process_mutations_randoms(cnvs, muts, output_path, k, i):
    '''
    Transforms false cnvs and mutations in the input format for PARADYS 

    Parameters
    ----------
    cnvs : DataFrame
        DataFrame containing copy number variations (CNVs).
        
    muts : DataFrame
        DataFrame containing mutations.
   
    output_path : String
        Path to the output directory where the transformed data will be saved.
    
    k : float
        Percentage of random additional mutations to add to the DataFrame.

    i : int
    ID of decoy mutations database

    '''
    # Create an empty list to store the transformed data
    transformed_data = []
    
    # Iterate over the rows of the original dataframe
    for _, row in cnvs.iterrows():
        gene_id = row["Gene Symbol"]
        
        # Iterate over the columns representing patient ids
        for col_name, value in row.items():
            # Skip the "gene Symbol" column
            if col_name == "gene Symbol":
                continue
            
            # Extract the patient id and value for the gene
            patient_id = col_name
            gene_value = value
            
            # Check if the value is 1 or -1, then add to the transformed_data
            if gene_value == 1 or gene_value == -1:
                transformed_data.append({
                    "patient": extract_patient_id(patient_id),
                    "gene": gene_id
                })
    
    # Filter the matrix to include only non-synonymous_variant effects
    transformed_data = muts[muts["effect"] != "synonymous_variant"]
    
    # Extract the first 12 characters from the "Sample ID" column
    transformed_data["patient"] = transformed_data["Sample_ID"].apply(extract_patient_id)
    
    # Create the final dataframe with selected columns
    transformed_dataframe2 = transformed_data[["patient", "gene"]]
    all_genes=transformed_dataframe2['gene']
    
    result_df=pd.DataFrame()
    for patient, group in transformed_dataframe2.groupby('patient'):
        # Extract the list of genes for this patient
        original_genes = group['gene'].tolist()
    
        # Calculate the number of random genes to select for this patient
        num_genes_to_select = int(k * len(original_genes))
    
        # Randomly select unique genes from all_genes
        selected_genes = random.sample(set(all_genes) - set(original_genes), num_genes_to_select)
    
        # Create a new DataFrame for this patient
        patient_df = pd.DataFrame({f'gene_{i}': selected_genes[i] for i in range(num_genes_to_select)})
    
        # Add the 'patient' column and concatenate it with the new DataFrame
        patient_df['patient'] = patient
        result_df = pd.concat([result_df, patient_df], axis=0, ignore_index=True)
        
    # Save the combined dataframe to a CSV file
    result_df.to_csv(output_path + 'mutations '+k+' '+str(i)+' .csv', index=False)

if __name__ == '__main__':    
    
    # Specify input files paths.
    input_path = './'
    output_path = './'
    fea_file_path = input_path+'PRAD_nets_SSN.fea'
    cnvs_file_path = input_path+'TCGA-PRAD_CNVS.tsv'
    snv_file_path = input_path+'TCGA-PRAD.muse_snv.tsv'
    
    # Load input data.
    df = pd.read_feather(fea_file_path)
    cnvs=pd.read_csv(cnvs_file_path,sep='\t')
    muts = pd.read_csv(snv_file_path,sep='\t')
    
    print("Processing mutations...")
    process_mutations(cnvs, muts, output_path)
    print("Processing networks...")
    process_networks(df, output_path)
