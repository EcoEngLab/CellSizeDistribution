import pandas as pd
import os
import glob
from pathlib import Path

def load_merged_file(file_path):
    """
    Load a merged long-format TSV file.
    
    Parameters:
        file_path (str): Path to the merged TSV file
        
    Returns:
        pandas.DataFrame: Loaded dataframe
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        return None

def remove_columns(df, columns_to_remove):
    """
    Remove specified columns from the dataframe if they exist.
    
    Parameters:
        df (pandas.DataFrame): Input dataframe
        columns_to_remove (list): List of column names to remove
        
    Returns:
        pandas.DataFrame: Dataframe with specified columns removed
    """
    existing_columns = [col for col in columns_to_remove if col in df.columns]
    if existing_columns:
        df = df.drop(columns=existing_columns)
    return df

def get_grouping_columns(df):
    """
    Get the columns to use for grouping (all columns except run_accession and read_count).
    
    Parameters:
        df (pandas.DataFrame): Input dataframe
        
    Returns:
        list: List of column names for grouping
    """
    grouping_cols = [col for col in df.columns if col not in ['run_accession', 'read_count']]
    return grouping_cols

def pool_data(df):
    """
    Pool the data by summing read_count for each unique combination of grouping columns.
    
    Parameters:
        df (pandas.DataFrame): Input dataframe
        
    Returns:
        pandas.DataFrame: Pooled dataframe
    """
    # Get columns for grouping
    grouping_cols = get_grouping_columns(df)
    
    # Group by all columns except run_accession and read_count, then sum read_count
    pooled_df = df.groupby(grouping_cols, as_index=False)['read_count'].sum()
    
    return pooled_df

def process_merged_file(file_path, output_dir):
    """
    Process a single merged TSV file.
    
    Parameters:
        file_path (str): Path to the merged TSV file
        output_dir (str): Output directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Load the merged file
        df = load_merged_file(file_path)
        if df is None:
            return False
        
        # Remove specified columns
        columns_to_remove = ['sample_accession', 'experiment_accession']
        df = remove_columns(df, columns_to_remove)
        
        # Pool the data
        pooled_df = pool_data(df)
        
        # Create output filename
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_name}_pooled_by_study.tsv"
        output_path = os.path.join(output_dir, output_filename)
        
        # Save pooled data
        pooled_df.to_csv(output_path, sep='\t', index=False)
        
        return True
        
    except Exception as e:
        return False

def main():
    """
    Main function to orchestrate the pooling process.
    """
    # Configuration
    input_dir = "output"
    output_dir = "pooled_output"
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # Find all merged TSV files
    pattern = os.path.join(input_dir, "*_with_metadata.tsv")
    merged_files = glob.glob(pattern)
    
    if not merged_files:
        return
    
    # Process each merged file
    successful = 0
    failed = 0
    
    for file_path in merged_files:
        if process_merged_file(file_path, output_dir):
            successful += 1
        else:
            failed += 1

if __name__ == "__main__":
    main()
