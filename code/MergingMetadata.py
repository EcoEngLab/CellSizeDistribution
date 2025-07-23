import pandas as pd
import os
import glob
from pathlib import Path

def load_metadata(metadata_file):
    """
    Load the metadata TSV file containing run_accession and study_accession columns.
    
    Parameters:
        metadata_file (str): Path to the metadata TSV file
        
    Returns:
        pandas.DataFrame: Metadata dataframe
    """
    try:
        metadata = pd.read_csv(metadata_file, sep='\t')
        return metadata
    except Exception as e:
        return None

def identify_sample_columns(df):
    """
    Identify sample columns that start with 'SRR', 'ERR', or 'DRR'.
    
    Parameters:
        df (pandas.DataFrame): Input dataframe
        
    Returns:
        list: List of sample column names
    """
    sample_cols = [col for col in df.columns if col.startswith(('SRR', 'ERR', 'DRR'))]
    return sample_cols

def reshape_to_long_format(df, sample_cols):
    """
    Reshape data from wide to long format.
    
    Parameters:
        df (pandas.DataFrame): Input dataframe
        sample_cols (list): List of sample column names
        
    Returns:
        pandas.DataFrame: Long format dataframe
    """
    # Get non-sample columns (taxonomic info)
    id_cols = [col for col in df.columns if col not in sample_cols]
    
    # Reshape to long format
    long_df = df.melt(
        id_vars=id_cols,
        value_vars=sample_cols,
        var_name='run_accession',
        value_name='read_count'
    )
    
    return long_df

def merge_with_metadata(long_df, metadata_df):
    """
    Merge long format data with metadata on run_accession.
    
    Parameters:
        long_df (pandas.DataFrame): Long format dataframe
        metadata_df (pandas.DataFrame): Metadata dataframe
        
    Returns:
        pandas.DataFrame: Merged dataframe
    """
    # Merge on run_accession, keeping all rows from long_df
    merged_df = pd.merge(
        long_df,
        metadata_df[['run_accession', 'study_accession', 'sample_accession', 'experiment_accession']],
        on='run_accession',
        how='left'
    )
    
    # Fill missing values with empty string
    merged_df['study_accession'] = merged_df['study_accession'].fillna('')
    merged_df['sample_accession'] = merged_df['sample_accession'].fillna('')
    merged_df['experiment_accession'] = merged_df['experiment_accession'].fillna('')
    
    return merged_df

def process_host_species_file(file_path, metadata_df, output_dir):
    """
    Process a single host species TSV file.
    
    Parameters:
        file_path (str): Path to the host species TSV file
        metadata_df (pandas.DataFrame): Metadata dataframe
        output_dir (str): Output directory path
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Load the host species file
        host_df = pd.read_csv(file_path, sep='\t')
        
        # Identify sample columns
        sample_cols = identify_sample_columns(host_df)
        if not sample_cols:
            return False
        
        # Reshape to long format
        long_df = reshape_to_long_format(host_df, sample_cols)
        
        # Merge with metadata
        merged_df = merge_with_metadata(long_df, metadata_df)
        
        # Create output filename
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_filename = f"{base_name}_with_metadata.tsv"
        output_path = os.path.join(output_dir, output_filename)
        
        # Save merged data
        merged_df.to_csv(output_path, sep='\t', index=False)
        
        return True
        
    except Exception as e:
        return False

def main():
    """
    Main function to orchestrate the metadata merging process.
    """
    # Configuration
    metadata_file = "filtered_metadata.tsv"
    host_species_dir = "qiime_taxon_tables_split (1)"
    output_dir = "output"
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # Load metadata
    metadata_df = load_metadata(metadata_file)
    if metadata_df is None:
        return
    
    # Check required columns
    required_cols = ['run_accession', 'study_accession', 'sample_accession', 'experiment_accession']
    missing_cols = [col for col in required_cols if col not in metadata_df.columns]
    if missing_cols:
        return
    
    # Find all host species TSV files
    pattern = os.path.join(host_species_dir, "*_L*.tsv")
    host_files = glob.glob(pattern)
    
    if not host_files:
        return
    
    # Process each host species file
    successful = 0
    failed = 0
    
    for file_path in host_files:
        if process_host_species_file(file_path, metadata_df, output_dir):
            successful += 1
        else:
            failed += 1

if __name__ == "__main__":
    main()
