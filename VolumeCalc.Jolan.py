import pandas as panda
import numpy as npy
import math


bacFile = panda.read_csv('bacteria.bacdive_cell_data.ver2 - edited.csv')
archaeaFile = panda.read_csv('archaea.bacdive_cell_data.ver2.csv')

currentrownum = 0
pi = math.pi

def volumecalc(thisrow):
    '''
    Purpose:
        Calculate the minimum, maximum and average volume of the current species of bacteria or archea being iterated

    Parameters:
        thisrow: this is the current row in the csv file

    Returns:
        minVolume: minimum volume of cell
        maxVolume: maximum volume of cell
        avgVolume: average of the two 
    '''
    shape = thisrow["Cell shape"]
    minVolume = -1
    maxVolume = -1
    surfaceArea = -1
    minlengthnocaps = 0
    maxlengthnocaps = 0
    try:
        minlengthtotal = thisrow["LengthMin"]
        maxlengthtotal = thisrow["LengthMax"]
        minwidth = thisrow["WidthMin"]
        maxwidth = thisrow["WidthMax"]
    except:
        return panda.Series({"VolumeMin":npy.nan, "VolumeMax":npy.nan, "SurfaceArea":npy.nan})
            
    minradius = minwidth/2
    maxradius = maxwidth/2
    if (shape in ["rod-shaped", "filament-shaped", "vibrio-shaped", "helical-shaped", "spiral-shaped", "curved-shaped", "spore-shaped"]):
        #Im unsure about shape of spore shape it seems to be either rod or ellipse. assumed rod
        #Assumed that caps have same radius as Diameter
        #Assumed that height is equal to width
        #Assumed that filament shape does not taper
        #Approximated vibrio as a rod
        #Assumed that helical and spiral widths are width of cell not width of helix
        #Assumed that helical and spiral lengths are actual cell lenth not spiral length
        minlengthnocaps = minlengthtotal - minwidth #needed for calculating cylinder not including caps
        maxlengthnocaps = maxlengthtotal - maxwidth
        if (minlengthnocaps < 0): # if length is shorter than width, width and length are switched
            minlengthnocaps = minlengthnocaps*-1
            minwidth = minlengthtotal
        if (maxlengthnocaps < 0):
            maxlengthtotal = maxlengthtotal*-1
            maxwidth = maxlengthtotal
        if (minlengthtotal == minwidth):
            minlengthnocaps = 0
        if (maxlengthtotal == maxwidth):
            maxlengthnocaps == 0
        minVolume = ((pi*minradius*minradius)*(((4/3)*minradius)+minlengthnocaps))
        maxVolume = ((pi*maxradius*maxradius)*(((4/3)*maxradius)+maxlengthnocaps))
        surfaceArea = npy.sqrt(((2*pi*minradius)*((2*minradius) + minlengthnocaps))*((2*pi*maxradius)*((2*maxradius) + maxlengthnocaps)))
    elif (shape in ["ellipsoidal", "coccus-shaped", "ovoid-shaped", "oval-shaped", "sphere-shaped", "crescent-shaped"]): #all approximated as ellipse
        minVolume = ((4/3)*pi*minradius*minradius*(minlengthtotal/2)) #Again assuming height and width are equal
        maxVolume = ((4/3)*pi*maxradius*maxradius*(maxlengthtotal/2))
        surfaceArea = npy.sqrt((4*pi*minradius*minradius)*(4*pi*maxradius*maxradius))
    else:
        # For unknown shapes, default to rod-shaped calculation
        minlengthnocaps = minlengthtotal - minwidth
        maxlengthnocaps = maxlengthtotal - maxwidth
        if (minlengthnocaps < 0):
            minlengthnocaps = minlengthnocaps*-1
            minwidth = minlengthtotal
        if (maxlengthnocaps < 0):
            maxlengthtotal = maxlengthtotal*-1
            maxwidth = maxlengthtotal
        if (minlengthtotal == minwidth):
            minlengthnocaps = 0
        if (maxlengthtotal == maxwidth):
            maxlengthnocaps == 0
        minVolume = ((pi*minradius*minradius)*(((4/3)*minradius)+minlengthnocaps))
        maxVolume = ((pi*maxradius*maxradius)*(((4/3)*maxradius)+maxlengthnocaps))
        surfaceArea = npy.sqrt(((2*pi*minradius)*((2*minradius) + minlengthnocaps))*((2*pi*maxradius)*((2*maxradius) + maxlengthnocaps)))
    
    if (surfaceArea > 0 and minVolume > 0 and maxVolume > 0):
        return panda.Series({"VolumeMin":minVolume, "VolumeMax":maxVolume, "SurfaceArea":surfaceArea})
    else:
        return panda.Series({"VolumeMin":npy.nan, "VolumeMax":npy.nan, "SurfaceArea":npy.nan})

def get_dominant_shape_by_genus(df, genus_col="Genus", shape_col="Cell shape"):
    """
    Find the most dominant cell shape for each genus
    
    Parameters:
        df: DataFrame containing the data
        genus_col: Column name for genus
        shape_col: Column name for cell shape
    
    Returns:
        Dictionary mapping genus to dominant shape
    """
    dominant_shapes = {}
    
    for genus in df[genus_col].dropna().unique():
        genus_data = df[df[genus_col] == genus]
        # count shapes, excluding NaN/empty values
        shape_counts = genus_data[shape_col].dropna().value_counts()
        if len(shape_counts) > 0:
            dominant_shapes[genus] = shape_counts.index[0]  # most common shape
        else:
            dominant_shapes[genus] = "rod-shaped"  # default to rod-shaped
    
    return dominant_shapes

def fill_missing_volumes_with_dominant_shape(df, dominant_shapes_dict, genus_col="Genus"):
    """
    Fill missing volumes using the dominant shape from the same genus
    
    Parameters:
        df: DataFrame containing the data
        dominant_shapes_dict: Dictionary mapping genus to dominant shape
        genus_col: Column name for genus
    
    Returns:
        DataFrame with filled volumes
    """
    # identify rows with missing volumes
    missing_volumes_mask = df["VolumeMin"].isna() | df["VolumeMax"].isna()
    
    df_filled = df.copy()
    
    filled_count = 0
    for idx in df_filled[missing_volumes_mask].index:
        row = df_filled.loc[idx]
        genus = row[genus_col]
        
        if genus in dominant_shapes_dict:
            temp_row = row.copy()
            temp_row["Cell shape"] = dominant_shapes_dict[genus]
            
            volume_result = volumecalc(temp_row)
            
            # update the row if calculation was successful
            if not (panda.isna(volume_result["VolumeMin"]) or panda.isna(volume_result["VolumeMax"])):
                df_filled.loc[idx, "VolumeMin"] = volume_result["VolumeMin"]
                df_filled.loc[idx, "VolumeMax"] = volume_result["VolumeMax"]
                df_filled.loc[idx, "SurfaceArea"] = volume_result["SurfaceArea"]
                filled_count += 1
        else:
            # if no dominant shape found, use rod-shaped as default
            temp_row = row.copy()
            temp_row["Cell shape"] = "rod-shaped"
            
            volume_result = volumecalc(temp_row)
            
            if not (panda.isna(volume_result["VolumeMin"]) or panda.isna(volume_result["VolumeMax"])):
                df_filled.loc[idx, "VolumeMin"] = volume_result["VolumeMin"]
                df_filled.loc[idx, "VolumeMax"] = volume_result["VolumeMax"]
                df_filled.loc[idx, "SurfaceArea"] = volume_result["SurfaceArea"]
                filled_count += 1
    
    return df_filled

# initial volume calculation
bacFile[["VolumeMin", "VolumeMax", "SurfaceArea"]] = bacFile.apply(volumecalc, axis=1, result_type='expand')
archaeaFile[["VolumeMin", "VolumeMax", "SurfaceArea"]] = archaeaFile.apply(volumecalc, axis=1, result_type='expand')

# get dominant shapes for each genus
bac_dominant_shapes = get_dominant_shape_by_genus(bacFile)
arch_dominant_shapes = get_dominant_shape_by_genus(archaeaFile)

# fill missing volumes using dominant shapes
bacFile = fill_missing_volumes_with_dominant_shape(bacFile, bac_dominant_shapes)
archaeaFile = fill_missing_volumes_with_dominant_shape(archaeaFile, arch_dominant_shapes)

# calculate average and geometric mean volumes
bacFile["AvgVolume"] = bacFile[["VolumeMin", "VolumeMax"]].mean(axis=1, skipna=True)
validbac = (bacFile["VolumeMin"] >0) & (bacFile["VolumeMax"]>0)
bacFile.loc[validbac, "GeoMeanVolume"] = npy.sqrt(bacFile["VolumeMin"] * bacFile["VolumeMax"])
bacFile.loc[~validbac, "GeoMeanVolume"] = npy.nan

archaeaFile["AvgVolume"] = archaeaFile[["VolumeMin", "VolumeMax"]].mean(axis=1, skipna=True)
validarc = (archaeaFile["VolumeMin"] >0) & (archaeaFile["VolumeMax"]>0)
archaeaFile.loc[validarc, "GeoMeanVolume"] = npy.sqrt(archaeaFile["VolumeMin"] * archaeaFile["VolumeMax"])
archaeaFile.loc[~validarc, "GeoMeanVolume"] = npy.nan

bacFile.to_csv("VolumeOutputBacteria.csv", index=False)
archaeaFile.to_csv("VolumeOutputArchaea.csv", index=False)