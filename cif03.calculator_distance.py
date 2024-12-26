import numpy as np
import argparse
import math
import pandas as pd

# Dictionary mapping spectral types to their absolute magnitudes
spectral_types = {
    70.0: 7.98,
    70.5: 8.33,
    71.0: 8.61,
    71.5: 8.86,
    72.0: 9.23,
    72.5: 9.65,
    73.0: 10.08,
    73.5: 10.51,
    74.0: 10.88,
    74.5: 11.44,
    75.0: 12.25,
    75.5: 12.88,
    76.0: 13.64,
    76.5: 14.42,
    77.0: 14.55,
    77.5: 15.01,
    78.0: 15.19,
    78.5: 15.76,
    79.0: 16.11,
    79.5: 16.39,
    80.0: 16.51
}

def photometric_distance(apparent_magnitude, spectral_type, magnitude_uncertainty):
    """
    Calculate the photometric distance of a star and its uncertainty.
    
    Parameters:
    - apparent_magnitude (float): The apparent magnitude of the star.
    - spectral_type (float): The spectral type of the star in numerical format.
    - magnitude_uncertainty (float): The uncertainty in the apparent magnitude.
    
    Returns:
    - distance (float): The distance in parsecs.
    - distance_error (float): The uncertainty in the distance.
    """
    if spectral_type not in spectral_types:
        print(f"Error: Spectral type {spectral_type} is outside the valid range.")
        return np.nan, np.nan
    
    absolute_magnitude = spectral_types[spectral_type]
    distance = 10 ** ((apparent_magnitude - absolute_magnitude + 5) / 5)
    
    # Calculate the distance uncertainty
    ln10_div_5 = math.log(10) / 5
    distance_error = distance * ln10_div_5 * magnitude_uncertainty
    
    return distance, distance_error

def process_csv(input_file, output_file):
    # Read the input CSV file
    df = pd.read_csv(input_file)
    
    # Ensure the necessary columns are present
    required_columns = ['ID_star', 'SpT', 'G_mag', 'eG_mag']
    if not all(column in df.columns for column in required_columns):
        raise ValueError(f"Input CSV must contain the columns: {', '.join(required_columns)}")
    
    # Calculate d_pc and ed_pc for each row
    distances = df.apply(lambda row: photometric_distance(row['G_mag'], row['SpT'], row['eG_mag']), axis=1)
    df['d_pc'], df['ed_pc'] = zip(*distances)
    
    # Write the results to the output CSV file
    df.to_csv(output_file, index=False)
    print(f"Output saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description='Calculate photometric distances for stars from a CSV file.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to save the output CSV file')
    
    args = parser.parse_args()
    
    try:
        process_csv(args.input, args.output)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
