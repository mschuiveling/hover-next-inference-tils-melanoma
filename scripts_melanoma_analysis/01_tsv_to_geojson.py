import os
import pandas as pd
import json
import uuid
import argparse
from tqdm import tqdm

# Define the colors for each cell type
COLORS = {
    "nuclei_lymphocyte": [0, 255, 255],  # Cyan
    "nuclei_tumor": [255, 0, 0],         # Red
    "nuclei_other": [255, 255, 0],        # Yellow
}

# Function to create a polygon from a center point
def create_square_polygon(x, y, size=15):
    half_size = size / 2
    return [
        [x - half_size, y - half_size],
        [x + half_size, y - half_size],
        [x + half_size, y + half_size],
        [x - half_size, y + half_size],
        [x - half_size, y - half_size]
    ]

# Function to process TSV files in a subfolder and create a GeoJSON file
def process_subfolder(subfolder_path, output_file):
    geojson = {
        "type": "FeatureCollection",
        "features": []
    }
    
    tsv_files = [f for f in os.listdir(subfolder_path) if f.endswith('.tsv')]
    
    # Skip processing if no TSV files are found
    if not tsv_files:
        print(f"No TSV files found in {subfolder_path}. Skipping...")
        return
    
    # Process each TSV file in the subfolder
    for file_name in tqdm(tsv_files, desc=f'Processing {subfolder_path}', unit='file'):
        file_path = os.path.join(subfolder_path, file_name)
        data = pd.read_csv(file_path, sep='\t')
        
        for _, row in data.iterrows():
            x, y = row['x'], row['y']
            polygon = create_square_polygon(x, y)
            feature = {
                "type": "Feature",
                "id": str(uuid.uuid4()),
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [polygon]
                },
                "properties": {
                    "objectType": "annotation",
                    "classification": {
                        "name": row['name'],
                        "color": COLORS.get(row['name'], [0, 0, 0])
                    }
                }
            }
            geojson['features'].append(feature)
    
    # Save the GeoJSON to a file in the same folder as the TSV files
    # Save one folder above the subfolder_path
    output_file_path = os.path.join(os.path.dirname(subfolder_path), output_file)
    
    with open(output_file_path, 'w') as f:
        json.dump(geojson, f, indent=2)

    print(f"GeoJSON file created: {output_file_path}")

def main(input_dir):
    subfolders = [d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d))]

    subfolders = sorted(subfolders, reverse=True)
    
    # Process each subfolder in the main directory
    for subfolder_name in tqdm(subfolders, desc='Processing subfolders', unit='folder'):
        subfolder_path = os.path.join(input_dir, subfolder_name)
        output_file = f'{subfolder_name}.geojson'
        output_file_path = os.path.join(os.path.dirname(subfolder_path), output_file)
        if os.path.exists(output_file_path):
            print(f"Output file {output_file} already exists. Skipping...")
            continue
        process_subfolder(subfolder_path, output_file)

    print('Done!')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process TSV files into GeoJSON.')
    parser.add_argument('--input_dir', type=str, required=True, help='Path to the main directory containing subfolders with TSV files.')

    args = parser.parse_args()
    main(args.input_dir)
