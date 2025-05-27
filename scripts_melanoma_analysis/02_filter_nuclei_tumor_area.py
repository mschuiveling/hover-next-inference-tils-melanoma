import os
import json
import argparse
from pathlib import Path
import pandas as pd
import geopandas as gpd
from tqdm import tqdm
from cellseg_gsontools import read_gdf
from cellseg_gsontools.grid import grid_classify, grid_overlay

def get_total_cell_cnt_per_class(gdf: gpd.GeoDataFrame, class_name: str, **kwargs) -> int:
    if gdf is None:
        return 0
    return len(gdf[gdf['class_name'] == class_name])

def add_percentage_columns(grid, nuclei):
    if nuclei is None or nuclei.empty:
        print("Nuclei data is None or empty. Skipping percentage calculation.")
        return None

    try:
        unique_class_names = nuclei['class_name'].unique()
    except KeyError:
        print("Key 'class_name' not found in nuclei data. Skipping this file.")
        return None
    
    for class_name in unique_class_names:
        grid = grid_classify(
            grid=grid,
            objs=nuclei,
            metric_func=lambda gdf, **kwargs: get_total_cell_cnt_per_class(gdf, class_name, **kwargs),
            predicate="intersects",
            new_col_names=f"total_cnt_{class_name}",
            parallel=False,
        )
    
    if grid is None or grid.empty:
        print("Grid classification resulted in an empty grid. Skipping further processing.")
        return None

    try:
        grid['total_cell_cnt'] = grid[[f"total_cnt_{class_name}" for class_name in unique_class_names]].sum(axis=1)
        for class_name in unique_class_names:
            percentage_col = f"{class_name}_percentage"
            grid[percentage_col] = (grid[f"total_cnt_{class_name}"] / grid['total_cell_cnt']) * 100
            grid[percentage_col].fillna(0, inplace=True)
    except KeyError as e:
        print(f"Error during percentage calculation: {e}")
        return None

    return grid

def process_geojson_files(nuclei_files, tissue_files, nuclei_geojson_folder, tissue_geojson_folder, results_path, geojson_output_path, pixel_grid_size=880.6693086745927):
    for geojson in tqdm(nuclei_files):
        if geojson in tissue_files:
            output_filepath = Path(results_path) / f"{Path(geojson).stem}.csv"
            lock_file = Path(results_path) / f"{Path(geojson).stem}.lock"

            if os.path.exists(output_filepath):
                print(f"Skipping {geojson} as it already exists")
                continue

            if os.path.exists(lock_file):
                print(f"Skipping {geojson} as it is locked")
                continue
            else:
                open(lock_file, 'w').close()

            tissue_geojson_path = Path(tissue_geojson_folder) / geojson
            if os.path.getsize(tissue_geojson_path) < 2 * 1024:
                print(f"Skipping {geojson} as tissue file is smaller than 2KB")
                continue

            print(f'Processing {geojson}')

            try:
                tissue = read_gdf(Path(tissue_geojson_folder) / geojson)
            except Exception as e:
                print(f"Skipping {geojson} due to error loading data: {e}")
                continue

            try:
                nuclei = read_gdf(Path(nuclei_geojson_folder) / geojson)
            except Exception as e:
                print(f"Skipping {geojson} due to error loading data: {e}")
                continue
            
            try:
                tissue = tissue[tissue['classification'].notnull()]
                def parse_classification(x):
                    if isinstance(x, str):
                        try:
                            return json.loads(x)
                        except json.JSONDecodeError:
                            print(f"Warning: Failed to parse classification: {x}")
                            return {"name": None}
                    elif isinstance(x, dict):
                        return x
                    else:
                        return {"name": None}
                tissue['classification'] = tissue['classification'].apply(parse_classification)
                tissue['class_name'] = tissue['classification'].apply(lambda x: x.get("name", "").lower() if x else None)
                nuclei['classification'] = nuclei['classification'].apply(parse_classification)
                nuclei['class_name'] = nuclei['classification'].apply(lambda x: x.get("name", ""))
            except KeyError as e:
                print(f"Skipping {geojson} due to missing key in classification: {e}")
                continue

            tumor_filter = tissue[tissue["class_name"] == "tumor"]
            necrosis_filter = tissue[tissue["class_name"] == "necrosis"]

            try:
                nuclei_within_tumor = gpd.sjoin(nuclei, tumor_filter, how='inner', predicate='within')
                nuclei_within_tumor = nuclei_within_tumor.drop(columns=['index_right', 'id_right', 'objectType_right', 'classification_right', 'class_name_right'])
                nuclei_within_tumor.columns = [col.replace('_left', '') for col in nuclei_within_tumor.columns]
                nuclei_within_necrosis = gpd.sjoin(nuclei_within_tumor, necrosis_filter, how='inner', predicate='within')
                nuclei_within_necrosis_ids = set(nuclei_within_necrosis.index)
                nuclei_within_tumor_not_necrosis = nuclei_within_tumor[~nuclei_within_tumor.index.isin(nuclei_within_necrosis_ids)]
                nuclei_within_tumor_not_necrosis = nuclei_within_tumor_not_necrosis.drop_duplicates(subset=['geometry'])
            except Exception as e:
                print(f"Error during spatial join or data cleanup: {e}")
                continue

            nuclei_within_tissue_output_path = Path(geojson_output_path) / f"{Path(geojson).stem}_nuclei_within_tumor_not_necrosis.geojson"
            nuclei_within_tumor_not_necrosis.to_file(nuclei_within_tissue_output_path, driver='GeoJSON')

            try:
                grid = grid_overlay(nuclei_within_tumor_not_necrosis, patch_size=(pixel_grid_size, pixel_grid_size), stride=(pixel_grid_size, pixel_grid_size))
                grid = add_percentage_columns(grid, nuclei_within_tumor_not_necrosis)
            except Exception as e:
                print(f"Error during grid overlay or percentage calculation: {e}")
                continue

            if grid is not None:
                grid.to_csv(output_filepath, index=False)

            os.remove(lock_file)

    print('No more files to analyze; done')

def main():
    parser = argparse.ArgumentParser(description="Process nuclei and tissue GeoJSON files.")
    parser.add_argument("--nuclei_geojson_folder", type=str, required=True, help="Folder containing nuclei GeoJSON files.")
    parser.add_argument("--tissue_geojson_folder", type=str, required=True, help="Folder containing tissue GeoJSON files.")
    parser.add_argument("--results_path", type=str, required=True, help="Folder to store results CSVs.")
    parser.add_argument("--geojson_output_path", type=str, required=True, help="Folder to save output GeoJSONs.")
    parser.add_argument("--pixel_grid_size", type=float, default=880.6693086745927, help="Grid size in pixels.") # Default is 880.6693086745927 --> 200 microns at 0.227 microns per pixel

    args = parser.parse_args()

    if not os.path.exists(args.results_path):
        os.makedirs(args.results_path)
    if not os.path.exists(args.geojson_output_path):
        os.makedirs(args.geojson_output_path)

    nuclei_files = os.listdir(args.nuclei_geojson_folder)
    tissue_files = os.listdir(args.tissue_geojson_folder)

    process_geojson_files(
        nuclei_files=nuclei_files,
        tissue_files=tissue_files,
        nuclei_geojson_folder=args.nuclei_geojson_folder,
        tissue_geojson_folder=args.tissue_geojson_folder,
        results_path=args.results_path,
        geojson_output_path=args.geojson_output_path,
        pixel_grid_size=args.pixel_grid_size,
    )

if __name__ == "__main__":
    main()
