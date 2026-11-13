import glob
import yaml
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata

config = yaml.safe_load(open('test/config.yaml', 'r'))

# Find celltypist output files
files = glob.glob('test/out/celltype_prediction/dataset~*/file_id~*.zarr')
assert len(files) > 0, 'No celltypist output files found'

for file in files:
    matching_datasets = [x for x in config['DATASETS'] if f'dataset~{x}/' in file]
    if not matching_datasets:
        logging.info(f'No matching dataset found for {file}, skipping...')
        continue
    
    dataset = matching_datasets[0]
    logging.info(f'Check celltypist dataset "{dataset}", file {file}')
    
    adata_out = read_anndata(file, verbose=False)
    
    try:
        # Check celltypist columns exist (should start with celltypist_)
        celltypist_cols = [col for col in adata_out.obs.columns if col.startswith('celltypist_')]
        assert len(celltypist_cols) > 0, "No celltypist columns found (should start with 'celltypist_')"
        
        # Check for predicted_labels column (at minimum)
        expected_cols = [col for col in celltypist_cols if 'predicted_labels' in col]
        assert len(expected_cols) > 0, "No predicted_labels column found"
        
        logging.info(f'Found celltypist columns: {celltypist_cols}')
        logging.info(f"✅ All celltypist assertions passed for dataset={dataset}")
        
    except AssertionError as e:
        logging.error(f"❌ Assertion failed for {dataset}: {e}")
        logging.error(adata_out.obs.columns.tolist())
        raise e

for file in files:
    matching_datasets = [x for x in config['DATASETS'] if f'dataset~{x}/' in file]
    if not matching_datasets:
        logging.info(f'No matching dataset found for {file}, skipping...')
        continue
    
    dataset = matching_datasets[0]
    logging.info(f'Check predict_sex dataset "{dataset}", file {file}')
    
    adata_out = read_anndata(file, verbose=False)
    predict_config = config['DATASETS'][dataset]['celltype_prediction'].get('predict_sex')
    
    try:
        if predict_config:
            predict_column = predict_config.get('predict_column', 'sex')
            
            # Check prediction column exists
            assert predict_column in adata_out.obs.columns, f"Prediction column '{predict_column}' not found"
            
            # Check it has valid values
            assert adata_out.obs[predict_column].notna().any(), f"Prediction column '{predict_column}' is all NaN"
            
            logging.info(f'Found prediction column: {predict_column}')
            logging.info(f'Value counts:\n{adata_out.obs[predict_column].value_counts()}')
            logging.info(f"✅ All predict_sex assertions passed for dataset={dataset}")
    
    except AssertionError as e:
        logging.error(f"❌ Assertion failed for {dataset}: {e}")
        logging.error(adata_out.obs.columns.tolist())
        raise e

logging.info("🎉 All celltype_prediction assertion tests passed!")
