import yaml


config_file = "../nguyenb_config.yaml"
with open(config_file) as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    configs = yaml.load(file, Loader=yaml.FullLoader)
    
# Run on server:
root = Path(configs['root']['server'])
scratchDir = Path(configs['scratchDir']['server'])
mapDir = root/configs['mapDir']
countDir = root/configs['libraryCountsDir']
resultDir = root/configs['resultDir']
sample_data_file = root/configs['sampleData']
