import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pathlib import Path
import seaborn as sns
import sys
import plotly.express as px
import plotly.io as pio
import yaml

sns.set_context("notebook", font_scale=1.4)
pd.set_option("display.max_columns", 100)
pd.set_option("display.max_rows", 100)
plt.rcParams["figure.figsize"] = (16, 12)
plt.rcParams['savefig.dpi'] = 200
plt.rcParams['figure.autolayout'] = False
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 20
plt.rcParams['font.size'] = 16
plt.rcParams['lines.linewidth'] = 2.0
plt.rcParams['lines.markersize'] = 8
plt.rcParams['legend.fontsize'] = 14
pd.set_option('display.float_format', lambda x: '{:,.4f}'.format(x))


config_file = "manuscript_config.yaml"
with open(config_file) as file:
    # The FullLoader parameter handles the conversion from YAML
    # scalar values to Python the dictionary format
    configs = yaml.load(file, Loader=yaml.FullLoader)
    
# Run on server:
run_on = "server"
root = Path(configs['root'][run_on])
scratchDir = Path(configs['scratchDir'][run_on])
figuresDir = Path(configs['figuresDir'][run_on])

