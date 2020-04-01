from pathlib import Path
import numpy as np
import pandas as pd


data_dir_raw=Path("../data")


brk_rh= pd.read_csv(data_dir_raw/"brk_rh.csv")["x"]
col_rh= pd.read_csv(data_dir_raw/"col_rh.csv")["x"]


brk_tem= pd.read_csv(data_dir_raw/"brk_tem.csv")["x"]
col_tem= pd.read_csv(data_dir_raw/"col_tem.csv")["x"]

brk_ws= pd.read_csv(data_dir_raw/"brk_ws.csv")["x"]
col_ws= pd.read_csv(data_dir_raw/"col_ws.csv")["x"]






