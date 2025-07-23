import pandas as pd

# 文件路径
file_left_path = "AnaEx01_nt_PhotonLeft.csv"
file_right_path = "AnaEx01_nt_PhotonRight.csv"

# 读取 CSV 文件（保留表头）
df_left = pd.read_csv(file_left_path, comment='#', header=0)
df_right = pd.read_csv(file_right_path, comment='#', header=0)

# 获取数据列名
col_left = df_left.columns[0]
col_right = df_right.columns[0]

# 创建掩码：保留不是同时为0的行
mask = ~((df_left[col_left] == 0) & (df_right[col_right] == 0))

# 过滤数据
df_left_filtered = df_left[mask]
df_right_filtered = df_right[mask]

# 保存结果到新的 CSV 文件（保留表头，不保存索引）
df_left_filtered.to_csv("AnaEx01_nt_PhotonLeft_filtered.csv", index=False)
df_right_filtered.to_csv("AnaEx01_nt_PhotonRight_filtered.csv", index=False)
