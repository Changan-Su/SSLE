import pandas as pd
import uproot

# 读取数据
df_left = pd.read_csv("AnaEx01_nt_PhotonLeft.csv", comment='#')
df_right = pd.read_csv("AnaEx01_nt_PhotonRight.csv", comment='#')

# 合并为一个 DataFrame
df_combined = pd.DataFrame({
    "L": df_left.iloc[:, 0],
    "R": df_right.iloc[:, 0],
})

# 写入 ROOT 文件
with uproot.recreate("PhotonData.root") as f:
    f["PhotonTree"] = df_combined
