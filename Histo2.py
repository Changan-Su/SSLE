import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ==== 文件名 ====
file_left = r'.\build\AnaEx01_nt_PhotonLeft.csv'
file_right = r'.\build\AnaEx01_nt_PhotonRight.csv'

# 跳过注释行并读取第一列数据
L = pd.read_csv(file_left, comment='#', header=0).iloc[:, 0].astype(float).to_numpy()
R = pd.read_csv(file_right, comment='#', header=0).iloc[:, 0].astype(float).to_numpy()

# ==== 基本处理 ====
eps = 1e-6  # 防止除零
total_photons = L + R

# ==== 两种 DOI 计算方式 ====
asymmetry = (L - R) / (L + R + eps)
log_ratio = np.log((L + eps) / (R + eps))
doi2 = L / (L+R+eps) 

# ==== Asymmetry 直方图 ====
plt.figure(figsize=(8, 5))
plt.hist(asymmetry, weights=total_photons, bins=1000, color='blue', alpha=0.7)
plt.xlabel("DOI (Asymmetry: (L - R)/(L + R))")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs DOI (Asymmetry)")
plt.grid(True)
plt.tight_layout()

# ==== Log-Ratio 直方图 ====
plt.figure(figsize=(8, 5))
plt.hist(log_ratio, weights=total_photons, bins=1000, color='green', alpha=0.7)
plt.xlabel("DOI (Log Ratio: log(L / R))")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs DOI (Log Ratio)")
plt.grid(True)
plt.tight_layout()

#doi 2
plt.figure(figsize=(8, 5))
plt.hist(doi2, weights=total_photons, bins=1000, color='red', alpha=0.7)
plt.xlabel("DOI : L/(L+R)")
plt.ylabel("Total Photon Counts")
plt.title("Photon Counts vs DOI (L/L+R)")
plt.grid(True)
plt.tight_layout()

# ==== 显示图像 ====
plt.show()