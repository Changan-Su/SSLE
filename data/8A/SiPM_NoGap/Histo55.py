# -*- coding: utf-8 -*-
"""
多峰高斯拟合（支持弱峰）：prominence找峰 + 受限拟合 + FWHM与P/V输出
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks

# ===================== 参数区（按需微调） =====================
file_left  = './AnaEx01_nt_PhotonLeft.csv'
file_right = './AnaEx01_nt_PhotonRight.csv'

n_bins = 1000                 # 直方图bin数
sigma_smooth = 1.0            # 平滑强度（过大容易并峰；0.8~1.2常用）
distance_ratio = 0.035        # 峰间最小距离占比（相对bin数）
prominence_ratio = 0.01       # 初始突显度阈值占比（相对全局最大值）
target_n = 11                 # 期望的峰数；不想限定就设为None

sigma_init_guess = 0.008      # 每个峰初始sigma
mu_window = 0.015             # mu拟合允许偏移窗口（左右各±）
sigma_bounds = (0.003, 0.030) # sigma的约束范围
# ============================================================

# ==== 数据读取 ====
L = pd.read_csv(file_left,  comment='#', header=0).iloc[:, 0].astype(float).to_numpy()
R = pd.read_csv(file_right, comment='#', header=0).iloc[:, 0].astype(float).to_numpy()

# ==== DOI 计算 ====
eps = 1e-6
total_photons = L + R
asymmetry = R / (L + R + eps)

# 可选：去除非数与无穷
mask = np.isfinite(asymmetry) & np.isfinite(total_photons)
asymmetry = asymmetry[mask]
total_photons = total_photons[mask]

# ==== 构建直方图 ====
hist_y, bins = np.histogram(asymmetry, weights=total_photons, bins=n_bins, range=(0.0, 1.0))
hist_x = (bins[:-1] + bins[1:]) / 2.0
hist_y_smooth = gaussian_filter1d(hist_y.astype(float), sigma=sigma_smooth)

# ==== 自动寻找峰：使用 prominence 而不是 height ====
peak_distance = max(1, int(distance_ratio * len(hist_x)))
prom0 = prominence_ratio * float(np.max(hist_y_smooth))

def detect_peaks(prom):
    return find_peaks(
        hist_y_smooth,
        distance=peak_distance,
        prominence=(prom, None),     # 关键
        width=(2, None),             # 过滤尖刺
        wlen=int(0.15 * len(hist_x)) # 计算prominence的窗口
    )

# 自适应放宽阈值以捕获弱峰
prom_candidates = [prom0, 0.008*np.max(hist_y_smooth), 0.006*np.max(hist_y_smooth),
                   0.004*np.max(hist_y_smooth), 0.003*np.max(hist_y_smooth),
                   0.002*np.max(hist_y_smooth), 0.001*np.max(hist_y_smooth)]

peaks, props = np.array([], dtype=int), {}
for prom in prom_candidates:
    peaks, props = detect_peaks(prom)
    if target_n is None:
        if len(peaks) >= 3:  # 至少抓到几个明显峰就停
            break
    else:
        if len(peaks) >= target_n:
            # 如果抓到的比期望多，按prominence从大到小截取
            if len(peaks) > target_n:
                order = np.argsort(props["prominences"])[::-1][:target_n]
                peaks = peaks[order]
            break

# 若仍没抓够，最后一次使用最低prom
if len(peaks) == 0:
    peaks, props = detect_peaks(prom_candidates[-1])

# 按x从小到大排序
order = np.argsort(peaks)
peaks = peaks[order]

# ==== 多峰高斯模型 ====
def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        A, mu, sigma = params[i], params[i+1], params[i+2]
        y += A * np.exp(-(x - mu) ** 2 / (2.0 * sigma ** 2))
    return y

# ==== 初始猜测与边界 ====
init_params = []
lb, ub = [], []
for pk in peaks:
    A_init = max(hist_y_smooth[pk], 1.0)
    mu_init = hist_x[pk]
    s_init = sigma_init_guess
    init_params += [A_init, mu_init, s_init]

    lb += [0.0,            mu_init - mu_window, sigma_bounds[0]]
    ub += [np.inf,         mu_init + mu_window, sigma_bounds[1]]

# ==== 拟合 ====
try:
    popt, pcov = curve_fit(
        multi_gaussian, hist_x, hist_y_smooth,
        p0=init_params, bounds=(lb, ub),
        maxfev=20000
    )
except Exception as e:
    # 兜底：放宽sigma上界并再试一次
    ub2 = ub[:]
    for i in range(2, len(ub2), 3):
        ub2[i] = max(0.05, ub2[i])
    popt, pcov = curve_fit(
        multi_gaussian, hist_x, hist_y_smooth,
        p0=init_params, bounds=(lb, ub2),
        maxfev=40000
    )

fitted_y = multi_gaussian(hist_x, *popt)

# ==== 提取参数并计算 FWHM ====
peaks_params = []
for i in range(0, len(popt), 3):
    A = popt[i]
    mu = popt[i+1]
    sigma = popt[i+2]
    fwhm = 2.354820045 * sigma
    peaks_params.append((mu, A, sigma, fwhm))

# 按mu从小到大排序，保证输出顺序一致
peaks_params.sort(key=lambda t: t[0])

# ==== 峰谷比（P/V） ====
# 在平滑曲线上找“谷”
valleys, _ = find_peaks(-hist_y_smooth, distance=max(1, int(0.02*len(hist_x))))
pv_results = []
# 用拟合的mu定位到最近的离散索引
mu_indices = [np.clip(np.searchsorted(hist_x, mu)-1, 0, len(hist_x)-1) for (mu, *_ ) in peaks_params]

for idx in mu_indices:
    # 最近的左/右谷
    left_candidates = valleys[valleys < idx]
    right_candidates = valleys[valleys > idx]
    if len(left_candidates) == 0 or len(right_candidates) == 0:
        continue
    left_v = left_candidates[-1]
    right_v = right_candidates[0]
    valley_min = float(min(hist_y_smooth[left_v], hist_y_smooth[right_v]))
    peak_height = float(hist_y_smooth[idx])
    if valley_min > 0:
        pv = peak_height / valley_min
        pv_results.append((hist_x[idx], pv))

# ==== 绘图 ====
plt.figure(figsize=(16, 6))
plt.plot(hist_x, hist_y,        label="Raw Histogram", alpha=0.35)
plt.plot(hist_x, hist_y_smooth, label="Smoothed",      linewidth=1.2)
plt.plot(hist_x, fitted_y,      label="Total Fit",     linewidth=2.0)
plt.scatter(hist_x[peaks], hist_y_smooth[peaks], s=25, zorder=5, label="Detected Peaks")

plt.xlabel("DOI (Asymmetry: R / (L + R))")
plt.ylabel("Total Photon Counts")
plt.title("Multi-Gaussian Fit with FWHM and P/V Ratios (prominence-based)")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# ==== 打印结果 ====
print("\n== Gaussian Peak Fit Results ==")
for i, (mu, A, sigma, fwhm) in enumerate(peaks_params, 1):
    print(f"Peak {i:2d}: μ = {mu:.6f}, A = {A:.1f}, σ = {sigma:.6f}, FWHM = {fwhm:.6f}")

print("\n== Peak-to-Valley Ratios ==")
for x, pv in pv_results:
    print(f"Peak near x = {x:.6f} → P/V = {pv:.2f}")

# ==== 提示：如果中间峰仍缺失，可适当再做三件事 ====
# 1) sigma_smooth ↓ 到 0.8 或 0.6
# 2) prominence_ratio ↓ 到 0.005 或更低
# 3) 增大 target_n 或放宽 mu_window、sigma_bounds 上界
