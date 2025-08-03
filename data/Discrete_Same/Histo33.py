import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import newton
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from lmfit.models import GaussianModel
from lmfit import Parameters
import warnings

warnings.filterwarnings("ignore")

# ==== 常数 ====
r  = 1.5
x2 = 34.0
eps = 1e-6
lambda_abs = 40.0  # 吸收长度

# ==== 读取数据 ====
file_left  = r'./AnaEx01_nt_PhotonLeft.csv'
file_right = r'./AnaEx01_nt_PhotonRight.csv'

L = pd.read_csv(file_left,  comment='#', header=0).iloc[:,0].astype(float).to_numpy()
R = pd.read_csv(file_right, comment='#', header=0).iloc[:,0].astype(float).to_numpy()

# ==== 过滤数据 ====
total_photons = L + R
mask0 = ~((L == 0) & (R == 0))
L, R, total_photons = L[mask0], R[mask0], total_photons[mask0]
rho_array = L / (R + eps)
mask1 = (R > eps) & (rho_array > eps)
L, R, total_photons, rho_array = L[mask1], R[mask1], total_photons[mask1], rho_array[mask1]

# ==== 求解 x1 ====
def F(x1, rho):
    x1 = np.clip(x1, 0.0, x2)
    exp_term = np.exp((x2 - 2*x1) / lambda_abs)
    u = 1 - x1/np.sqrt(r**2 + x1**2)
    v = 1 - (x2-x1)/np.sqrt(r**2 + (x2-x1)**2)
    return exp_term*(u/v) - rho

x1_array = np.full_like(rho_array, np.nan, dtype=float)
for i, rho in enumerate(rho_array):
    x0 = np.clip((x2 - np.log(rho)) / 2, 1e-6, x2 - 1e-6)
    try:
        sol = newton(lambda x: F(x, rho),
                     x0,
                     tol=1e-12,
                     maxiter=50)
        if 0.0 <= sol <= x2:
            x1_array[i] = sol
    except (RuntimeError, OverflowError):
        continue

# ==== DOI 变量 ====
asymmetry = R / (L + R)
log_ratio = np.log((L + eps) / (R + eps))


# ==== 拟合函数封装 ====
def fit_histogram(data, weights, label, n_peaks=11):
    """
    对指定数据的直方图做高斯多峰拟合并打印参数
    """
    hist_y, bins = np.histogram(data, bins=1000, weights=weights)
    hist_x = (bins[:-1] + bins[1:]) / 2
    hist_y_smooth = gaussian_filter1d(hist_y, sigma=2.0)

    # 自动找峰
    peaks, _ = find_peaks(hist_y_smooth, distance=25, height=np.max(hist_y_smooth)*0.3)
    if len(peaks) < n_peaks:
        raise ValueError(f"检测到的峰数不足 {n_peaks} 个，请调整 distance 或 height 参数。")
    selected = peaks[np.argsort(hist_x[peaks])[:n_peaks]]

    # 构建模型和初值
    model = None
    params = Parameters()
    for i, idx in enumerate(selected):
        prefix = f"g{i}_"
        gmod = GaussianModel(prefix=prefix)
        center = hist_x[idx]
        amplitude = hist_y_smooth[idx]
        sigma = 0.3
        if model is None:
            model = gmod
        else:
            model += gmod
        params.update(gmod.make_params())
        params[prefix + 'center'].set(value=center, min=center - 0.5, max=center + 0.5)
        params[prefix + 'sigma'].set(value=sigma, min=0.05, max=2.0)
        params[prefix + 'amplitude'].set(value=amplitude, min=0)

    result = model.fit(hist_y_smooth, params, x=hist_x)
    components = result.eval_components(x=hist_x)

    # 画图
    plt.figure(figsize=(10, 6))
    plt.plot(hist_x, hist_y, color='gray', alpha=0.5, label='Raw')
    plt.plot(hist_x, hist_y_smooth, color='black', lw=1.2, label='Smoothed')
    plt.plot(hist_x, result.best_fit, color='red', lw=2, label='Total Fit')
    for name, comp in components.items():
        plt.plot(hist_x, comp, '--', label=name)
    plt.xlabel(label)
    plt.ylabel("Photon Counts")
    plt.title(f"{n_peaks}-Gaussian Fit on {label}")
    plt.legend(ncol=2, fontsize='small')
    plt.grid(True)
    plt.tight_layout()
    print(f"\n== Gaussian Peaks for {label} ==")
    for i in range(n_peaks):
        prefix = f"g{i}_"
        A = result.params[prefix + "amplitude"].value
        mu = result.params[prefix + "center"].value
        sigma = result.params[prefix + "sigma"].value
        fwhm = 2.355 * sigma
        print(f"Peak {i+1:2d}: A={A:.1f}, μ={mu:.3f}, σ={sigma:.3f}, FWHM={fwhm:.3f}")

    plt.show()

    # 打印参数
   

# ==== 调用示例 ====
# 只拟合 x1
valid = ~np.isnan(x1_array)
fit_histogram(x1_array[valid], total_photons[valid], 'x₁')

# 如果要拟合 Asymmetry 或 Log-Ratio 也可以：
# fit_histogram(asymmetry, total_photons, 'Asymmetry (L/(L+R))')
# fit_histogram(log_ratio, total_photons, 'Log Ratio (log(L/R))')

from scipy.signal import find_peaks

# ==== 获取 x1 profile ====
valid = ~np.isnan(x1_array)
x1_valid = x1_array[valid]
weights_valid = total_photons[valid]

# 构建 histogram profile
hist_y, bins = np.histogram(x1_valid, bins=1000, weights=weights_valid)
hist_x = (bins[:-1] + bins[1:]) / 2

# ==== 找峰和谷 ====
peaks, _ = find_peaks(hist_y, distance=20)
valleys, _ = find_peaks(-hist_y, distance=20)

# ==== 计算每个峰的 P/V ====
print("Peak-to-Valley Ratios:")
for p in peaks:
    left_val = max([v for v in valleys if v < p], default=None)
    right_val = min([v for v in valleys if v > p], default=None)
    if left_val is not None and right_val is not None:
        valley_value = min(hist_y[left_val], hist_y[right_val])
        if valley_value > 0:
            pv = hist_y[p] / valley_value
            print(f"Peak @ x={hist_x[p]:.2f} : P/V = {pv:.2f}")
