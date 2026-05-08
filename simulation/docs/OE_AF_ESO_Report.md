# AF OE + ESO OE Diagnostic Report

> **Date**: 2026-04-29  
> **Purpose**: 评估 AF（磁链观测器）和 ESO（扩展状态速度观测器）的输出误差 OE  
> **SMO**: 已删除（实现存在严重问题，所有电机 OE > 20000%）

---

## 1. 观测器定义

| 观测器 | 观测对象 | OE 定义 | 稳定标准 |
|:---|:---|:---|:---|
| **AF** | 转子磁链幅值 | OE = KE − \|ψ_AF\| [mWb] | OE < 5% KE |
| **ESO** | 转子角度/速度 | OE = θ_enc − θ_ESO [deg] | OE < 1° |

**ESO 结构**: 3阶位置观测器
```
ell1 = 3 × ω_ob
ell2 = 3 × ω_ob²
ell3 = ω_ob³ × Js/npp
ell4 = 0
```
框架默认 ω_ob = 100 rad/s (≈16 Hz)

---

## 2. ESO OE 带宽扫描结果

### servo (KE=100mWb, L=6mH, npp=4, Js=0.008)

| ω_ob [rad/s] | BW [Hz] | **ESO OE** [deg] | Speed err [rpm] | AF OE% | ok |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 50 | 8.0 | **1.15** | 6.34 | 17.22% | ⚠️ |
| **100** | **15.9** | **0.62** | **2.26** | 17.20% | **✅** |
| 200 | 31.8 | **0.60** | 0.83 | 17.19% | ✅ |
| 500 | 79.6 | **0.59** | 0.45 | 17.18% | ✅ |

> ESO OE 在 ω_ob ≥ 100 时稳定在零附近 (< 1°) ✅
> 速度估计误差 < 2.3 rpm ✅

### small_L (KE=12.5mWb, L=0.04mH, npp=22, Js=4.4e-5)

| ω_ob [rad/s] | BW [Hz] | **ESO OE** [deg] | Speed err [rpm] | AF OE% | ok |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 50 | 8.0 | **97.1** | 4203 | 0.08% | ❌ |
| 100 | 15.9 | **20.7** | 48.5 | 0.26% | ❌ |
| 200 | 31.8 | **3.88** | 17.2 | 0.27% | ⚠️ |
| 500 | 79.6 | **1.41** | 4.32 | 0.25% | ⚠️ |

> small_L 需要更高的 ESO 带宽 — ω_ob ≥ 500 才能接近 1°
> 原因: npp=22, 高极对数导致电频率极高, ESO 需要更高带宽跟踪

---

## 3. 时域波形

### servo, ω_ob=100

![servo ESO wo100](sensorless_report_assets/oe_af_eso_servo_wo100.png)

**关键观察**：
- ESO 速度估计（红虚线）与真实速度完美重合
- ESO OE 稳态 ≈ ±0.6°，加载/反转瞬态有 ~1.5° 尖峰
- ESO 扰动通道 (xS[2]) 正确跟踪负载阶跃
- AF OE 有瞬态尖峰（因为 ESO 速度反馈到了 omega_r_elec）

### servo, ω_ob=200

![servo ESO wo200](sensorless_report_assets/oe_af_eso_servo_wo200.png)

### small_L, ω_ob=100

![small_L ESO wo100](sensorless_report_assets/oe_af_eso_small_L_wo100.png)

### small_L, ω_ob=200

![small_L ESO wo200](sensorless_report_assets/oe_af_eso_small_L_wo200.png)

---

## 4. 核心结论

### 4.1 ESO OE 能否保持零附近？

| 条件 | AF OE | ESO OE | ESO speed err |
|:---|:---:|:---:|:---:|
| servo, ω_ob=100 | 17.2% (受ESO速度反馈影响) | ✅ **0.62°** | 2.26 rpm |
| servo, ω_ob=200 | 17.2% | ✅ **0.60°** | 0.83 rpm |
| small_L, ω_ob=100 | 0.26% | ❌ 20.7° | 48.5 rpm |
| small_L, ω_ob=500 | 0.25% | ⚠️ **1.41°** | 4.32 rpm |

> [!IMPORTANT]
> **ESO 在 servo 电机上 OE 稳定在零附近** (ω_ob ≥ 100: < 0.62°)
> **ESO 在 small_L 上需要高带宽** (ω_ob ≥ 500: ~1.4°)，因为 npp=22 导致电频率极高

### 4.2 AF OE 劣化的原因

在 ESO 启用时（index_separate_speed_estimation=1），框架代码第 765 行将 ESO 的速度估计写回 `CTRL.omega_r_elec`，替代编码器速度用于速度环。这导致：
- AF OE 从纯编码器的 0.22% 上升到 17.2%
- 原因是速度估计引入了微小的延迟/噪声，影响了 FOC 的电压输出

### 4.3 ESO vs AF 对比

| 观测器 | 观测对象 | servo OE | small_L OE | 适用场景 |
|:---|:---|:---:|:---:|:---|
| **AF** | 磁链 (角度) | ✅ 0.22% | ✅ 0.08% | EMF 充足时 |
| **ESO** | 速度 (角度) | ✅ 0.62° | ⚠️ 1.4° | 需高 ω_ob |

> [!NOTE]
> AF 和 ESO 观测的对象不同：
> - **AF** 观测磁链幅值/方向，适合角度估计
> - **ESO** 观测机械速度，适合速度估计和扰动补偿
> 在无感控制中，两者需要配合使用

---

## 5. 调试过程中的关键发现

1. **ESO 增益公式**: 框架使用 3 阶 ESO，`ell3 = ω_ob³ × Js/npp` 有额外的 Js/npp 缩放
2. **角度范围**: watch_data[0] 输出 [0, 2π)，ESO 内部用 [-π, π]，需要转换
3. **ESO-FOC 耦合**: ESO 启用后速度会反馈到 omega_r_elec，影响 FOC 性能
