# Open-Loop Observer OE (Output Error) Analysis Report

> **Date**: 2026-04-29  
> **Purpose**: 在开环（编码器控制）条件下，诊断 AF/PLL/SMO 三种观测器的 OE 时域表现  
> **核心问题**: OE = KE − |ψ_obs| 能否在稳态下保持零附近波动？

---

## 1. 实验方法

**Stage 0 开环观测**：速度和角度均由编码器真值控制（FOC 闭环），AF/PLL/SMO 仅作为**纯观测器**运行，不反馈到控制器。

**速度剖面**：0→cmd_rpm (斜坡0.3s) → 稳态 → 加载 → 卸载 → 反转 → 稳态

**参数条件**：
- **Ideal**: 电阻无失配 (R_mis=1.0)，无电压偏移
- **Mismatch**: R×1.5, 电压偏移 (servo: 0.5V/-0.3V, big_L: 1.0V/-0.8V)

---

## 2. 三台电机参数

| 参数 | servo | small_L | big_L |
|:---|:---:|:---:|:---:|
| KE [mWb] | 100.0 | 12.5 | 74.5 |
| L [mH] | 6.00 | 0.04 | 106.3 |
| R [Ω] | 1.1 | 0.035 | 1.97 |
| **KE/L** | **16.7** | **347** | **0.7** |
| npp | 4 | 21 | 5 |
| cmd_rpm | 500 | 200 | 150 |

---

## 3. OE 稳态统计 (t > 0.4s)

### 3.1 Ideal Parameters (No Mismatch)

| Motor | Observer | OE_RMS [mWb] | OE_Max [mWb] | OE/%KE | Angle RMS | OE≈0? |
|:---|:---|:---:|:---:|:---:|:---:|:---:|
| **servo** | AF | 0.222 | 2.377 | **0.2%** | **0.93°** | ✅ |
| | PLL | 0.414 | 5.301 | **0.4%** | 2.01° | ✅ |
| | SMO | 20406 | 25855 | 20406% | 107° | ❌ |
| **small_L** | AF | 0.010 | 0.072 | **0.1%** | **1.55°** | ✅ |
| | PLL | 0.010 | 0.072 | **0.1%** | 4.35° | ✅ |
| | SMO | 3531 | 4411 | 28246% | 101° | ❌ |
| **big_L** | AF | 12.2 | 63.6 | **16.4%** | 93.3° | ⚠️ |
| | PLL | 27.2 | 67.7 | **36.5%** | 67.7° | ❌ |
| | SMO | 14627 | 27161 | 19634% | 91.4° | ❌ |

### 3.2 With Parameter Mismatch (R×1.5 + Voltage Offset)

| Motor | Observer | OE_RMS [mWb] | OE/%KE | Angle RMS | OE≈0? |
|:---|:---|:---:|:---:|:---:|:---:|
| **servo** | AF | 7.990 | **8.0%** | 15.6° | ⚠️ |
| | PLL | 8.181 | 8.2% | 15.8° | ⚠️ |
| | SMO | 20431 | 20431% | 108° | ❌ |
| **small_L** | AF | 0.080 | **0.6%** | 1.34° | ✅ |
| | PLL | 0.080 | 0.6% | 3.95° | ✅ |
| | SMO | 3477 | 27817% | 102° | ❌ |
| **big_L** | AF | 105 | **141%** | 115° | ❌ |
| | PLL | 116 | 156% | 115° | ❌ |
| | SMO | 25475 | 34195% | 128° | ❌ |

---

## 4. OE 时域波形

### 4.1 servo — Ideal

![servo OE ideal](sensorless_report_assets/oe_timeseries_servo_ideal.png)

**关键观察**：
- **AF OE (红)**: 完全在零附近，RMS 仅 0.222 mWb = 0.2% KE ✅
- **PLL OE (绿)**: 与 AF 几乎重合，RMS 0.414 mWb = 0.4% KE ✅
- **SMO OE (紫)**: 完全不在零附近，偏移约 -20000 mWb ❌
- AF 修正量远小于 EMF 信号（底部子图），说明观测器运行在"轻载"模式

### 4.2 small_L — Ideal

![small_L OE ideal](sensorless_report_assets/oe_timeseries_small_L_ideal.png)

**关键观察**：
- AF/PLL OE 均极小 (0.010 mWb = 0.1% KE)，**OE 完美在零** ✅
- 角度误差也极小 (AF 1.55°, PLL 4.35°)
- KE/L=347 → EMF 远大于电感压降 → 积分器收敛极快

### 4.3 big_L — Ideal

![big_L OE ideal](sensorless_report_assets/oe_timeseries_big_L_ideal.png)

**关键观察**：
- AF OE 在稳态约 ±10 mWb 波动，**但角度误差 93°！** ⚠️
- **OE 小不等于角度准**：OE 只反映磁链幅值误差，不反映方向误差
- big_L 的 OE 在加载/反转瞬态时有较大波动 (±64 mWb)
- KE/L=0.7 → 电感压降淹没 EMF → 积分器无法收敛到正确方向

### 4.4 servo — Parameter Mismatch

![servo OE mismatch](sensorless_report_assets/oe_timeseries_servo_mismatch.png)

**关键观察**：
- R 失配导致 OE 从 0.2% 上升到 8.0% KE
- 角度误差从 0.93° 劣化到 15.6°
- AF 修正量增大（因为 R_obs 错误导致电压方程残差增大）

---

## 5. AF 增益对 OE 的影响

![Gain Sweep](sensorless_report_assets/fig_gain_sweep_3motors.png)

**三台电机的增益扫描 (Kp 从 0 到 10000)**：

### servo (KE/L=16.7)
- **Kp=0 (纯积分)**: OE = 4.0 mWb, 角度 2.6° — 纯积分即可工作
- **Kp=50~100**: OE 降到 ~0.2 mWb, 角度最优 0.75~0.8°
- **Kp>500**: OE 被压到 ~0, 但角度反而劣化到 0.85° — **过高增益的修正方向误差**

### small_L (KE/L=347)
- OE 始终极小 (<0.1 mWb)，增益几乎无影响
- 角度随增益增大而减小（1.5°→1.5°），EMF 占主导

### big_L (KE/L=0.7)
- **无论增益如何，角度始终 91~96°** — 纯积分方向就错了
- OE 可以被高增益压小，但**方向不对 = OE小但角度大**

---

## 6. 核心结论

### 6.1 OE 能否保持零附近？

| 条件 | AF | PLL | SMO |
|:---|:---:|:---:|:---:|
| servo ideal | ✅ **0.2%** | ✅ 0.4% | ❌ |
| servo mismatch | ⚠️ 8.0% | ⚠️ 8.2% | ❌ |
| small_L ideal | ✅ **0.1%** | ✅ 0.1% | ❌ |
| small_L mismatch | ✅ 0.6% | ✅ 0.6% | ❌ |
| big_L ideal | ⚠️ 16.4% | ❌ 36.5% | ❌ |
| big_L mismatch | ❌ 141% | ❌ 156% | ❌ |

> [!IMPORTANT]
> **回答核心问题**：
> - **AF 和 PLL 在 servo/small_L 上 OE 稳定在零附近** (理想参数 <1% KE)
> - **SMO 在所有电机上 OE 都无法收敛** — 实现有严重问题
> - **big_L 即使 OE 不大 (16% KE)，角度也完全不对 (93°)**

### 6.2 OE≈0 是角度精度的必要条件但非充分条件

> [!WARNING]
> **关键发现**：big_L 的 AF OE = 16.4% KE（不算很大），但角度误差 93°！
>
> 原因：OE = KE − |ψ_AF| 只反映磁链**幅值**误差，不反映**方向**误差。
> 当 KE/L < 1 时，电感电流项 L·i 在 ψ_s = ψ_AF + L·i 中占主导，
> 即使 |ψ_AF| ≈ KE（OE≈0），ψ_AF 的方向可以与真实磁链方向差 90°。

### 6.3 各电机的根本限制

| 电机 | 根本限制 | 解决方案 |
|:---|:---|:---|
| **small_L** | ✅ 无限制 | AF 即可全无感 (已验证 S3 稳定) |
| **servo** | ⚠️ 闭环 NSO-AF 耦合 | 需打破 NSO→FOC→AF 正反馈环路 |
| **big_L** | ❌ EMF 信号不足 (KE/L<1) | 必须 HFI，AF/PLL/SMO 均不可行 |

---

## 7. SMO 失败原因分析

SMO 在所有电机上 OE 均约为 -KE（即 |ψ_SMO| ≈ 2×KE），说明：
1. SMO 的滑模面定义有误 — 电流误差的符号函数方向可能反了
2. 或者 SMO 的积分器没有正确的磁链幅值约束
3. 需要重新检查 `observers_alt.py` 中 `SlidingModeObserver` 的实现

---

## 8. 下一步工作

1. **修复 SMO 实现** — 检查符号和积分方向
2. **servo 闭环解耦** — 用 αβ 系速度估计替代 dq 系 NSO，打破正反馈
3. **big_L 的 HFI** — 利用凸极性 (Ld≠Lq) 注入高频信号提取角度
4. **参数在线辨识** — servo 的 R 失配导致 OE 从 0.2%→8.0%，需要在线 R 估计
