# Coding Project 2 & 3 参考答案

## 说明

这份参考答案对应 `ACMSimPy/course/codingProject23.html`，风格参考 `homework1_solution.md`。  
重点不是给出某一次机器上的唯一数值结果，而是给出：

- 正确的实验设置和代码修改位置。
- 每个实验应观察到的定性结论。
- 报告中应写出的解释、公式和问答答案。
- 可用于批改的评分细则。

## GPT 阅卷说明和评分原则

这份文件是给 GPT 批改学生报告用的参考标准。批改时应先判断学生是否完成了项目要求，再按后面的 Project 2 / Project 3 细则给分。不要把本文当成唯一措辞模板；学生只要实验设置正确、图像证据充分、解释与代码自洽，就可以给相应分数。

建议 GPT 按下面流程阅卷：

1. **先确认提交内容。** 检查是否有 PDF/报告正文、关键图、代码修改说明、参数设置和问答答案。若缺少报告主体，只给能从代码或图中确认的部分分。
2. **先读图和代码，再读结论。** Coding Project 2 & 3 的核心证据是波形图、Bode 图、代码截图、参数表和仿真结果。文字结论若与图或代码冲突，应以图和代码为准，并在评语中指出冲突。
3. **按任务分别评分。** Project 2 和 Project 3 各按 `100` 分制独立评分；如果成绩表只有一个合并列，再按课程要求合并。没有明确要求时，可取两个项目平均值，但评语中应保留两个子分数。
4. **看实验是否可比。** ON/OFF、matched/mismatched、3rd/4th observer、不同带宽等对比必须使用相同速度命令、负载、带宽、时间窗口和作图变量。实验条件不一致时，不能用图直接支持结论。
5. **允许等价实现。** 学生可以用不同脚本结构、不同函数名或不同作图方式完成任务；只要注入点、观测量、控制逻辑和符号约定自洽，就不要因为代码形式不同扣分。
6. **过程分优先于关键词。** 公式有小符号错误但推理路线正确，可以保留过程分；只堆关键词但没有实验、图像或推导证据，不能给高分。
7. **符号自洽比固定正负号更重要。** 特别是 `xS[2]`、`-xS[2]`、`ACM.TLoad` 和扰动前馈的符号，必须结合学生定义判断。若学生明确把估计负载定义为 `-xS[2]`，并且图和代码一致，应视为正确。
8. **不要从低清图片中猜分。** 若图片无法读清轴、图例、单位或关键曲线，应标记为证据不足。可给文字分析分，但图像/实验验证分应扣。
9. **评语要可执行。** 最终批改评论应指出缺失的实验、错误的符号或不充分的图，而不是只写“理解不够”。建议给出 `Project 2: xx/100`、`Project 3: xx/100`，再列出主要扣分原因。

## 阅卷提示：读图片非常重要

Coding Project 2 & 3 的批改不能只看文字结论。学生提交中最有价值的证据通常在波形图、Bode 图、代码截图和参数截图里。读图时建议按下面顺序检查：

1. 先看标题、图例、坐标轴和单位。确认横轴是时间、带宽还是频率；纵轴是 rpm、A、V、Nm、dB 还是 degree。很多错误不是公式错，而是把 `xS[2]`、`-xS[2]`、`ACM.TLoad`、`CTRL.omega_r_mech` 画混了。
2. 再看实验场景是否一致。比较 ON/OFF、matched/mismatched、3rd/4th observer 时，带宽、负载、速度命令、PI 参数和时间窗口必须相同。若场景不同，图像不能直接支持结论。
3. 读图要抓趋势而不是只看截图是否漂亮。带宽扫描应看到上升时间随带宽下降、超调或噪声风险上升；解耦 ON 应在参数匹配时改善电流跟踪；参数失配时错误解耦可能比 OFF 更差；四阶观测器应在斜坡扰动中优于三阶。
4. 对扰动估计图要特别检查符号。当前代码通常有 `xS[2] ≈ -ACM.TLoad`，所以如果图中直接比较 `xS[2]` 与 `ACM.TLoad`，两者应符号相反；若比较 `-xS[2]` 与 `ACM.TLoad`，两者应同向并在低频/稳态处重合。
5. 对四阶观测器图要看 `xS[3]`。斜坡负载下，`xS[3]` 应收敛到扰动变化率；如果报告只画 `xS[2]` 而完全没有 `xS[3]`，不能证明学生真正验证了四阶观测器的核心优势。
6. 对 Bode 图要看低频极限和相位。若画的是 `T_L_hat/T_L = -xS[2]/T_L`，低频幅值应接近 `0 dB`，相位接近 `0 deg`；若画的是 `xS[2]/T_L`，低频相位会接近 `180 deg` 或 `-180 deg`。符号解释自洽即可。
7. 对代码截图要读关键行而不是只看“有改代码”。重点看 `FOC()` 中扰动前馈是否真正加到 `cmd_idq[1]`，是否有 `KA` 零值保护；看 `get_global_objects()` 是否覆盖 `ell1` 到 `ell4`，而不是修改 `@jitclass` spec 新增字段。
8. 手写或低清图片不要只依赖 OCR。需要放大原图，逐项确认题号、公式符号、上下标、正负号和图例。尤其是 `+/- omega Lq iq`、`xS[2]` 符号、`observer_order`、`ell4` 这些位置，OCR 很容易读错。

读图的基本技巧是：先确认“这张图在证明哪个任务”，再确认“横纵轴和图例是否能支持这个证明”，最后确认“图中的趋势是否和理论预期一致”。如果图和文字结论冲突，应以图和代码证据为主，并在批改评论中指出冲突。

注意：当前仓库中的 `tutorials_ep6_svpwm.py` 已经把扰动前馈接入了 `FOC()`，而且使用的是减号：

```python
CTRL.cmd_idq[1] -= CTRL.total_disrubance_feedforward / (1.5 * CTRL.npp * CTRL.KA)
```

这是与当前代码符号约定一致的写法。因为观测器方程中

$$
\dot{\hat \omega}_e
=
\ell_2 e_\theta
+
(T_{em}+\hat d_{to})\frac{n_{pp}}{J_s},
$$

而真实机械方程为

$$
J_s\dot\omega_m=T_{em}-T_L.
$$

所以正负载转矩 `T_L>0` 时，观测器扰动状态通常收敛到

$$
\hat d_{to}\approx -T_L.
$$

因此若要用它抵消负载，应令

$$
i_{q,ff}
=
-\frac{\hat d_{to}}{1.5\,n_{pp}K_A}.
$$

如果某份代码选择把 `total_disrubance_feedforward` 定义成 `-xS[2]`，则在 `cmd_idq[1]` 上用加号也是等价的。批改时应看符号闭环是否自洽，而不是只看加号或减号。

---

# Project 2 参考答案

## Part 1 带宽扫描

### 1.1 自动整定逻辑

当 `CL_SERIES_KP is None` 时，`Simulation_Benchmark` 会调用 `tuner.tunner_wrapper(d)` 自动计算电流环和速度环 PI 参数。核心关系是：

$$
K_{p,i}=2\pi L\,\mathrm{BW}_c,
\qquad
K_{i,i}=\frac{R}{L}.
$$

速度环整定中：

$$
K_{i,\omega}
=
\frac{\omega_{ci}}{\delta^2},
\qquad
K_{p,\omega}
=
\frac{J_s/n_{pp}}{1.5\,n_{pp}K_A}\,
\delta K_{i,\omega}.
$$

其中 `FOC_delta` 类似阻尼调节参数，`FOC_desired_VLBW_HZ` 是目标速度环带宽。

### 1.2 指标计算

速度阶跃响应可以用以下指标评价：

$$
t_r=t_{90\%}-t_{10\%},
$$

$$
M_p
=
\frac{y_{\max}-y_{ss}}{|y_{ss}|}\times 100\%,
$$

$$
e_{ss}=|y_{ss}-y_{ref}|.
$$

参考代码骨架：

```python
import copy
import numpy as np
import matplotlib.pyplot as plt
from tutorials_ep6_svpwm import Simulation_Benchmark

def step_metrics(t, y, ref, t0=0.2, t1=1.0):
    mask = (t >= t0) & (t <= t1)
    tt = t[mask]
    yy = y[mask]
    y0 = yy[0]
    yss = np.mean(yy[-max(10, len(yy)//10):])
    lo = y0 + 0.1 * (yss - y0)
    hi = y0 + 0.9 * (yss - y0)

    if yss >= y0:
        i10 = np.argmax(yy >= lo)
        i90 = np.argmax(yy >= hi)
    else:
        i10 = np.argmax(yy <= lo)
        i90 = np.argmax(yy <= hi)

    rise_time = tt[i90] - tt[i10]
    overshoot = (np.max(yy) - yss) / max(abs(yss), 1e-9) * 100
    ss_error = abs(yss - ref)
    return rise_time, overshoot, ss_error

bandwidths = [20, 40, 60, 80, 100, 120]
rise, overshoot, ess = [], [], []

for bw in bandwidths:
    d_run = copy.deepcopy(d)
    d_run["FOC_desired_VLBW_HZ"] = bw
    d_run["CL_SERIES_KP"] = None
    d_run["CL_SERIES_KI"] = None
    d_run["VL_SERIES_KP"] = None
    d_run["VL_SERIES_KI"] = None
    sim = Simulation_Benchmark(d_run)

    t = sim.global_machine_times
    y = sim.gdd["CTRL.omega_r_mech"]
    r = sim.gdd["CTRL.cmd_rpm"]
    tr, os, e = step_metrics(t, y, ref=np.mean(r[-len(r)//10:]))
    rise.append(tr)
    overshoot.append(os)
    ess.append(e)

fig, ax = plt.subplots(3, 1, sharex=True)
ax[0].plot(bandwidths, rise, "o-"); ax[0].set_ylabel("rise time [s]")
ax[1].plot(bandwidths, overshoot, "o-"); ax[1].set_ylabel("overshoot [%]")
ax[2].plot(bandwidths, ess, "o-"); ax[2].set_ylabel("ss error [rpm]")
ax[2].set_xlabel("FOC_desired_VLBW_HZ [Hz]")
fig.tight_layout()
```

### 1.3 预期结论

带宽增大时，上升时间通常变短，速度响应更快。代价是超调更容易增大，电压/电流命令更激烈，更容易遇到限幅，并且噪声敏感性和稳定裕度会变差。

`FOC_delta` 增大时，响应更保守，超调和振荡减小，但上升时间增加，负载扰动后的恢复速度变慢。

### 1.4 Q1 和 Q2

**Q1: When the system responds faster, what is the cost?**

更快的响应通常需要更高的环路增益，因此会带来更大的超调、更高的电压和电流峰值、更强的噪声放大，以及更小的稳定裕度。若电压或电流限幅被触发，名义带宽再高也无法体现，甚至会造成积分饱和和振荡。

**Q2: When the system is more stable, what is lost?**

更大的阻尼或更低的带宽会提高鲁棒性和稳定裕度，但损失跟踪速度和抗扰速度。阶跃命令到来时上升时间变长；负载扰动到来时，速度误差持续时间变长，暂态误差面积增大。

---

## Part 2 前馈解耦

### 2.1 dq 电压方程

PMSM 在 `dq` 坐标系下可写为：

$$
u_d
=
R i_d
+
L_d\frac{di_d}{dt}
-
\omega_e L_q i_q,
$$

$$
u_q
=
R i_q
+
L_q\frac{di_q}{dt}
+
\omega_e(K_A+L_q i_d).
$$

交叉耦合项为

$$
-\omega_e L_q i_q,
\qquad
\omega_e(K_A+L_q i_d).
$$

前馈解耦的目标是在 PI 输出后直接补偿这些项：

$$
u_{d,ff}=-\omega_e L_q i_q^*,
\qquad
u_{q,ff}=+\omega_e(K_A+L_q i_d^*).
$$

代码中对应：

```python
if CTRL.bool_apply_decoupling_voltages_to_current_regulation:
    decoupled_M_axis_voltage = -CTRL.omega_syn * CTRL.Lq * CTRL.cmd_idq[1]
    decoupled_T_axis_voltage =  CTRL.omega_syn * (CTRL.KA + CTRL.Lq * CTRL.cmd_idq[0])
    CTRL.cmd_udq[0] += decoupled_M_axis_voltage
    CTRL.cmd_udq[1] += decoupled_T_axis_voltage
```

### 2.2 预期实验现象

在参数匹配、速度较高、交叉耦合明显时，开启解耦通常会让 `i_d` 和 `i_q` 更接近给定值，电流动态更干净，速度响应也更接近设计带宽。

关闭解耦时，PI 必须把交叉耦合当作扰动处理。因此代码会在解耦关闭时把电流环积分增益乘以

```python
FOC_CL_KI_factor_when__bool_apply_decoupling_voltages_to_current_regulation__is_False
```

默认值为 `10`，目的是增强低频扰动抑制，弥补没有前馈解耦时的电流误差。

### 2.3 Q3 和 Q4

**Q3: Does feedforward decoupling have side effects?**

有。前馈解耦依赖 `L_d`、`L_q`、`K_A`、速度和电角度等模型信息。参数准确时，它能提前抵消耦合项；参数错误时，它会主动注入错误电压。此时解耦不仅不能抵消扰动，还会制造额外扰动，表现为电流波形畸变、速度超调增加，甚至限幅更频繁。

**Q4: Why does the code boost Ki when decoupling is OFF?**

解耦关闭时，交叉耦合项没有被前馈抵消，只能由 PI 的反馈作用消除。增加 `K_i` 可以提高电流环对常值或低频扰动的抑制能力，让稳态电流误差更小。但积分增益过大也会增加超调、振荡和饱和风险。

---

## Part 3 参数失配

### 3.1 正确实验方法

先用标称参数自动整定 PI，然后固定这组 PI 参数，再改变电机参数运行仿真。这样才是在测试“控制器设计基于旧模型，但实际对象变了”。

参考代码：

```python
import copy

d_base = copy.deepcopy(d)
d_base["CL_SERIES_KP"] = None
d_base["CL_SERIES_KI"] = None
d_base["VL_SERIES_KP"] = None
d_base["VL_SERIES_KI"] = None
sim_base = Simulation_Benchmark(d_base)

saved_gains = {
    "CL_SERIES_KP": d_base["CL_SERIES_KP"],
    "CL_SERIES_KI": d_base["CL_SERIES_KI"],
    "VL_SERIES_KP": d_base["VL_SERIES_KP"],
    "VL_SERIES_KI": d_base["VL_SERIES_KI"],
}

d_mismatch = copy.deepcopy(d)
d_mismatch.update(saved_gains)
d_mismatch["init_R"]  = d["init_R"]  * 1.5
d_mismatch["init_Ld"] = d["init_Ld"] * 0.7
d_mismatch["init_Lq"] = d["init_Lq"] * 0.7
d_mismatch["init_Js"] = d["init_Js"] * 2.0
```

三组必须比较：

| 场景 | 参数 | 解耦 |
|---|---|---|
| A | matched | ON |
| B | mismatched | ON |
| C | mismatched | OFF |

### 3.2 预期结论

场景 A 通常最好。参数匹配且解耦开启时，交叉耦合项被较准确抵消。

场景 B 可能明显变差。因为前馈解耦使用错误参数，会把错误补偿电压加进系统，导致电流环看到更大的等效扰动。

场景 C 在某些失配下反而可能比 B 稳健。关闭解耦后，控制器不再主动注入错误模型项，而是让 PI 把耦合和失配统一当作扰动处理。性能可能慢一些，但更不容易被错误前馈带偏。

### 3.3 Q5 和 Q6

**Q5: Why does parameter mismatch hurt decoupling more than no-decoupling?**

无解耦时，模型误差主要体现为被控对象动态变化，PI 仍然基于反馈误差修正。解耦开启时，控制器会根据错误参数直接计算补偿电压，这等价于把模型误差主动注入控制输入。因此参数失配对前馈解耦的伤害通常比对纯反馈控制更直接。

**Q6: Which parameters are most sensitive?**

`L_d`、`L_q` 对电流环和解耦电压最敏感，因为解耦项中直接含有 `omega_syn * L_q * i_q` 和 `omega_syn * L_q * i_d`。高速时这些误差被速度放大。`K_A` 影响转矩常数和 q 轴 back-EMF 补偿，因此会影响速度环和扰动前馈的单位转换。`J_s` 主要改变机械惯性和速度环设计匹配，惯量变大时响应变慢、抗扰恢复变慢。`R` 对低速和稳态电流误差更敏感，但高速下通常不如电感和 active flux 相关参数显著。

---

## Part 4 基础扰动观测器

### 4.1 三阶观测器

代码中的三阶速度/扰动观测器可理解为：

$$
\dot{\hat\theta}_d
=
\ell_1 e_\theta+\hat\omega_e,
$$

$$
\dot{\hat\omega}_e
=
\ell_2 e_\theta
+
(T_{em}+\hat d_{to})\frac{n_{pp}}{J_s},
$$

$$
\dot{\hat d}_{to}
=
\ell_3 e_\theta.
$$

其中 `xS[2]` 是总扰动估计。它不仅包含负载转矩，也包含摩擦、惯量误差、解耦误差和其他未建模项。

### 4.2 扰动前馈

若使用当前代码符号约定，推荐写法是：

```python
CTRL.cmd_idq[1] = reg_speed.Out
if CTRL.use_disturbance_feedforward_rejection > 0:
    if CTRL.KA > 1e-12:
        CTRL.cmd_idq[1] -= CTRL.total_disrubance_feedforward / (1.5 * CTRL.npp * CTRL.KA)
```

当 `TLoad>0` 时，通常有 `xS[2]≈-TLoad`。所以 `-xS[2]/Kt` 是正的 q 轴电流补偿，可以产生额外电磁转矩抵消负载。

### 4.3 Q7 和 Q8

**Q7: Why can an observer compensate for modeling errors?**

因为观测器不需要知道扰动的精确来源，只需要从位置/速度误差中反推出“必须存在一个等效扰动，才能解释实际输出和模型预测的差异”。负载转矩、摩擦、参数失配、错误解耦电压都会进入这个等效扰动状态。局限是观测器带宽有限，变化太快的扰动会滞后；噪声会通过高增益注入估计状态；若扰动通道与模型假设不一致，补偿也只能近似有效。

**Q8: Mode 1 vs Mode 2: which is better?**

`mode 1` 直接使用 `xS[2]`，估计更平滑，噪声较小，适合常规负载补偿。`mode 2` 使用 `xS[2] + ell2 * output_error`，等于加入更强的瞬时校正，响应更快，但更容易把编码器噪声、角度误差和数值抖动放大到 q 轴电流命令中。低噪声、扰动突变明显时 `mode 2` 可能更快；噪声较大或带宽较高时 `mode 1` 更稳健。

---

# Project 3 参考答案

## Part 1 启用四阶观测器

### 1.1 三阶与四阶的区别

三阶观测器状态为：

$$
[\hat\theta_d,\hat\omega_e,\hat d_{to}].
$$

它的内部模型是假设扰动为常值：

$$
\dot d_{to}=0.
$$

四阶观测器增加一个扰动变化率状态：

$$
[\hat\theta_d,\hat\omega_e,\hat d_{to},\hat p_{to}],
$$

并假设

$$
\dot d_{to}=p_{to},
\qquad
\dot p_{to}=0.
$$

因此四阶观测器可以无稳态误差地跟踪斜坡扰动。

### 1.2 增益配置

三阶极点配置在

$$
(s+\omega_{ob})^3
=
s^3+3\omega_{ob}s^2+3\omega_{ob}^2s+\omega_{ob}^3.
$$

所以：

$$
\ell_1=3\omega_{ob},
\quad
\ell_2=3\omega_{ob}^2,
\quad
\ell_3=\omega_{ob}^3\frac{J_s}{n_{pp}},
\quad
\ell_4=0.
$$

四阶极点配置在

$$
(s+\omega_{ob})^4
=
s^4+4\omega_{ob}s^3+6\omega_{ob}^2s^2+4\omega_{ob}^3s+\omega_{ob}^4.
$$

所以：

$$
\ell_1=4\omega_{ob},
\quad
\ell_2=6\omega_{ob}^2,
\quad
\ell_3=4\omega_{ob}^3\frac{J_s}{n_{pp}},
\quad
\ell_4=\omega_{ob}^4.
$$

### 1.3 推荐代码修改

由于 `The_Motor_Controller` 是 `@jitclass`，不要新增字段。可以在 `get_global_objects()` 创建 `CTRL` 后直接覆盖已有的 `ell1` 到 `ell4`。

```python
CTRL = The_Motor_Controller(...)

omega_ob = d.get("omega_ob", 100.0)
observer_order = d.get("observer_order", 3)

if observer_order == 4:
    CTRL.ell1 = 4 * omega_ob
    CTRL.ell2 = 6 * omega_ob**2
    CTRL.ell3 = 4 * omega_ob**3 * d["init_Js"] / d["init_npp"]
    CTRL.ell4 = omega_ob**4
else:
    CTRL.ell1 = 3 * omega_ob
    CTRL.ell2 = 3 * omega_ob**2
    CTRL.ell3 = omega_ob**3 * d["init_Js"] / d["init_npp"]
    CTRL.ell4 = 0.0

CTRL.index_separate_speed_estimation = d.get(
    "CTRL.index_separate_speed_estimation", 0
)
CTRL.use_disturbance_feedforward_rejection = d.get(
    "use_disturbance_feedforward_rejection", 0
)
```

同时在 `d` 中加入：

```python
"observer_order": 4,
"omega_ob": 100.0,
"CTRL.index_separate_speed_estimation": 1,
"use_disturbance_feedforward_rejection": 1,
```

### 1.4 阶跃负载验证

阶跃负载：

```python
"user_system_input_code": (
    "CTRL.cmd_rpm = 50\n"
    "if ii >= 2: ACM.TLoad = 0.2"
)
```

预期结果：

- 三阶和四阶观测器都能让 `xS[2]` 收敛到常值扰动。
- 若画的是 `xS[2]`，它应约等于 `-ACM.TLoad`。
- 若画的是估计负载 `T_L_hat=-xS[2]`，它应收敛到 `ACM.TLoad`。
- 四阶可能有稍快或稍振荡的暂态，但阶跃稳态能力不是它相对三阶的主要优势。

---

## Part 2 接入扰动前馈

扰动估计是转矩量，q 轴电流与电磁转矩之间的关系为：

$$
T_{em}=1.5\,n_{pp}K_A i_q.
$$

所以扰动补偿电流为：

$$
i_{q,ff}
=
-\frac{\hat d_{to}}{1.5\,n_{pp}K_A}.
$$

当前代码符号下推荐：

```python
CTRL.cmd_idq[1] = reg_speed.Out
if CTRL.use_disturbance_feedforward_rejection > 0:
    if CTRL.KA > 1e-12:
        CTRL.cmd_idq[1] -= CTRL.total_disrubance_feedforward / (
            1.5 * CTRL.npp * CTRL.KA
        )
```

预期结果：

- 前馈关闭时，负载扰动先造成速度下跌，再由速度 PI 慢慢恢复。
- 前馈开启时，观测器一旦估计出扰动，q 轴电流命令会提前增加，速度下跌幅度明显减小。
- 若符号接反，表现会非常明显：负载一加，速度误差变大，q 轴电流朝错误方向变化。

---

## Part 3 斜坡负载与阶跃负载对比

### 3.1 四种场景

| 场景 | 观测器 | 前馈 | 斜坡负载下预期表现 |
|---|---|---|---|
| A | 3rd | OFF | 只有 PI，速度误差最大，负载增加时速度持续偏离 |
| B | 3rd | ON | 前馈有帮助，但 `xS[2]` 对斜坡有滞后，仍有稳态跟踪误差 |
| C | 4th | OFF | `xS[2]` 估计更准，但不接入补偿时速度仍主要靠 PI 恢复 |
| D | 4th | ON | 最好，`xS[2]` 跟踪斜坡，`xS[3]` 收敛到扰动斜率，速度误差最小 |

斜坡负载设置：

```python
"user_system_input_code": (
    "CTRL.cmd_rpm = 50\n"
    "if ii >= 2: ACM.TLoad = 0.05 * (ii - 2) * d['TIME_SLICE']"
)
```

### 3.2 需要画的图

每个场景建议画三行子图：

```python
ax[0].plot(t, gdd["CTRL.cmd_rpm"], label="cmd rpm")
ax[0].plot(t, gdd["CTRL.omega_r_mech"], label="actual rpm")

ax[1].plot(t, gdd["ACM.TLoad"], label="TLoad")
ax[1].plot(t, -gdd["CTRL.xS[2]"], label="-xS[2] = estimated load")

ax[2].plot(t, gdd["CTRL.xS[3]"], label="xS[3] disturbance rate")
```

如果直接画 `xS[2]` 而不是 `-xS[2]`，需要在图注中说明符号相反。

### 3.3 阶跃负载对比

阶跃负载：

```python
"user_system_input_code": (
    "CTRL.cmd_rpm = 50\n"
    "if ii >= 2: ACM.TLoad = 0.2"
)
```

预期结论：

- 阶跃负载是常值扰动，三阶观测器的内部模型已经包含它，所以三阶和四阶最终都可零稳态误差估计。
- 四阶的优势不在常值扰动，而在斜坡或更慢变化的时变扰动。
- 若四阶带宽太高，阶跃时可能比三阶更嘈杂或更振荡。

---

## Part 4 性能分析

闭环输出可分解为：

$$
\Omega(s)
=
\Phi_r(s)\Omega^*(s)
+
\Phi_d(s)d_n(s)
+
\Phi_n(s)(-s\delta_p(s)).
$$

### 4.1 Q1: Command tracking

三阶和四阶观测器在相同速度环 PI、相同带宽、相同前馈连接方式下，对指令跟踪本身影响不应很大。指令跟踪主要由速度环 PI 和电流环带宽决定。观测器阶数主要影响扰动通道，而不是参考输入通道。

实际仿真中若四阶的阶跃跟踪略有差异，通常来自两个原因：一是四阶估计状态更快、更容易参与前馈；二是更高阶观测器对角度误差和噪声更敏感，可能间接影响 q 轴电流命令。

### 4.2 Q2: Disturbance rejection

斜坡扰动可写为

$$
d_{to}(t)=a t.
$$

三阶观测器假设扰动常值，因此对斜坡扰动存在稳态估计滞后，误差量级与斜率成正比、与观测器带宽成反比：

$$
e_d \propto \frac{a}{\omega_{ob}}.
$$

四阶观测器包含 `p_to`，也就是扰动斜率的内部模型。对理想斜坡扰动，四阶可使

$$
\hat p_{to}\to a,
\qquad
\hat d_{to}\to d_{to},
$$

因此斜坡扰动估计稳态误差为零。前馈开启时，四阶的速度稳态误差也应显著小于三阶。

### 4.3 Q3: Noise sensitivity

四阶观测器在相同 `omega_ob` 下通常比三阶更容易放大噪声，因为它多估计了一个微分性质更强的状态 `xS[3]`。增大 `omega_ob` 会提高扰动跟踪速度，但也会把角度测量噪声、离散积分误差和数值抖动更强地注入 `xS[2]`、`xS[3]`，进而影响 q 轴电流前馈。

合理选择是：扰动变化慢、传感器噪声大时降低 `omega_ob`；负载变化快、噪声较小时可以提高 `omega_ob`。四阶不是无条件更好，它用噪声敏感性换取对斜坡扰动的更高跟踪阶次。

### 4.4 与参数失配结合

高阶观测器对参数失配也有帮助，因为参数失配会进入总扰动项。例如错误的 `L_q` 解耦、电机惯量误差、转矩常数误差都会表现为模型预测与实际速度之间的偏差，观测器可以把它们折算进 `d_to`。

但是四阶只对“可被低阶多项式近似的慢变扰动”更有优势。若失配引起的是高频振荡、饱和、强非线性或错误符号前馈，四阶观测器不能完全补救，甚至可能因噪声和错误估计导致更差表现。

---

## 频域实验参考答案

### 1. 传递路径 `T_L_hat / T_L`

建议定义估计负载为：

$$
\hat T_L=-xS[2].
$$

这样低频下应有：

$$
\frac{\hat T_L}{T_L}\approx 1,
$$

即 Bode 幅值约 `0 dB`，相位约 `0 deg`。

若直接画 `xS[2]/T_L`，低频幅值仍约 `0 dB`，但相位会多出约 `180 deg`，因为 `xS[2]≈-T_L`。

### 2. 频率响应预期

低频区域：三阶和四阶都能较好估计常值或慢变负载，`|T_L_hat/T_L|` 接近 `1`。

中频区域：估计开始出现幅值衰减和相位滞后。观测器带宽越高，截止频率越高，但噪声越大。

高频区域：观测器无法及时跟踪快速负载变化，估计幅值下降，相位滞后增加。四阶在扰动斜率估计方面更强，但高频噪声也更明显。

### 3. `omega=0 rpm`、负载幅值 `0.15 Nm` 扫频

实验设置：

```python
CTRL.cmd_rpm = 0
ACM.TLoad = 0.15 * np.sin(2 * np.pi * f * t)
```

报告中应给出：

- 每个扫频点的时域波形：`ACM.TLoad` 与 `-CTRL.xS[2]`。
- Bode 幅值：`20 log10(|T_L_hat/T_L|)`。
- Bode 相位：`angle(T_L_hat)-angle(T_L)`。

预期结论：

- 低频时估计负载与实际负载几乎重合。
- 频率接近或超过观测器带宽后，估计负载开始滞后并衰减。
- 增大 `omega_ob` 可改善带宽内跟踪，但会让 `xS[2]` 和 `xS[3]` 更嘈杂。

---

# 最终报告建议结构

1. Project 2 Part 1：带宽扫描方法、指标定义、结果图和 Q1-Q2。
2. Project 2 Part 2：解耦公式、ON/OFF 波形、Q3-Q4。
3. Project 2 Part 3：参数失配方法、三场景对比、Q5-Q6。
4. Project 2 Part 4：三阶观测器、扰动前馈、负载阶跃与失配补偿、Q7-Q8。
5. Project 3 Part 1：三阶/四阶观测器公式、增益、代码修改。
6. Project 3 Part 2：扰动前馈接线和符号验证。
7. Project 3 Part 3：阶跃/斜坡负载四场景对比。
8. Project 3 Part 4：Q1-Q3、噪声权衡、频域扫频与 Bode 图。

---

# 评分细则建议

## Project 2 总分 100

| 部分 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| Part 1 带宽扫描 | 30 | 正确扫描多个 `FOC_desired_VLBW_HZ`，计算上升时间、超调、稳态误差，并解释带宽/阻尼权衡 | 只给波形不算指标扣 8-12 分；没有强制重新整定扣 5 分；Q1-Q2 解释空泛扣 4 分 |
| Part 2 解耦分析 | 30 | 写出 dq 耦合项和前馈补偿项，完成 ON/OFF 对比，解释副作用与 `Ki` 放大原因 | 只说“解耦更好”不解释参数依赖扣 6 分；未比较电流/电压波形扣 6-10 分 |
| Part 3 参数失配 | 30 | 先保存标称整定增益，再改变参数；完成 A/B/C 三场景对比，解释为什么错误解耦可能更差 | 直接改参数并重新整定扣 8-12 分；没有比较 mismatched ON/OFF 扣 8 分 |
| Part 4 基础观测器 | 10 | 正确说明 `xS[2]` 总扰动含义，展示负载阶跃估计和前馈补偿效果，回答 Q7-Q8 | 符号解释错误扣 3 分；只开观测器不接前馈扣 3 分 |

## Project 3 总分 100

| 部分 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| Part 1 四阶观测器 | 20 | 正确配置 `observer_order`，覆盖 `ell1`-`ell4`，不修改 jitclass spec，阶跃负载验证三阶/四阶均能估计常值扰动 | 新增 jitclass 字段导致运行失败扣 8 分；四阶增益系数错误扣 5 分 |
| Part 2 扰动前馈 | 20 | 在 `FOC()` 中正确把扰动估计转换为 q 轴电流，符号自洽，并展示前馈开启后速度跌落减小 | 符号接反扣 8-12 分；没有 `KA` 零值保护扣 2 分 |
| Part 3 斜坡/阶跃对比 | 30 | 完成 A-D 四场景，分别测试斜坡和阶跃负载，说明四阶优势只在时变扰动中明显 | 只做阶跃不做斜坡扣 10 分；没有 `xS[3]` 变化率图扣 5 分 |
| Part 4 性能与频域分析 | 30 | 回答 Q1-Q3，说明跟踪/抗扰/噪声权衡，给出 `T_L_hat/T_L` Bode 图和 `0 rpm, 0.15 Nm` 扫频时域波形 | 没有 Bode 图扣 8 分；没有 0 rpm 扫频扣 8 分；未说明 `xS[2]` 符号扣 3 分 |

---

# 核心结论

Project 2 的核心是：高带宽和前馈解耦能提高响应速度，但鲁棒性受参数准确性限制；纯反馈慢一些，但在模型不准时更稳健。

Project 3 的核心是：三阶观测器内置“常值扰动”模型，因此能处理阶跃负载；四阶观测器额外估计扰动变化率，因此能处理斜坡负载。更高阶和更高带宽可以提升时变扰动跟踪能力，但代价是更强的噪声敏感性。
