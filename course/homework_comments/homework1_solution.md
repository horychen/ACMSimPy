# Homework 1 参考答案

## From Coordinates to the Unified Machine ODE

> 说明  
> 1. 以下答案按题目顺序作答，可直接作为电子版作业初稿。  
> 2. 若老师要求保留一般初相位，可把文中的 `\omega t` 统一替换为 `\omega t + \phi_0`；不影响主结论。  
> 3. 记 active flux 为 `\psi_\text{AF}`，代码中对应状态 `KA = x[2]`。  

---

## GPT 阅卷说明和评分原则

这份文件是给 GPT 批改 Homework 1 用的参考答案和评分依据。批改时不要只做字符串匹配；应判断学生是否理解坐标变换、旋转系求导、磁链建模、active-flux 统一模型和五维 ODE 的推导逻辑。学生的符号可以和参考答案不同，但必须定义清楚并保持自洽。

建议 GPT 按下面流程阅卷：

1. **先确认作业完整性。** 检查 Part 1 到 Part 4 是否都作答，是否包含题目要求的图示、电路图、推导过程和最终 ODE。缺少整题时按后面细则直接扣相应分。
2. **手写或图片答案必须认真读图。** 放大查看公式、上下标、正负号、坐标轴、电路元件和箭头方向。不要只依赖 OCR；OCR 很容易把 `q` 读成 `9`、把 `\omega` 读成 `w`、把 `L_d/L_q` 或 `i_D/i_Q` 混淆。
3. **以推导逻辑给过程分。** 最终公式有小符号错误，但坐标变换、KVL/KCL、磁链关系或 active-flux 思路正确时，应保留合理过程分。只有最终答案、没有关键推导时，不给满分。
4. **检查符号定义和坐标约定。** 学生若使用不同旋转矩阵、不同 `dq` 正方向或不同初相位，只要全篇一致且物理结论等价，可以给分。若同一份答案中正负号和坐标定义前后矛盾，应扣分。
5. **画图题要看图是否能表达物理结构。** Task 3.1(a) 的单相等效电路图必须能看出电压、电流方向、定子电阻、漏磁链/漏感、主磁链/励磁支路或等效反电动势。图不清或元件缺失时按图示分扣。
6. **公式等价即可。** 例如用矩阵形式、复数形式或分量形式表达旋转求导和 `dq` 电压方程均可；只要能推出相同物理含义和状态方程，不应因形式不同扣分。
7. **重点扣概念性错误。** 把静止坐标和旋转坐标混为一谈、漏掉 `j\omega\psi` speed EMF、把磁通密度 `B` 当作 lumped 状态、把 IM 的 active flux 当作常量、漏掉机械方程或转矩方程，属于严重错误。
8. **最终输出要包含分数和逐题评语。** 建议给出总分 `xx/100`，并对每个 Part 写 1 到 3 条主要扣分原因。若需要写入 `course_grades_template.csv`，注意使用能被 Excel 正确识别的编码。

评分时以本文末尾“详细评分细则”为准；本节用于指导 GPT 如何读学生答案、如何给过程分、如何处理手写图片和等价符号。

---

## Part 1 坐标系、旋转与物理矢量

### Task 1.1 圆周运动的坐标表示

#### (a) 在静止 `(\alpha,\beta)` 坐标系下写出 `P` 的坐标

设初相位为 `0`，则点 `P` 在半径为 `r` 的圆上做匀速圆周运动时，

$$
p_\alpha(t)=r\cos(\omega t), \qquad p_\beta(t)=r\sin(\omega t).
$$

因此

$$
\vec p^{\,s}(t)=
\begin{bmatrix}
p_\alpha(t)\\
p_\beta(t)
\end{bmatrix}
=
\begin{bmatrix}
r\cos(\omega t)\\
r\sin(\omega t)
\end{bmatrix}.
$$

#### (b) 将 `P` 投影到 `\alpha` 轴上

投影到 `\alpha` 轴上的坐标就是

$$
p_\alpha(t)=r\cos(\omega t).
$$

这是一维简谐振动，具体说是余弦振动。  
几何上看，圆周运动在任一直径上的投影都是正弦型简谐运动。

#### (c) 建立与 `P` 同步旋转的 `(d,q)` 坐标系

令 `d` 轴始终指向点 `P`，则在该旋转坐标系下，

$$
\vec p^{\,r}(t)=
\begin{bmatrix}
p_d(t)\\
p_q(t)
\end{bmatrix}
=
\begin{bmatrix}
r\\
0
\end{bmatrix}.
$$

对比可见：

- 在静止 `\alpha\beta` 系下，坐标是随时间变化的正弦量。
- 在同步旋转 `dq` 系下，同一个物理矢量的坐标变成常数。

所以 `(d,q)` 坐标系更简单。这正是电机建模里要做坐标变换的核心原因。

---

### Task 1.2 物理矢量与坐标无关性

#### (a) 写出旋转矩阵并验证 `\mathbf R(-\theta)=\mathbf R^T(\theta)`

取标准二维旋转矩阵

$$
\mathbf R(\theta)=
\begin{bmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{bmatrix}.
$$

则

$$
\mathbf R(-\theta)=
\begin{bmatrix}
\cos\theta & \sin\theta\\
-\sin\theta & \cos\theta
\end{bmatrix}.
$$

另一方面，

$$
\mathbf R^T(\theta)=
\begin{bmatrix}
\cos\theta & \sin\theta\\
-\sin\theta & \cos\theta
\end{bmatrix}.
$$

因此

$$
\mathbf R(-\theta)=\mathbf R^T(\theta).
$$

这说明旋转矩阵是正交矩阵，逆变换等于转置。

#### (b) 验证从 `\alpha\beta` 变换到 `dq` 后得到 Task 1.1(c) 的结果

由 Task 1.1(a) 得

$$
\vec p^{\,s}(t)=
\begin{bmatrix}
r\cos(\omega t)\\
r\sin(\omega t)
\end{bmatrix}.
$$

变换到同步旋转坐标系：

$$
\vec p^{\,r}(t)=\mathbf R(-\omega t)\vec p^{\,s}(t).
$$

代入

$$
\vec p^{\,r}(t)=
\begin{bmatrix}
\cos\omega t & \sin\omega t\\
-\sin\omega t & \cos\omega t
\end{bmatrix}
\begin{bmatrix}
r\cos\omega t\\
r\sin\omega t
\end{bmatrix}
=
\begin{bmatrix}
r(\cos^2\omega t+\sin^2\omega t)\\
r(-\sin\omega t\cos\omega t+\cos\omega t\sin\omega t)
\end{bmatrix}
=
\begin{bmatrix}
r\\
0
\end{bmatrix}.
$$

与 Task 1.1(c) 完全一致。

---

### Task 1.3 旋转坐标系中的求导算子

#### (a) 麦克风与桌子的平移类比

设麦克风相对墙的绝对位置为

$$
x_\text{abs}=X+Y,
$$

其中：

- `X` 是桌子相对墙的位置。
- `Y` 是麦克风相对桌子边缘的位置。

若麦克风在桌子参考系中静止，则

$$
\dot Y=0,\qquad \ddot Y=0.
$$

于是绝对速度和绝对加速度分别为

$$
\dot x_\text{abs}=\dot X+\dot Y=\dot X=v,
$$

$$
\ddot x_\text{abs}=\ddot X+\ddot Y=\ddot X=\dot v.
$$

结论：

- 若桌子匀速运动，`v` 为常数，则 `\ddot x_\text{abs}=0`。
- 若桌子加速运动，则 `\ddot x_\text{abs}=\dot v`，这个额外项来自参考系本身的运动。

在桌子参考系里，这个“由参考系运动带来的额外项”表现为惯性力。  
严格地说，平移坐标系只有在加速度不为零时才出现惯性力；而旋转坐标系即使角速度恒定，也会因为基向量本身在转动而出现额外项。

#### (b) 推广到旋转坐标系，写出修正求导算子

设物理矢量 `\hat{\vec p}` 在旋转坐标系中的复数表示为

$$
\vec p^{\,r}=p_d+j p_q.
$$

若坐标系以角速度 `\omega` 旋转，则矢量在惯性系中的绝对导数与在旋转系中的相对导数满足

$$
\left.\frac{d\vec p}{dt}\right|_{\text{inertial}}
=
\left.\frac{d\vec p}{dt}\right|_{\text{rotating}}
+j\omega\,\vec p.
$$

因此修正后的求导算子为

$$
\boxed{
\left.\frac{d}{dt}\right|_{\text{inertial}}
=
\left.\frac{d}{dt}\right|_{\text{rotating}}
+j\omega
}
$$

这就是课上反复强调的

$$
s \;\longrightarrow\; s+j\omega.
$$

如果写成实数二维向量形式，则有

$$
\left.\frac{d}{dt}
\begin{bmatrix}
p_d\\
p_q
\end{bmatrix}
\right|_{\text{inertial}}
=
\left.\frac{d}{dt}
\begin{bmatrix}
p_d\\
p_q
\end{bmatrix}
\right|_{\text{rotating}}
+
\omega
\begin{bmatrix}
0 & -1\\
1 & 0
\end{bmatrix}
\begin{bmatrix}
p_d\\
p_q
\end{bmatrix}.
$$

---

### Task 1.4 反电动势的例子

#### (a) 解释“力 ↔ 反电动势”的类比

更严格的能量观点下，机械系统和电磁系统的类比如下：

- 力 `F` 对应电压 `u`，都是“广义驱动力/努力变量”。
- 速度 `v` 对应电流 `i`，都是“流变量”。
- 动量 `p=mv` 对应磁链 `\psi=Li`，都是储能变量。

于是：

$$
F=\frac{dp}{dt}, \qquad u=\frac{d\psi}{dt}.
$$

若再把 `p=mv`、`v=\dot x` 代入，就得到

$$
F=m\ddot x.
$$

因此“力 ↔ 反电动势”的类比本质上是：

- 机械里，质量储存动能，状态变化会产生惯性反作用。
- 电磁里，电感/磁链储存磁场能，状态变化会产生感应电压反作用。

两者都体现了“储能元件反对其状态突变”的规律。

#### (b) 写出旋转系中的反电动势表达式，并说明 `j\omega\vec\psi` 的物理意义

法拉第定律写成矢量形式：

$$
\vec e=-\left.\frac{d\vec\psi}{dt}\right|_{\text{inertial}}.
$$

利用 Task 1.3 的修正算子，

$$
\vec e^{\,r}
=
-\left(
\left.\frac{d\vec\psi^{\,r}}{dt}\right|_{\text{rotating}}
+j\omega \vec\psi^{\,r}
\right).
$$

若

$$
\vec\psi^{\,r}=\psi_d+j\psi_q,
$$

则

$$
\vec e^{\,r}
=
-\Big[
(\dot\psi_d-\omega\psi_q)+j(\dot\psi_q+\omega\psi_d)
\Big].
$$

所以分量形式为

$$
e_d=-\dot\psi_d+\omega\psi_q,
\qquad
e_q=-\dot\psi_q-\omega\psi_d.
$$

其中额外项

$$
j\omega\vec\psi
$$

代表的就是速度电动势 speed EMF。  
它不是因为磁链幅值本身在变，而是因为磁链矢量相对于观察坐标系在转动，从而“切割导体”产生感应电压。

---

## Part 2 从磁场到电路：磁通与磁链

### Task 2.1

#### (a) 为什么在电路分析中更偏好磁链 `\psi` 而不是磁通密度 `\vec B`

原因如下：

- `\vec B` 是场量，依赖空间位置 `\vec r`，要通过面积分才能得到总磁通。
- 电路微分方程希望状态量只依赖时间 `t`，而不想直接处理空间分布。
- 磁链 `\psi=N\Phi` 把空间分布信息“浓缩”为单个绕组可感受到的 lumped quantity。
- 法拉第定律可以直接写成 `u=d\psi/dt`，因此 `\psi(t)` 很自然地成为电路状态变量。

所以：

- 场论里更基础的是 `\vec B`。
- 电路建模里更实用的是 `\psi`。

#### (b) 若 `\psi=Li` 且 `L` 为常数，写出 `u` 关于 `i` 的表达式

由法拉第定律

$$
u=\frac{d\psi}{dt}.
$$

又因为

$$
\psi=Li,\qquad L=\text{const},
$$

所以

$$
u=\frac{d(Li)}{dt}=L\frac{di}{dt}.
$$

这就是理想电感电压方程。

#### (c) 若 `L` 不是常数，`u=\frac{d(Li)}{dt}` 有几项，哪一项对应反电动势

乘积求导：

$$
u=\frac{d(Li)}{dt}
=L\frac{di}{dt}+i\frac{dL}{dt}.
$$

因此一共有两项：

1. `L \, di/dt`：电流变化引起的感应电压。
2. `i \, dL/dt`：电感变化引起的感应电压。

若 `L` 随转子位置变化，则

$$
\frac{dL}{dt}=\frac{dL}{d\theta}\frac{d\theta}{dt}
=\frac{dL}{d\theta}\omega.
$$

所以第二项

$$
i\frac{dL}{dt}
$$

就是典型的速度相关反电动势项。

---

## Part 3 电机等效电路与 KVL/KCL

### Task 3.1 单相等效电路

#### (a) 画出单相等效电路图

可用最简串联 lumped model 表示为：

```text
      i
 + ──▶──[ R_s ]──[ L_sigma , psi_sigma ]──[ e_m = d psi_m / dt ]── −
 |                                                                  |
 └──────────────────────────── phase return ────────────────────────┘
```

也可把最后一个感应电压源写成励磁电感 `L_m`：

```text
      i
 + ──▶──[ R_s ]──[ L_sigma , psi_sigma ]──[ L_m , psi_m ]── −
 |                                                            |
 └────────────────────── phase return ────────────────────────┘
```

其中：

- `u`：相电压。
- `i`：定子相电流。
- `R_s`：定子铜阻。
- `L_\sigma`：漏感，对应漏磁链 `\psi_\sigma=L_\sigma i`。
- `L_m`：励磁电感，对应主磁链 `\psi_m`。

对这道题，画成“`R_s + 漏感 + 主磁链感应电压`”的串联等效最容易写出 KVL。

#### (b) 写出 KVL，并说明每一项物理含义

由基尔霍夫电压定律，

$$
u=R_s i+\frac{d\psi_\sigma}{dt}+\frac{d\psi_m}{dt}.
$$

若线性近似下 `\psi_\sigma=L_\sigma i`，则

$$
u=R_s i+L_\sigma\frac{di}{dt}+\frac{d\psi_m}{dt}.
$$

各项物理意义：

- `R_s i`：绕组铜损造成的电阻压降。
- `d\psi_\sigma/dt`：漏磁链变化引起的感应电压。
- `d\psi_m/dt`：主磁链变化引起的感应电压，也就是与转子耦合相关的 back-EMF。

#### (c) 解释 `L_\sigma` 和 `L_m` 的物理区别

- `L_\sigma` 对应漏磁通路径，磁通主要局限在定子附近，不与转子有效交链。
- `L_m` 对应主磁通路径，磁通会穿过气隙并与转子磁场或转子导体有效交链。

所以：

- “不与转子交链”的是 `L_\sigma`。
- “与永磁体/转子感应绕组交链”的是 `L_m`。

---

### Task 3.2 从三相到矢量：KCL 的约束

#### (a) 星形接法下写出 KCL 约束

若中性点相连且无中性线电流，则

$$
i_A+i_B+i_C=0.
$$

#### (b) 为什么可以用二维矢量 `\vec i^s` 等价表示三相电流

三相电流原本位于三维空间：

$$
\begin{bmatrix}
i_A\\
i_B\\
i_C
\end{bmatrix}.
$$

但由于约束

$$
i_A+i_B+i_C=0,
$$

这三个分量并不独立，而是被限制在三维空间中的一个二维平面上。  
因此真实自由度只有 `2` 个。

这时可以用 Clarke 变换把该二维平面映射为 `\alpha\beta` 平面：

$$
\begin{bmatrix}
i_\alpha\\
i_\beta
\end{bmatrix}
=
\frac{2}{3}
\begin{bmatrix}
1 & -\frac12 & -\frac12\\
0 & \frac{\sqrt3}{2} & -\frac{\sqrt3}{2}
\end{bmatrix}
\begin{bmatrix}
i_A\\
i_B\\
i_C
\end{bmatrix}.
$$

因此二维矢量

$$
\vec i^{\,s}=
\begin{bmatrix}
i_\alpha\\
i_\beta
\end{bmatrix}
$$

就足以无损表示无零序分量的三相电流系统。

#### (c) `\alpha\beta \to dq` 变换与 Part 1 圆周运动的联系

联系完全相同：

- Part 1 中，圆周运动的矢量在静止 `\alpha\beta` 系里表现为正弦变化。
- 若坐标系跟着它一起转，则在 `dq` 系里变成常数。

电机稳态三相电流经 Clarke 变换后，在 `\alpha\beta` 平面中正好构成一个旋转电流矢量。  
再做 Park 变换到同步旋转 `dq` 系后，稳态时电流分量会接近常数：

$$
\vec i^{\,r}=\mathbf R(-\theta)\vec i^{\,s}.
$$

这样做的好处是：

- 把正弦量变成直流量。
- 把含有三角函数调制的时变模型，变成系数更简单的微分方程。
- 便于控制器设计，因为 PI 更擅长调节“常值”而不是“正弦”。

---

## Part 4 综合推导：五条统一 ODE

### Task 4.1 电压方程推导

#### (a) 从 KVL 写成磁链形式

对定子绕组写 KVL：

$$
u=R_s i+\frac{d\psi}{dt}.
$$

这是电机电压方程最基本的 lumped form。

#### (b) 变换到以 `\omega_\text{syn}` 旋转的 `dq` 坐标系

在旋转坐标系中，矢量形式可写成

$$
\vec u^{\,r}
=
R_s\vec i^{\,r}
+
\left(
\frac{d}{dt}+j\omega_\text{syn}
\right)\vec\psi^{\,r}.
$$

令

$$
\vec\psi^{\,r}=\psi_d+j\psi_q,
\qquad
\vec i^{\,r}=i_D+j i_Q,
\qquad
\vec u^{\,r}=u_D+j u_Q.
$$

展开得

$$
u_D=R_s i_D+\frac{d\psi_d}{dt}-\omega_\text{syn}\psi_q,
$$

$$
u_Q=R_s i_Q+\frac{d\psi_q}{dt}+\omega_\text{syn}\psi_d.
$$

又因为

$$
\psi_d=\psi_\text{AF}+L_q i_D,
\qquad
\psi_q=L_q i_Q,
$$

所以

$$
u_D
=
R_s i_D+\frac{d\psi_\text{AF}}{dt}+L_q\frac{di_D}{dt}-\omega_\text{syn}L_q i_Q,
$$

$$
u_Q
=
R_s i_Q+L_q\frac{di_Q}{dt}+\omega_\text{syn}(L_q i_D+\psi_\text{AF}).
$$

#### (c) 用 `(\psi_\text{AF},i_D,i_Q)` 作为状态量，解出 `di_D/dt` 和 `di_Q/dt`

由上两式直接整理：

$$
\boxed{
\frac{di_D}{dt}
=
\frac{u_D-R_s i_D+\omega_\text{syn}L_q i_Q-\frac{d\psi_\text{AF}}{dt}}{L_q}
}
$$

$$
\boxed{
\frac{di_Q}{dt}
=
\frac{u_Q-R_s i_Q-\omega_\text{syn}L_q i_D-\omega_\text{syn}\psi_\text{AF}}{L_q}
}
$$

这就是代码里的 `\dot x[3]` 和 `\dot x[4]` 的统一来源。

---

### Task 4.2 Active flux 动力学

#### (a) PMSM 情况：`R_\text{req}=0`

PMSM 中

$$
\psi_\text{AF}=(L_d-L_q)i_D+\psi_\text{PM}.
$$

因为 `\psi_\text{PM}` 为常数，所以

$$
\frac{d\psi_\text{AF}}{dt}
=(L_d-L_q)\frac{di_D}{dt}.
$$

这就是代码中 `\dot x[2]` 的 PMSM 分支。

再代回 `d` 轴电压方程：

$$
u_D
=
R_s i_D+(L_d-L_q)\frac{di_D}{dt}+L_q\frac{di_D}{dt}-\omega_\text{syn}L_q i_Q
$$

于是

$$
u_D=R_s i_D+L_d\frac{di_D}{dt}-\omega_\text{syn}L_q i_Q,
$$

得到

$$
\boxed{
\frac{di_D}{dt}
=
\frac{u_D-R_s i_D+\omega_\text{syn}L_q i_Q}{L_d}
}
$$

并且

$$
\boxed{
\frac{d\psi_\text{AF}}{dt}=(L_d-L_q)\frac{di_D}{dt}
}
$$

这与代码分支完全一致。

#### (b) IM 情况：`R_\text{req}>0`

感应电机没有永磁体，转子磁场不是常量，而是由励磁电流建立并受转子电阻衰减。  
因此 `\psi_\text{AF}` 不再是 `i_D` 的代数量，而是一个独立的一阶动态状态。

其动力学为

$$
\boxed{
\frac{d\psi_\text{AF}}{dt}
=
R_\text{req}i_D
-\frac{R_\text{req}}{L_d-L_q}\psi_\text{AF}
}
$$

物理含义：

- `R_\text{req}i_D`：d 轴励磁电流建立 active flux。
- `-\dfrac{R_\text{req}}{L_d-L_q}\psi_\text{AF}`：由于等效转子电阻，磁链会衰减。

于是 IM 的 `d` 轴电流方程写成

$$
\boxed{
\frac{di_D}{dt}
=
\frac{u_D-R_s i_D+\omega_\text{syn}L_q i_Q-\frac{d\psi_\text{AF}}{dt}}{L_q}
}
$$

其中 `d\psi_\text{AF}/dt` 由上面的独立状态方程给出。

---

### Task 4.3 机械方程

#### (a) 转矩方程

统一转矩表达式为

$$
\boxed{
T_\text{em}=\frac{3}{2}n_\text{pp}\psi_\text{AF}i_Q
}
$$

这说明在统一 active-flux 框架下，所有电机的转矩形式都相同，差别只体现在 `\psi_\text{AF}` 的来源。

#### (b) 牛顿第二定律的旋转形式

机械方程为

$$
J_s\frac{d\omega_\text{mech}}{dt}=T_\text{em}-T_\text{load},
$$

即

$$
\boxed{
\frac{d\omega_\text{mech}}{dt}
=
\frac{T_\text{em}-T_\text{load}}{J_s}
}
$$

这就是代码中的 `\dot x[1]`。

#### (c) 转角方程与滑差项

机械角位置状态满足

$$
\boxed{
\frac{d\theta_\text{mech}}{dt}
=
\omega_\text{mech}
+\frac{\omega_\text{slip}}{n_\text{pp}}
}
$$

其中

$$
\omega_\text{slip}
=
\frac{R_\text{req} i_Q}{\psi_\text{AF}}.
$$

对 PMSM 而言：

- `R_\text{req}=0`，转子不是通过感应电流建立磁场，而是由永磁体直接提供磁链。
- 因此不需要滑差去“感应出”转子磁场。

所以

$$
\omega_\text{slip}=0.
$$

于是 PMSM 有

$$
\frac{d\theta_\text{mech}}{dt}=\omega_\text{mech}.
$$

而 IM 则必须存在非零滑差，才能在转子回路中感应电流并建立磁场。

---

## 五条统一 ODE 总结

设状态向量为

$$
\mathbf x=
\begin{bmatrix}
\theta_\text{mech}\\
\omega_\text{mech}\\
\psi_\text{AF}\\
i_D\\
i_Q
\end{bmatrix},
$$

则统一模型为

$$
\dot x[0]
=
\omega_\text{mech}
+\frac{\omega_\text{slip}}{n_\text{pp}},
$$

$$
\dot x[1]
=
\frac{T_\text{em}-T_\text{load}}{J_s},
$$

$$
T_\text{em}
=
\frac{3}{2}n_\text{pp}\psi_\text{AF}i_Q,
$$

$$
\omega_\text{syn}=n_\text{pp}\omega_\text{mech}+\omega_\text{slip}.
$$

其中

$$
\omega_\text{slip}
=
\begin{cases}
\dfrac{R_\text{req}i_Q}{\psi_\text{AF}}, & R_\text{req}>0 \quad (\text{IM})\\[6pt]
0, & R_\text{req}=0 \quad (\text{PMSM/IPMSM/SynRM})
\end{cases}
$$

### 1. IM 分支 `R_\text{req}>0`

$$
\dot x[2]
=
\frac{d\psi_\text{AF}}{dt}
=
R_\text{req}i_D-\frac{R_\text{req}}{L_d-L_q}\psi_\text{AF},
$$

$$
\dot x[3]
=
\frac{di_D}{dt}
=
\frac{
u_D-R_s i_D+\omega_\text{syn}L_q i_Q-\dot x[2]
}{L_q},
$$

$$
\dot x[4]
=
\frac{di_Q}{dt}
=
\frac{
u_Q-R_s i_Q-\omega_\text{syn}L_q i_D-\omega_\text{syn}\psi_\text{AF}
}{L_q}.
$$

### 2. PMSM / IPMSM / SynRM 分支 `R_\text{req}=0`

$$
\dot x[3]
=
\frac{di_D}{dt}
=
\frac{
u_D-R_s i_D+\omega_\text{syn}L_q i_Q
}{L_d},
$$

$$
\dot x[2]
=
\frac{d\psi_\text{AF}}{dt}
=(L_d-L_q)\dot x[3],
$$

$$
\dot x[4]
=
\frac{di_Q}{dt}
=
\frac{
u_Q-R_s i_Q-\omega_\text{syn}L_q i_D-\omega_\text{syn}\psi_\text{AF}
}{L_q}.
$$

这正是课程中 `DYNAMICS_MACHINE` 的数学结构。

---

## 详细评分细则

总分 `100` 分，与题面分配一致：

- Part 1: `30` 分
- Part 2: `20` 分
- Part 3: `30` 分
- Part 4: `20` 分

评分原则：

- 公式正确但缺少解释，只给该小问的 `60%` 到 `80%`。
- 只有结论没有推导，通常不给满分。
- 符号定义混乱、坐标系方向未说明、正负号错误，会按严重程度扣分。
- 图题要求“画图”的，若完全不画图，相关小问最多给一半分。

### Part 1 评分细则（30 分）

| 小问 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| 1.1(a) | 4 | 正确写出 `p_\alpha=r\cos\omega t, p_\beta=r\sin\omega t`，或写出含初相位的一般式 | 漏一个分量扣 2 分；把正弦余弦写反但能自洽扣 1 分；无时间变量扣 1 分 |
| 1.1(b) | 2 | 明确投影为 `p_\alpha(t)`，并指出是余弦型简谐运动 | 只写“简谐”不给具体表达式扣 1 分；把投影写成 `p_\beta` 扣 1 分 |
| 1.1(c) | 3 | 写出 `dq` 系坐标为 `[r,0]^T`，并说明比静止系更简单 | 只写“是常数”不写具体坐标扣 1 分；未比较简化意义扣 1 分 |
| 1.2(a) | 4 | 正确写出旋转矩阵并验证 `R(-\theta)=R^T(\theta)` | 矩阵符号错一处扣 1 到 2 分；只写结论不验证扣 1 分 |
| 1.2(b) | 3 | 完成矩阵乘法并得到 `[r,0]^T` | 只写“可验证”无计算扣 2 分；结果对但过程缺失扣 1 分 |
| 1.3(a) | 4 | 写出 `x_abs=X+Y`，推得 `a_abs=\ddot X=\dot v`；说明匀速时为 0，加速时出现惯性项 | 只写“有惯性力”但没算加速度扣 2 分；忽略“匀速时为 0”扣 1 分 |
| 1.3(b) | 4 | 正确写出 `d/dt|_inertial = d/dt|_rotating + j\omega` 或等价向量式 | 少写 `j\omega p` 扣 2 分；只写文字不写公式扣 2 分 |
| 1.4(a) | 3 | 清楚说明“储能变量变化引起反作用 effort”的类比，并给出合适对应关系 | 只说“都和导数有关”扣 1 到 2 分；类比对象完全错位扣 2 分 |
| 1.4(b) | 3 | 写出旋转系反电动势表达式，并指出 `j\omega\psi` 是 speed EMF | 少掉额外项扣 2 分；只说“与速度有关”但不解释物理意义扣 1 分 |

### Part 2 评分细则（20 分）

| 小问 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| 2.1(a) | 6 | 明确区分 `\vec B` 是空间分布场量、`\psi` 是适合时域 ODE 的电路 lumped quantity，并联系 `u=d\psi/dt` | 只说“磁链更方便”扣 2 到 3 分；未提时域建模扣 2 分 |
| 2.1(b) | 6 | 正确推得 `u=L\,di/dt` | 直接写答案但无推导扣 1 分；把 `L` 当变量扣 2 分 |
| 2.1(c) | 8 | 完整展开 `u=L di/dt + i dL/dt`，并指出第二项对应速度相关反电动势 | 只写出两项但不说明哪项是 back-EMF 扣 2 分；乘积求导错误扣 3 到 4 分 |

### Part 3 评分细则（30 分）

| 小问 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| 3.1(a) | 5 | 画出单相等效电路，标出 `u, R_s, L_\sigma, L_m` 或等效 `d\psi_m/dt`，并标出电流/磁链方向 | 不画图最多给 2 分；漏标两个以上元件扣 2 到 3 分 |
| 3.1(b) | 4 | 正确写出 `u=R_si+d\psi_\sigma/dt+d\psi_m/dt` 并说明每项物理意义 | 只写 KVL 不解释扣 1 分；漏掉主磁链项或漏感项扣 2 分 |
| 3.1(c) | 3 | 正确区分漏磁路和主磁路，指出谁与转子交链 | 只给结论不解释扣 1 分；两者颠倒扣 2 分 |
| 3.2(a) | 4 | 写出 `i_A+i_B+i_C=0` | 写成非零常数扣 4 分；少一相扣 2 分 |
| 3.2(b) | 6 | 明确说明三相电流受 KCL 约束落在二维平面内，因此可用 `\alpha\beta` 二维矢量无损表示；有 Clarke 变换更好 | 只说“可以降维”但未解释原因扣 2 到 3 分 |
| 3.2(c) | 8 | 准确联系 Part 1 的“旋转矢量变常数”思想，说明 `dq` 坐标如何把正弦量变成直流量并简化 ODE | 只说“更简单”但未说明为什么扣 3 分；未联系圆周运动例子扣 2 分 |

### Part 4 评分细则（20 分）

| 小问 | 分值 | 满分标准 | 常见扣分点 |
|---|---:|---|---|
| 4.1(a) | 2 | 写出 `u=R_si+d\psi/dt` | 漏掉电阻项或导数项扣 1 到 2 分 |
| 4.1(b) | 3 | 正确写出 `dq` 电压方程，至少应包含 `d` 轴耦合项 `-\omega L_q i_Q` 和 `q` 轴项 `+\omega(\psi_AF+L_q i_D)` | 漏任一关键耦合项扣 1 到 2 分；正负号错误视严重程度扣 1 到 2 分 |
| 4.1(c) | 3 | 正确整理出 `di_D/dt`、`di_Q/dt` 的状态方程 | `di_D/dt` 漏 `d\psi_AF/dt` 扣 2 分；`di_Q/dt` 漏 `\omega L_q i_D` 扣 1 分 |
| 4.2(a) | 3 | PMSM 下正确写出 `d\psi_AF/dt=(L_d-L_q)di_D/dt`，并能推到 `di_D/dt` 分母变为 `L_d` | 只写前式不说明 `L_d` 来源扣 1 分 |
| 4.2(b) | 3 | IM 下正确写出 `d\psi_AF/dt=R_req i_D-\frac{R_req}{L_d-L_q}\psi_AF`，并说明其是一阶独立动态 | 公式对但不解释物理意义扣 1 分；把 `i_D` 写成 `i_Q` 扣 2 分 |
| 4.3(a) | 2 | 写出 `T_em=\frac32 n_pp \psi_AF i_Q` | 漏 `3/2` 或漏 `n_pp` 扣 1 分 |
| 4.3(b) | 2 | 写出 `d\omega_mech/dt=(T_em-T_load)/J_s` | 分母写错或正负号错扣 1 到 2 分 |
| 4.3(c) | 2 | 写出 `d\theta_mech/dt=\omega_mech+\omega_slip/n_pp`，并说明 PMSM 因 `R_req=0` 故滑差为 0 | 只写结论不解释 PMSM 原因扣 1 分 |

---

## 建议的阅卷补充规则

### 一、过程分建议

- 若最终公式略有符号错误，但推导路线清晰，可保留 `50%` 以上过程分。
- 若最终答案正确但明显抄结论、缺少中间关键步骤，建议扣除该问 `20%` 到 `40%`。
- 要求英文回答，中文回答零分且提示。

### 二、图示分建议

- 题目要求“画图”的 Task 3.1(a)，若图中元件齐全、方向清楚、标注规范，可给满分。
- 若只写文字不画电路图，最多给该小问 `40%`。

### 三、符号与单位规范

- `\omega_\text{mech}`、`\omega_\text{syn}`、`\omega_\text{slip}` 混用，每处扣 `0.5` 到 `1` 分。
- `L_d`、`L_q`、`\psi_\text{AF}`、`\psi_\text{PM}` 含义不清，每个关键符号最多扣 `1` 分。
- 若全文完全没有说明状态变量对应关系 `x[0] \sim x[4]`，Part 4 可酌情再扣 `1` 分。

### 四、优秀答案特征

优秀答案通常同时满足以下几点：

- 不仅给出结论，还解释“为什么做坐标变换”。
- 能清楚地区分“绝对导数”和“相对导数”。
- 能把 `\psi_\text{AF}` 看成统一建模的核心状态，而不是孤立公式。
- 能把 Part 1 的几何直觉一路连到 Part 4 的五条 ODE。

---

## 最终可提交版结论

这份作业的核心结论可以压缩成一句话：

> 三相电机在 `\alpha\beta` 平面中本质上是一个旋转矢量系统；经过 `dq` 坐标变换、修正求导算子 `d/dt + j\omega`、磁链建模和 active-flux 重写后，就能把 IM、PMSM、IPMSM、SynRM 统一为同一组五维状态微分方程。

其中五个状态分别是：

$$
\theta_\text{mech},\quad
\omega_\text{mech},\quad
\psi_\text{AF},\quad
i_D,\quad
i_Q.
$$

这就是课程代码 `DYNAMICS_MACHINE` 的完整数学来源。

# 结果输出

1. 最终分数提交写入course_grades_template，注意写入的编码
2. 给每人的作业都逐题写错题批改评论。
