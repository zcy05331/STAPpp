# STAPpp Q4 平面四边形单元课程设计报告

## 1. 引言

本课程设计在 STAPpp 教学有限元程序中新增四节点双线性等参四边形单元 Q4，用于二维平面应力和平面应变线弹性静力分析。新增功能包括：

- Q4 单元刚度矩阵计算；
- 平面应力/平面应变本构矩阵；
- 四点 Gauss 积分；
- 四个 Gauss 点处的 $\sigma_x, \sigma_y, \tau_{xy}$ 应力输出；
- Q4 输入格式、材料参数和单元组注册；
- 单元分片试验、收敛性分析、独立解析解验证；
- 可视化结果（预留，待后续补充）。

实现保持 STAPpp 原有框架：不重构 `Node`、`Domain` 或 `Solver`，仅通过新增单元类和扩展已有单元组/输出分派完成开发。

## 2. Q4 单元基本原理

### 2.1 等参映射与形函数

Q4 单元在自然坐标 $(\xi, \eta) \in [-1,1]^2$ 中采用双线性形函数：

$$
N_1 = \frac{1}{4}(1-\xi)(1-\eta), \quad
N_2 = \frac{1}{4}(1+\xi)(1-\eta),
$$

$$
N_3 = \frac{1}{4}(1+\xi)(1+\eta), \quad
N_4 = \frac{1}{4}(1-\xi)(1+\eta).
$$

几何映射和位移插值均使用同一组形函数：

$$
x = \sum_{i=1}^{4} N_i x_i, \quad y = \sum_{i=1}^{4} N_i y_i,
$$

$$
u = \sum_{i=1}^{4} N_i u_i, \quad v = \sum_{i=1}^{4} N_i v_i.
$$

Jacobian 矩阵为：

$$
\mathbf{J} =
\begin{bmatrix}
\partial x / \partial \xi & \partial x / \partial \eta \\
\partial y / \partial \xi & \partial y / \partial \eta
\end{bmatrix}.
$$

程序在每个 Gauss 点检查 $\det \mathbf{J} > 0$。若节点顺序错误或单元退化导致 $\det \mathbf{J} \le 0$，程序立即报错终止。

### 2.2 应变-位移矩阵

平面问题中应变向量为：

$$
\boldsymbol{\varepsilon} =
\begin{bmatrix}
\varepsilon_x & \varepsilon_y & \gamma_{xy}
\end{bmatrix}^{\mathrm{T}}
= \mathbf{B}\mathbf{d}_e.
$$

对四节点 Q4 单元，单元自由度排列为：

$$
\mathbf{d}_e =
\begin{bmatrix}
u_1 & v_1 & u_2 & v_2 & u_3 & v_3 & u_4 & v_4
\end{bmatrix}^{\mathrm{T}}.
$$

矩阵 $\mathbf{B}$ 的结构为：

$$
\mathbf{B}_i =
\begin{bmatrix}
\partial N_i/\partial x & 0 \\
0 & \partial N_i/\partial y \\
\partial N_i/\partial y & \partial N_i/\partial x
\end{bmatrix}.
$$

### 2.3 本构矩阵

平面应力：

$$
\mathbf{D} = \frac{E}{1-\nu^2}
\begin{bmatrix}
1 & \nu & 0 \\
\nu & 1 & 0 \\
0 & 0 & (1-\nu)/2
\end{bmatrix}.
$$

平面应变：

$$
\mathbf{D} = \frac{E}{(1+\nu)(1-2\nu)}
\begin{bmatrix}
1-\nu & \nu & 0 \\
\nu & 1-\nu & 0 \\
0 & 0 & (1-2\nu)/2
\end{bmatrix}.
$$

其中 $E$ 为 Young 模量，$\nu$ 为 Poisson 比，$t$ 为单元厚度。

### 2.4 单元刚度矩阵

Q4 单元刚度矩阵为：

$$
\mathbf{K}_e = \int_{-1}^{1}\int_{-1}^{1}
\mathbf{B}^{\mathrm{T}}\mathbf{D}\mathbf{B} \, t \, \det\mathbf{J} \, d\xi d\eta.
$$

程序采用 $2\times2$ Gauss 积分，四个积分点依次为：

$$
(-1/\sqrt{3}, -1/\sqrt{3}),
( 1/\sqrt{3}, -1/\sqrt{3}),
( 1/\sqrt{3},  1/\sqrt{3}),
(-1/\sqrt{3},  1/\sqrt{3}).
$$

## 3. 程序实现方案

### 3.1 新增和修改文件

核心实现文件：

- `src/h/Q4.h`
- `src/cpp/Q4.cpp`
- `src/h/Material.h`
- `src/cpp/Material.cpp`
- `src/h/ElementGroup.h`
- `src/cpp/ElementGroup.cpp`
- `src/h/Outputter.h`
- `src/cpp/Outputter.cpp`

辅助验证文件：

- `make/validate_q4_cases.py`
- `data/q4_patch_single/`
- `data/q4_patch_multi/`
- `data/q4_convergence/`
- `data/q4_validation/`
- `data/q4_invalid/`

### 3.2 与原 STAPpp 架构的衔接

STAPpp 当前节点自由度数为 `CNode::NDF = 3`，即每个节点有 $x,y,z$ 三个平动自由度。Q4 平面单元只使用 $x,y$ 两个自由度。为避免未使用的 $z$ 自由度进入总方程形成零刚度行列，本实现采用以下约束：

1. Q4 输入文件中所有 Q4 节点必须设置 `bz=1`；
2. `CQ4::Read()` 检查所有关联节点的 `bcode[2]`，若存在活动 $z$ 自由度则报错；
3. `CQ4::GenerateLocationMatrix()` 覆盖基类默认实现，仅生成 8 个平面自由度。

Q4 单元定位矩阵顺序为：

```text
[n1 ux, n1 uy, n2 ux, n2 uy, n3 ux, n3 uy, n4 ux, n4 uy]
```

单元刚度矩阵按 STAPpp skyline 组装接口要求存储为上三角列优先格式。

### 3.3 应力输出约定

Q4 的 `ElementStress()` 输出 12 个浮点数，即 4 个 Gauss 点乘以每点 3 个应力分量：

```text
GP1: sigma_x sigma_y tau_xy
GP2: sigma_x sigma_y tau_xy
GP3: sigma_x sigma_y tau_xy
GP4: sigma_x sigma_y tau_xy
```

输出表头为：

```text
ELEMENT GP SIGMA_X SIGMA_Y TAU_XY
```

## 4. 输入数据格式

Q4 单元组控制行：

```text
2 NUME NUMMAT
```

其中 `2` 表示 `ElementTypes::Q4`。

Q4 材料行：

```text
nset E nu thickness analysisType
```

- `analysisType = 0`：平面应力；
- `analysisType = 1`：平面应变。

Q4 单元行：

```text
N n1 n2 n3 n4 mset
```

节点应按逆时针物理顺序输入。节点行仍采用 STAPpp 原格式：

```text
node bx by bz x y z
```

Q4 节点必须令 `bz=1`。

## 5. 分片试验

### 5.1 单单元分片试验

算例目录：`data/q4_patch_single/`

模型为单位正方形单个 Q4 单元，左边界约束 $u=0$，左下角约束 $v=0$，右边界施加等效节点力，使解析解为单向常应力：

$$
\sigma_x = 1, \quad \sigma_y = 0, \quad \tau_{xy} = 0.
$$

解析位移场为：

$$
u = \frac{\sigma_x}{E}x, \quad v = -\frac{\nu\sigma_x}{E}y.
$$

验证结果：

| 项目 | 最大误差 |
|---|---:|
| Gauss 点应力误差 | $1.429600\times10^{-16}$ |
| 节点位移误差 | $2.020440\times10^{-19}$ |

误差远小于验收阈值 $10^{-8}$，单单元分片试验通过。

### 5.2 多单元分片试验

算例目录：`data/q4_patch_multi/`

模型为 $2\times2$ Q4 网格，载荷和边界条件与单单元分片试验一致。

验证结果：

| 项目 | 最大误差 |
|---|---:|
| Gauss 点应力误差 | $3.728310\times10^{-16}$ |
| 节点位移误差 | $7.278450\times10^{-20}$ |

误差远小于验收阈值 $10^{-6}$，多单元分片试验通过。

## 6. 收敛性分析

算例目录：`data/q4_convergence/`

收敛性分析采用平面应力悬臂矩形膜问题。左边界全约束，右边界施加竖向总力 $P=-1$ 的等效节点力。采用 $32\times8$ 网格作为参考解，对较粗网格的右上角竖向位移进行比较。

| nx | ny | 右上角 $u_y$ | 相对误差 | 观测阶 |
|---:|---:|---:|---:|---:|
| 2 | 1 | $-1.0000000000\times10^{-1}$ | $1.650270\times10^{-1}$ | - |
| 4 | 2 | $-1.8661600000\times10^{-1}$ | $7.841100\times10^{-2}$ | 1.0736 |
| 8 | 2 | $-2.3717000000\times10^{-1}$ | $2.785700\times10^{-2}$ | 1.4930 |
| 16 | 4 | $-2.5867900000\times10^{-1}$ | $6.348000\times10^{-3}$ | 2.1337 |

误差随网格加密单调下降，说明 Q4 单元实现具有合理的网格收敛行为。由于本算例采用边界等效节点力和有限网格参考解，表中观测阶用于说明趋势，不将其过度解释为严格理论收敛阶。

## 7. 验证算例

算例目录：`data/q4_validation/`

验证算例采用 $4\times2$ Q4 网格的单向拉伸矩形膜，并与独立连续体解析解比较：

$$
\sigma_x = 1, \quad \sigma_y = 0, \quad \tau_{xy} = 0,
$$

$$
u = \frac{x}{E}, \quad v = -\frac{\nu y}{E}.
$$

验证结果：

| 项目 | 最大误差 | 验收标准 |
|---|---:|---:|
| Gauss 点应力误差 | $3.682050\times10^{-16}$ | $\le 1\%$ |
| 节点位移误差 | $1.730770\times10^{-19}$ | $\le 1\%$ |

验证算例通过。

## 8. 后处理与可视化（预留）

（本节预留空白，用于后续补充 TecPlot/ParaView 变形图、应力云图、图号、图注及结果说明。）

## 9. 构建与验证命令

由于仓库中已有 `build/` 缓存较旧，本项目采用临时目录进行新构建：

```bash
tmpbuild=$(mktemp -d /tmp/stappp-build-XXXXXX)
cmake -S src -B "$tmpbuild" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$tmpbuild" -- -j2
python3 make/validate_q4_cases.py --exe "$tmpbuild/stap++"
```

本报告中的数值结果由上述验证脚本生成。

## 10. 结论

本课程设计在 STAPpp 中完成了 Q4 平面四边形单元扩展，并保持原有 Bar 单元行为不变。Q4 单元通过了单单元分片试验、多单元分片试验、网格收敛性分析和独立解析解验证。可视化结果部分已预留，待后续补充变形图和应力云图。

后续可扩展方向包括：

- 增加体力项和边界分布力的原生输入；
- 增加更通用的二维/三维混合自由度管理；
- 实现 T3、H8、Beam、Plate 等更多单元；
- 增加自动化网格生成和更完整的 ABAQUS 对比算例。

## 参考文献

1. STAPpp 原始程序与课程说明文档。
2. 张雄、王天舒、刘岩，《计算动力学（第2版）》，清华大学出版社。
3. 有限元法基础课程讲义与 `Course Project - 2026.pdf`。
