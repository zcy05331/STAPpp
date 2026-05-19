# STAPpp 有限元求解程序代码精读总结

> **程序名称**: STAP++ (Structural Analysis Program in C++)  
> **版本**: Release 1.11, November 22, 2017  
> **开发单位**: 清华大学航天航空学院计算动力学实验室  
> **课程配套**: 有限元方法（张雄教授）  
> **输入输出格式**: 与 STAP90 (Fortran 90) 完全兼容

---

## 1. 项目概览

STAPpp 是一个用 C++ 实现的有限元方法（FEM）教学代码，采用面向对象设计，用于求解**线性静力学**问题。当前版本仅实现了**三维杆（Bar/Truss）单元**，但预留了 Q4、T3、H8、Beam、Plate、Shell 等单元类型的扩展接口。

### 1.1 文件结构

| 目录       | 内容                       |
| ---------- | -------------------------- |
| `src/h/`   | 头文件（类声明与模板实现） |
| `src/cpp/` | 源文件（类方法实现）       |
| `data/`    | 示例输入数据文件（`.dat`） |
| `docs/`    | Doxygen 生成的文档         |

### 1.2 源文件清单

| 文件                 | 职责                                     |
| -------------------- | ---------------------------------------- |
| `main.cpp`           | 主程序入口，控制整个求解流程             |
| `Domain.h/cpp`       | 问题域类（单例），管理全部数据与求解步骤 |
| `Node.h/cpp`         | 节点类：坐标、边界条件、方程号           |
| `Element.h`          | 单元基类（纯虚类），定义单元通用接口     |
| `Bar.h/cpp`          | 杆单元派生类，实现刚度矩阵和应力计算     |
| `Material.h/cpp`     | 材料基类与杆材料派生类                   |
| `ElementGroup.h/cpp` | 单元组：管理同类型单元的集合             |
| `LoadCaseData.h/cpp` | 荷载工况数据                             |
| `SkylineMatrix.h`    | Skyline 存储格式的稀疏矩阵模板类         |
| `Solver.h/cpp`       | LDLT 分解求解器                          |
| `Outputter.h/cpp`    | 输出管理类（单例），同时输出到文件和终端 |
| `Clock.h/cpp`        | 计时工具类                               |
| `CMakeLists.txt`     | CMake 构建配置                           |

---

## 2. 整体求解流程

`main.cpp` 控制的完整有限元分析流程如下：

```
┌─────────────────────────────────────────────┐
│  1. 读取输入数据 (FEMData->ReadData)         │
│     ├─ 读标题行、控制参数                     │
│     ├─ 读节点数据 → 计算方程编号               │
│     ├─ 读荷载数据                             │
│     └─ 读单元数据（含材料属性）                 │
├─────────────────────────────────────────────┤
│  2. 分配矩阵存储 (FEMData->AllocateMatrices) │
│     ├─ 分配力/位移向量                        │
│     ├─ 计算列高 (Column Heights)              │
│     ├─ 计算对角元地址                          │
│     └─ 分配 Skyline 刚度矩阵存储               │
├─────────────────────────────────────────────┤
│  3. 组装总刚度矩阵                            │
│     (FEMData->AssembleStiffnessMatrix)        │
├─────────────────────────────────────────────┤
│  4. LDLT 分解 (Solver->LDLT)                 │
├─────────────────────────────────────────────┤
│  5. 对每个荷载工况循环:                        │
│     ├─ 组装荷载向量 (AssembleForce)            │
│     ├─ 回代求解位移 (BackSubstitution)         │
│     ├─ 输出节点位移                            │
│     └─ 计算并输出单元应力                      │
├─────────────────────────────────────────────┤
│  6. 输出求解时间统计                           │
└─────────────────────────────────────────────┘
```

**关键设计**: 刚度矩阵只分解一次，多个荷载工况共用同一分解结果，仅需重复回代过程，大幅提高多工况分析效率。

---

## 3. 核心数据结构

### 3.1 CDomain —— 问题域（单例模式）

```cpp
class CDomain {
    static CDomain* _instance;   // 单例指针
    
    char Title[256];             // 问题标题
    unsigned int MODEX;          // 0=仅检查数据, 1=执行计算
    unsigned int NUMNP;          // 总节点数
    CNode* NodeList;             // 节点数组
    unsigned int NUMEG;          // 单元组数
    CElementGroup* EleGrpList;   // 单元组数组
    unsigned int NLCASE;         // 荷载工况数
    CLoadCaseData* LoadCases;    // 荷载工况数组
    unsigned int NEQ;            // 总方程数（自由度数）
    CSkylineMatrix<double>* StiffnessMatrix;  // 总刚度矩阵
    double* Force;               // 全局力/位移向量（复用同一数组）
};
```

**设计要点**:
- 采用**单例模式**确保全局仅一个 `CDomain` 实例
- `Force` 和 `Displacement` 指向同一数组——求解后力向量被位移值覆盖
- `GetDisplacement()` 和 `GetForce()` 返回同一指针

### 3.2 CNode —— 节点

```cpp
class CNode {
    static const unsigned int NDF = 3;  // 每节点自由度数（3D: x, y, z）
    unsigned int NodeNumber;            // 节点编号
    double XYZ[3];                      // 坐标 (x, y, z)
    unsigned int bcode[NDF];            // 边界条件码 / 方程号
};
```

**边界条件与方程编号的双重用途**:

| 阶段                           | `bcode[i]` 含义                                |
| ------------------------------ | ---------------------------------------------- |
| 读入阶段                       | `0` = 自由（活跃自由度）; `1` = 约束（非活跃） |
| `CalculateEquationNumber()` 后 | `0` = 被约束（无方程）; `>0` = 全局方程编号    |

编号逻辑：遍历所有节点的所有自由度，若 `bcode=1`（约束）则置 0；若 `bcode=0`（自由）则分配递增方程号。注意输入中的含义与内部表示是**取反**的。

### 3.3 CElement —— 单元基类

```cpp
class CElement {
    unsigned int NEN_;              // 单元节点数
    CNode** nodes_;                 // 节点指针数组
    CMaterial* ElementMaterial_;    // 材料指针
    unsigned int* LocationMatrix_;  // 定位矩阵（自由度 → 全局方程号）
    unsigned int ND_;               // 定位矩阵维度 = NEN_ × NDF

    // 纯虚函数 —— 必须由派生类实现
    virtual bool Read(...) = 0;
    virtual void Write(...) = 0;
    virtual void ElementStiffness(double* stiffness) = 0;
    virtual void ElementStress(double* stress, double* Displacement) = 0;
    
    // 通用实现
    virtual void GenerateLocationMatrix();     // 生成定位矩阵
    virtual unsigned int SizeOfStiffnessMatrix(); // 上三角存储大小 = ND*(ND+1)/2
};
```

**定位矩阵**（Location Matrix）是单元层面到全局层面的关键映射：
- 大小为 `ND_ = NEN_ × NDF`
- `LocationMatrix_[i]` 存储单元第 `i` 个自由度对应的全局方程号
- 值为 0 表示该自由度被约束

### 3.4 CBar —— 三维杆单元

杆单元是目前唯一实现的单元类型：
- **节点数** `NEN_ = 2`，每个节点 3 个自由度
- **定位矩阵维度** `ND_ = 6`
- **单元刚度矩阵大小** = 6×7/2 = 21 个元素（上三角，按列存储）

### 3.5 CSkylineMatrix —— Skyline 稀疏矩阵

```cpp
template <class T_>
class CSkylineMatrix {
    T_* data_;                      // 一维数组存储 skyline 数据
    unsigned int NEQ_;              // 矩阵维度
    unsigned int MK_;               // 最大半带宽
    unsigned int NWK_;              // 存储总元素数
    unsigned int* ColumnHeights_;   // 各列列高数组
    unsigned int* DiagonalAddress_; // 对角元在 data_ 中的地址（从1编号）
};
```

**Skyline 存储原理**:

对于对称正定刚度矩阵，只存储上三角中从 skyline（每列第一个非零元）到对角元的部分：

```
对于第 j 列: 存储从第 (j - ColumnHeight[j]) 行 到第 j 行 的元素
对角元地址: DiagonalAddress[0] = 1
            DiagonalAddress[j] = DiagonalAddress[j-1] + ColumnHeight[j-1] + 1
```

**关键操作**:

| 方法                         | 功能                                 |
| ---------------------------- | ------------------------------------ |
| `CalculateColumnHeight()`    | 根据单元定位矩阵计算各列列高         |
| `CalculateDiagnoalAddress()` | 计算对角元在一维数组中的地址         |
| `Allocate()`                 | 根据对角元地址计算总存储量并分配内存 |
| `Assembly()`                 | 将单元刚度矩阵组装到总刚度矩阵       |
| `operator()(i,j)`            | 通过行列号（从1编号）访问矩阵元素    |

---

## 4. 核心算法详解

### 4.1 方程编号与列高计算

**方程编号**（`CDomain::CalculateEquationNumber`）：

```
NEQ = 0
对每个节点 np:
    对每个自由度 dof:
        若 bcode[dof] == 1 (约束):  bcode[dof] = 0
        若 bcode[dof] == 0 (自由):  NEQ++; bcode[dof] = NEQ
```

**列高计算**（`CSkylineMatrix::CalculateColumnHeight`）：

对每个单元的定位矩阵：
1. 找到所有活跃自由度中最小的全局方程号 `nfirstrow`
2. 对每个活跃自由度的列号 `column`，列高至少为 `column - nfirstrow`
3. 取所有单元贡献中的最大值作为最终列高

这保证了 Skyline 存储能覆盖所有非零元素。

### 4.2 单元刚度矩阵组装

**杆单元刚度矩阵**（`CBar::ElementStiffness`）：

三维杆单元的 $6 \times 6$ 刚度矩阵为：

$$\mathbf{K}^e = \frac{EA}{L^3} \begin{bmatrix} \mathbf{D} & -\mathbf{D} \\ -\mathbf{D} & \mathbf{D} \end{bmatrix}$$

其中 $\mathbf{D}$ 为 $3\times3$ 方向余弦矩阵：

$$\mathbf{D} = \begin{bmatrix} dx^2 & dx \cdot dy & dx \cdot dz \\ dx \cdot dy & dy^2 & dy \cdot dz \\ dx \cdot dz & dy \cdot dz & dz^2 \end{bmatrix}$$

$dx = x_2 - x_1$, $dy = y_2 - y_1$, $dz = z_2 - z_1$, $L$ 为杆长。

代码中用 `k = EA / (L × L²)` 作为系数，矩阵按**列优先上三角**格式存储在长度为 21 的一维数组中。

**全局组装**（`CSkylineMatrix::Assembly`）：

```cpp
for j = 0 to ND-1:
    Lj = LocationMatrix[j]   // 第j个DOF的全局方程号
    if Lj == 0: skip         // 约束自由度
    DiagjElement = j*(j+1)/2 // 单元刚度矩阵中第j列对角元地址
    for i = 0 to j:
        Li = LocationMatrix[i]
        if Li == 0: skip
        K(Li, Lj) += Matrix[DiagjElement + j - i]
```

### 4.3 LDLT 分解求解器

**LDLT 分解**（`CLDLTSolver::LDLT`）：

将对称正定矩阵 $\mathbf{K}$ 分解为 $\mathbf{K} = \mathbf{L}\mathbf{D}\mathbf{L}^T$，其中 $\mathbf{L}$ 为下三角矩阵（对角元为 1），$\mathbf{D}$ 为对角矩阵。

分解过程中利用 Skyline 存储的特点，只对每列 skyline 范围内的元素进行运算：

```
对第 j 列 (j = 2, ..., N):
    mj = j - ColumnHeight[j]   // 第j列第一个非零元的行号
    
    // 第一步: 计算 U_ij (i = mj+1, ..., j-1)
    对 i = mj+1 到 j-1:
        mi = i - ColumnHeight[i]
        C = Σ(r = max(mi,mj) 到 i-1) L_ri × U_rj
        U_ij = K_ij - C
    
    // 第二步: 计算 L_rj 和 D_jj
    对 r = mj 到 j-1:
        L_rj = U_rj / D_rr
        D_jj -= L_rj × U_rj
    
    // 检查正定性
    若 |D_jj| ≤ ε: 报错退出
```

**回代求解**（`CLDLTSolver::BackSubstitution`）：

求解 $\mathbf{L}\mathbf{D}\mathbf{L}^T \mathbf{a} = \mathbf{R}$ 分三步完成，结果直接覆写力向量：

| 步骤     | 操作                                           | 方程                                                  |
| -------- | ---------------------------------------------- | ----------------------------------------------------- |
| 前代     | $\mathbf{L}\mathbf{V} = \mathbf{R}$            | $V_i = R_i - \sum_{j=m_i}^{i-1} L_{ji} V_j$           |
| 对角标定 | $\bar{\mathbf{V}} = \mathbf{D}^{-1}\mathbf{V}$ | $\bar{V}_i = V_i / D_{ii}$                            |
| 回代     | $\mathbf{L}^T\mathbf{a} = \bar{\mathbf{V}}$    | $a_i = \bar{V}_i - \sum_{j=i+1}^{N} L_{ij} \bar{V}_j$ |

### 4.4 应力计算

**杆单元应力**（`CBar::ElementStress`）：

杆的轴向应力为：

$$\sigma = \frac{E}{L^2} \begin{bmatrix} -dx & -dy & -dz & dx & dy & dz \end{bmatrix} \mathbf{d}^e$$

其中 $\mathbf{d}^e$ 是从全局位移向量中通过定位矩阵提取的单元节点位移。

输出时同时给出杆力 $F = \sigma \times A$ 和应力 $\sigma$。

---

## 5. 设计模式与编程技巧

### 5.1 设计模式

| 模式             | 应用                                             | 类                                  |
| ---------------- | ------------------------------------------------ | ----------------------------------- |
| **单例模式**     | 确保全局唯一实例                                 | `CDomain`, `COutputter`             |
| **模板方法模式** | 基类定义求解流程框架，派生类实现具体单元行为     | `CElement` → `CBar`                 |
| **工厂模式**     | `CElementGroup` 根据单元类型创建对应的派生类对象 | `CElementGroup::AllocateElements()` |

### 5.2 Skyline 存储的上三角按列存储

单元刚度矩阵和全局刚度矩阵均采用**上三角、按列存储**的一维数组格式。以 $4\times4$ 矩阵为例：

```
完整矩阵:         上三角按列存储:
[K11 K12 K13 K14]     data[0] = K11
[    K22 K23 K24]     data[1] = K22
[        K33 K34]     data[2] = K12
[            K44]     data[3] = K33
                      data[4] = K23
                      data[5] = K13
                      data[6] = K44
                      data[7] = K34
                      data[8] = K24
                      data[9] = K14
```

第 $j$ 列对角元在一维数组中的位置为 $j(j+1)/2$（从0编号）。

### 5.3 CElementGroup 的内存管理

`CElementGroup` 使用**手动指针算术**来管理派生类数组：

```cpp
// 通过字节偏移实现多态数组的索引访问
CElement& CElementGroup::operator[](unsigned int i) {
    return *(CElement*)((std::size_t)(ElementList_) + i * ElementSize_);
}
```

这里 `ElementSize_` 在运行时根据单元类型设置（如 `sizeof(CBar)`），允许用基类指针正确地定位到派生类对象——因为 `new CBar[size]` 分配的数组中元素按 `sizeof(CBar)` 对齐排列。

### 5.4 输出策略

`COutputter` 重载了 `operator<<`，使得所有输出**同时**写入文件和标准输出：

```cpp
template <typename T>
COutputter& operator<<(const T& item) {
    std::cout << item;     // 屏幕输出
    OutputFile << item;    // 文件输出
    return *this;
}
```

### 5.5 调试模式

通过 CMake 选项 `STAP++_DEBUG` 控制编译宏 `_DEBUG_`，开启后会输出：
- 定位矩阵
- 列高
- 对角元地址
- Skyline 格式与完整格式的刚度矩阵
- 位移向量

---

## 6. 输入数据文件格式

以 `truss.dat` 为例说明输入格式：

```
Cables to test STAP90              ← 标题行
    3    1    1    1                ← NUMNP=3, NUMEG=1, NLCASE=1, MODEX=1
    1    1    1    1  -0.3  0.5196   0.0   ← 节点1: bcode=(1,1,1), xyz
    2    1    1    1  0.5196 0.5196  0.0   ← 节点2: bcode=(1,1,1), xyz
    3    0    0    1   0.0    0.0    0.0   ← 节点3: bcode=(0,0,1), xyz
    1    1                         ← 荷载工况1, 1个集中荷载
    3    2    80.0E3               ← 节点3, 方向2(y), 荷载80kN
    1    2    1                    ← 单元类型=Bar, 2个单元, 1组材料
    1   207.0E9    120E-6          ← 材料组1: E=207GPa, A=120mm²
    1    1    3    1               ← 单元1: 节点1→3, 材料组1
    2    3    2    1               ← 单元2: 节点3→2, 材料组1
```

**数据块顺序**: 标题 → 控制行 → 节点数据 → 荷载数据 → 单元数据（含材料）

---

## 7. 类关系图

```
                    CDomain (单例)
                   /    |    \      \
              CNode  CElementGroup  CLoadCaseData
                        |      \
                   CElement*   CMaterial*
                   (基类)      (基类)
                      |           |
                    CBar      CBarMaterial
                   
               CSkylineMatrix<double>
                      |
                CLDLTSolver (引用)
                
                COutputter (单例)
                
                  Clock (计时)
```

**依赖关系**:
- `CDomain` 拥有并管理所有核心数据结构
- `CElementGroup` 持有 `CElement*` 和 `CMaterial*` 的多态数组
- `CLDLTSolver` 通过引用操作 `CSkylineMatrix`
- `COutputter` 通过 `CDomain::GetInstance()` 获取数据进行输出

---

## 8. 扩展新单元的步骤

STAPpp 的框架设计使得扩展新单元类型较为规范：

1. **定义材料类**: 从 `CMaterial` 派生，添加新材料参数（如惯性矩）
2. **定义单元类**: 从 `CElement` 派生，实现：
   - 构造函数：设置 `NEN_`、`ND_`
   - `Read()`: 从输入流读取单元数据
   - `Write()`: 输出单元信息
   - `ElementStiffness()`: 计算单元刚度矩阵
   - `ElementStress()`: 计算单元应力
3. **注册单元类型**: 在 `ElementTypes` 枚举中添加新类型
4. **更新 `CElementGroup`**: 在 `CalculateMemberSize()`、`AllocateElements()`、`AllocateMaterials()` 的 `switch` 分支中添加新类型的处理
5. **更新 `COutputter`**: 添加新单元的输出格式

### 8.1 添加梁（Beam）单元的详细步骤

#### 需要修改/新增的文件清单

| 文件                   | 操作         | 说明                                  |
| ---------------------- | ------------ | ------------------------------------- |
| `h/Beam.h`             | **新增**     | 梁单元类声明                          |
| `cpp/Beam.cpp`         | **新增**     | 梁单元刚度矩阵与应力实现              |
| `h/Material.h`         | **修改**     | 新增 `CBeamMaterial` 派生类           |
| `cpp/Material.cpp`     | **修改**     | 实现 `CBeamMaterial` 的 Read/Write    |
| `h/ElementGroup.h`     | **修改**     | `#include "Beam.h"`                   |
| `cpp/ElementGroup.cpp` | **修改**     | switch 分支添加 Beam 处理             |
| `h/Outputter.h`        | **修改**     | 声明 `OutputBeamElements()`           |
| `cpp/Outputter.cpp`    | **修改**     | 实现梁单元输出与应力输出              |
| `h/Node.h`             | **可能修改** | 若需要转动自由度，`NDF` 需从 3 改为 6 |

#### Step 1: 修改节点自由度（关键决策点）

三维梁单元每个节点有 **6 个自由度**（3 平动 + 3 转动）。当前 `CNode::NDF = 3`，需要决定处理策略：

**方案 A** — 统一改为 `NDF = 6`（简单但浪费）：
```cpp
// Node.h
const static unsigned int NDF = 6;  // 改为6
```
- 优点：修改最少
- 缺点：纯杆单元节点也会有多余的转动自由度，增大矩阵规模

**方案 B** — 让 `NDF` 可变（推荐但需更多重构）：
```cpp
// Node.h 中将 NDF 从 static const 改为成员变量
unsigned int NDF_;  // 由单元类型决定
```
- 需要修改所有使用 `CNode::NDF` 的地方

#### Step 2: 定义梁材料类

```cpp
// Material.h 中新增
class CBeamMaterial : public CMaterial
{
public:
    double Area;    // 截面面积 A
    double Iy;      // 绕 y 轴惯性矩
    double Iz;      // 绕 z 轴惯性矩
    double J;       // 扭转常数（极惯性矩）
    double nu;      // 泊松比（用于计算剪切模量 G = E/(2(1+ν))）
    
    virtual bool Read(ifstream& Input);
    virtual void Write(COutputter& output);
};
```

#### Step 3: 定义梁单元类

```cpp
// Beam.h
#pragma once
#include "Element.h"

class CBeam : public CElement
{
public:
    CBeam();
    ~CBeam();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
};
```

```cpp
// Beam.cpp 构造函数
CBeam::CBeam()
{
    NEN_ = 2;    // 梁单元: 2 个节点
    nodes_ = new CNode*[NEN_];
    ND_ = 12;    // 每个节点 6 个自由度 → 2×6 = 12
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
}
```

**梁单元刚度矩阵**为 $12 \times 12$，上三角存储大小为 $12 \times 13 / 2 = 78$ 个元素。

Euler-Bernoulli 梁的刚度矩阵（局部坐标系下）包含：
- 轴向: $EA/L$ 项（与杆相同）
- 弯曲: $EI/L^3$ 项（$12EI/L^3$, $6EI/L^2$, $4EI/L$, $2EI/L$）
- 扭转: $GJ/L$ 项

需要通过**坐标变换矩阵** $\mathbf{T}$ 将局部刚度矩阵转换到全局坐标系：$\mathbf{K}^g = \mathbf{T}^T \mathbf{K}^l \mathbf{T}$

#### Step 4: 注册到 ElementGroup

```cpp
// ElementGroup.cpp — CalculateMemberSize() 中添加:
case ElementTypes::Beam:
    ElementSize_ = sizeof(CBeam);
    MaterialSize_ = sizeof(CBeamMaterial);
    break;

// AllocateElements() 中添加:
case ElementTypes::Beam:
    ElementList_ = new CBeam[size];
    break;

// AllocateMaterials() 中添加:
case ElementTypes::Beam:
    MaterialList_ = new CBeamMaterial[size];
    break;
```

#### Step 5: 添加输出

```cpp
// Outputter.h 中声明:
void OutputBeamElements(unsigned int EleGrp);

// Outputter.cpp — OutputElementInfo() 中添加 switch 分支:
case ElementTypes::Beam:
    OutputBeamElements(EleGrp);
    break;

// OutputElementStress() 中添加 Beam 的应力输出分支
```

### 8.2 添加板（Plate）单元

板单元（如 Mindlin 板）与梁的扩展流程类似，但有以下不同：

| 对比项        | 梁单元             | 板单元                    |
| ------------- | ------------------ | ------------------------- |
| 节点数 `NEN_` | 2                  | 4（四节点）或 3（三角形） |
| 节点自由度    | 6 (u,v,w,θx,θy,θz) | 3~5 (w,θx,θy) 或 5~6      |
| `ND_`         | 12                 | 12~24                     |
| 刚度矩阵大小  | 78                 | 78~300                    |
| 需要数值积分  | 否                 | **是**（Gauss 积分）      |
| 坐标变换      | 3D 旋转矩阵        | 等参变换                  |

**板单元的额外要求**:
- 需实现**等参元**框架：形函数、Jacobi 矩阵、B 矩阵
- 需要 **Gauss 积分**来计算 $\mathbf{K}^e = \int \mathbf{B}^T \mathbf{D} \mathbf{B} \, dA$
- 材料矩阵 $\mathbf{D}$ 需扩展（弯曲刚度矩阵 $\mathbf{D}_b = \frac{Eh^3}{12(1-\nu^2)}$）
- 建议在 `CElement` 基类中增加数值积分的通用接口

### 8.3 添加壳（Shell）单元

壳单元是最复杂的扩展，通常有两种实现策略：

**方案 A** — 退化壳单元（Degenerated Shell）:
- 从三维实体单元退化得到
- 需要节点法向量信息
- 每节点 5~6 个自由度

**方案 B** — 平板壳单元（Flat Shell，推荐入门）:
- 将**膜单元**（面内）与**板弯曲单元**（面外）叠加
- 每节点 6 个自由度 = 膜 3 个 + 弯曲 3 个
- 刚度矩阵 = 膜刚度矩阵 + 板弯曲刚度矩阵

```
K_shell = [ K_membrane    0        ]
          [    0       K_bending    ]
```

**壳单元需要的额外基础设施**：
- 局部坐标系的建立（基于单元法向量）
- 从局部到全局的坐标变换
- 可能需要 drilling 自由度（θz）的稳定化处理

### 8.4 需要注意的框架层面修改

扩展梁板壳单元时，以下框架层面的问题需特别注意：

1. **`NDF` 的处理**:
   - 当前 `NDF=3` 是硬编码的 `static const`
   - 若混用杆和梁，节点自由度数不一致
   - 建议改为统一取最大值 `NDF=6`，对不需要转动自由度的节点约束 θx, θy, θz

2. **输入文件格式**:
   - 需要与 STAP90 格式兼容
   - 梁需要额外输入截面参数；板壳需要厚度
   - 节点数据中 `bcode` 数组长度需与 `NDF` 匹配

3. **应力输出**:
   - 杆：只有一个标量应力
   - 梁：弯矩、剪力、轴力（至少 6 个内力分量）
   - 板壳：应力张量（$\sigma_{xx}, \sigma_{yy}, \sigma_{xy}$）或广义内力

4. **`ElementStress` 接口**:
   - 当前接口 `ElementStress(double* stress, double* Displacement)` 只输出一个标量
   - 对梁板壳，应力是多分量的，建议修改接口或让 `stress` 指向数组

5. **CMakeLists.txt**:
   - 无需手动修改，新文件放在 `cpp/` 和 `h/` 目录下会被 CMake 自动发现

### 8.5 推荐的扩展顺序

```
Bar (已完成) → Beam (2D先) → Beam (3D) → T3膜 → Q4膜 → 板弯曲 → 平板壳 → 退化壳
```

每一步都建议用经典算例验证（悬臂梁、简支板等），确保刚度矩阵和应力计算正确后再进入下一步

---

## 9. 关键代码片段注释

### 9.1 main.cpp —— 主控流程

```cpp
int main(int argc, char *argv[])
{
    // ── 命令行参数处理 ──
    // 程序接受一个参数：输入文件名（可以带.dat后缀，也可以不带）
    string filename(argv[1]);
    string InFile = filename + ".dat";   // 输入文件
    string OutFile = filename + ".out";  // 输出文件

    // ── 获取全局唯一的问题域实例 ──
    CDomain* FEMData = CDomain::GetInstance();

    // ── 阶段1: 数据读入 ──
    // 读入全部数据：标题、控制参数、节点、荷载、单元
    // 同时进行方程编号和数据校验
    FEMData->ReadData(InFile, OutFile);

    // MODEX=0 时只检查数据不求解
    if (!FEMData->GetMODEX()) return 0;

    // ── 阶段2: 分配存储并计算列高 ──
    // 这一步确定了 Skyline 矩阵的存储结构
    FEMData->AllocateMatrices();

    // ── 阶段3: 组装总刚度矩阵 ──
    // 遍历所有单元，计算单元刚度矩阵并组装到全局矩阵
    FEMData->AssembleStiffnessMatrix();

    // ── 阶段4: LDLT 分解 ──
    // 对总刚度矩阵 K = L·D·L^T 进行原位分解
    // 分解完成后 K 的存储空间被 L 和 D 覆写
    CLDLTSolver* Solver = new CLDLTSolver(FEMData->GetStiffnessMatrix());
    Solver->LDLT();

    // ── 阶段5: 逐工况求解 ──
    // 【关键】：LDLT 只做一次，多个荷载工况只需重复组装力向量 + 回代
    for (unsigned int lcase = 0; lcase < FEMData->GetNLCASE(); lcase++)
    {
        FEMData->AssembleForce(lcase + 1);         // 组装第 lcase+1 个荷载工况的力向量
        Solver->BackSubstitution(FEMData->GetForce()); // 回代求解 → 力向量被位移覆写
        Output->OutputNodalDisplacement();          // 输出节点位移
        Output->OutputElementStress();              // 计算并输出单元应力
    }
}
```

### 9.2 CDomain::CalculateEquationNumber —— 方程编号

```cpp
void CDomain::CalculateEquationNumber()
{
    NEQ = 0;  // 全局方程计数器归零
    for (unsigned int np = 0; np < NUMNP; np++)         // 遍历所有节点
    {
        for (unsigned int dof = 0; dof < CNode::NDF; dof++)  // 遍历节点的 x,y,z 三个自由度
        {
            if (NodeList[np].bcode[dof])
                // 输入中 bcode=1 表示约束，置为 0（无对应方程）
                NodeList[np].bcode[dof] = 0;
            else
            {
                // 输入中 bcode=0 表示自由，分配递增的全局方程号
                NEQ++;
                NodeList[np].bcode[dof] = NEQ;
                // 此后 bcode 不再是边界码，而是全局方程编号
            }
        }
    }
    // 最终 NEQ 就是总刚度矩阵的维度
}
```

### 9.3 CBar::ElementStiffness —— 杆单元刚度矩阵计算

```cpp
void CBar::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());  // 清零 21 个元素

    // ── 计算杆的几何信息 ──
    double DX[3];  // 两节点坐标差 (dx, dy, dz)
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    // 计算方向余弦的二次多项式
    // DX2 = {dx², dy², dz², dx·dy, dy·dz, dx·dz}
    double DX2[6];
    DX2[0] = DX[0] * DX[0];  // dx²
    DX2[1] = DX[1] * DX[1];  // dy²
    DX2[2] = DX[2] * DX[2];  // dz²
    DX2[3] = DX[0] * DX[1];  // dx·dy
    DX2[4] = DX[1] * DX[2];  // dy·dz
    DX2[5] = DX[0] * DX[2];  // dx·dz

    double L2 = DX2[0] + DX2[1] + DX2[2]; // L² = dx² + dy² + dz²
    double L = sqrt(L2);                    // L = 杆长

    // ── 计算刚度系数 ──
    // k = EA / L³ (注意: EA/L³ 而非 EA/L，因为 DX2 中已含 L² 量级)
    CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);
    double k = material_->E * material_->Area / L / L2;

    // ── 填充上三角刚度矩阵 (按列存储) ──
    // 6×6 矩阵的上三角共 21 个元素
    // 按列顺序: 第1列对角元 → 第2列(对角元,上方元素) → ... → 第6列
    //
    // 物理含义: K^e = (EA/L³) * [  D  -D ]
    //                            [ -D   D ]
    // 其中 D = [ dx²      dx·dy   dx·dz ]
    //          [ dx·dy    dy²     dy·dz ]
    //          [ dx·dz    dy·dz   dz²   ]
    
    Matrix[0]  =  k*DX2[0];    // K(1,1) = k·dx²
    Matrix[1]  =  k*DX2[1];    // K(2,2) = k·dy²
    Matrix[2]  =  k*DX2[3];    // K(1,2) = k·dx·dy
    Matrix[3]  =  k*DX2[2];    // K(3,3) = k·dz²
    Matrix[4]  =  k*DX2[4];    // K(2,3) = k·dy·dz
    Matrix[5]  =  k*DX2[5];    // K(1,3) = k·dx·dz
    Matrix[6]  =  k*DX2[0];    // K(4,4) = k·dx²
    Matrix[7]  = -k*DX2[5];    // K(3,4) = -k·dx·dz
    Matrix[8]  = -k*DX2[3];    // K(2,4) = -k·dx·dy
    Matrix[9]  = -k*DX2[0];    // K(1,4) = -k·dx²
    Matrix[10] =  k*DX2[1];    // K(5,5) = k·dy²
    Matrix[11] =  k*DX2[3];    // K(4,5) = k·dx·dy  (注意无负号)
    // ... 以此类推，Matrix[12]~Matrix[20] 为第5、6列的元素
}
```

### 9.4 CSkylineMatrix::Assembly —— 总刚度矩阵组装

```cpp
void CSkylineMatrix<T_>::Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND)
{
    // 遍历单元刚度矩阵的每一列 j
    for (unsigned int j = 0; j < ND; j++)
    {
        unsigned int Lj = LocationMatrix[j];  // 单元第j个DOF → 全局方程号
        if (!Lj) continue;                     // Lj=0 表示该自由度被约束，跳过

        // 单元刚度矩阵按列存储，第j列对角元在一维数组中的偏移
        unsigned int DiagjElement = (j+1)*j/2;

        // 遍历第j列的第0~j行（上三角部分）
        for (unsigned int i = 0; i <= j; i++)
        {
            unsigned int Li = LocationMatrix[i];  // 单元第i个DOF → 全局方程号
            if (!Li) continue;

            // 单元矩阵 Matrix[DiagjElement + j - i] 对应 K^e(i,j)
            // 累加到全局矩阵 K(Li, Lj)
            (*this)(Li, Lj) += Matrix[DiagjElement + j - i];
        }
    }
}
```

### 9.5 CLDLTSolver::LDLT —— LDLT 分解

```cpp
void CLDLTSolver::LDLT()
{
    unsigned int N = K.dim();                        // 方程总数
    unsigned int* ColumnHeights = K.GetColumnHeights();

    for (unsigned int j = 2; j <= N; j++)  // 从第2列开始逐列处理
    {
        // mj: 第j列第一个非零元的行号 (skyline 边界)
        unsigned int mj = j - ColumnHeights[j-1];

        // ── 步骤1: 修正第j列中 i=mj+1..j-1 行的元素 ──
        // K(i,j) → U(i,j) = K(i,j) - Σ L(r,i)·U(r,j)
        for (unsigned int i = mj+1; i <= j-1; i++)
        {
            unsigned int mi = i - ColumnHeights[i-1];
            double C = 0.0;
            for (unsigned int r = max(mi, mj); r <= i-1; r++)
                C += K(r,i) * K(r,j);  // L(r,i) × U(r,j)
            K(i,j) -= C;
        }

        // ── 步骤2: 计算 L(r,j) 并更新对角元 D(j,j) ──
        for (unsigned int r = mj; r <= j-1; r++)
        {
            double Lrj = K(r,j) / K(r,r);  // L(r,j) = U(r,j) / D(r,r)
            K(j,j) -= Lrj * K(r,j);         // D(j,j) -= L(r,j) × U(r,j)
            K(r,j) = Lrj;                   // 用 L(r,j) 覆写 U(r,j) 的存储位置
        }

        // 对角元为零或负说明矩阵不正定
        if (fabs(K(j,j)) <= FLT_MIN) {
            cerr << "Stiffness matrix is not positive definite!" << endl;
            exit(4);
        }
    }
}
```

### 9.6 CLDLTSolver::BackSubstitution —— 回代求解

```cpp
void CLDLTSolver::BackSubstitution(double* Force)
{
    unsigned int N = K.dim();
    unsigned int* ColumnHeights = K.GetColumnHeights();

    // ── 第一步: 前代求解 LV = R ──
    // V(i) = R(i) - Σ_{j=mi}^{i-1} L(j,i) × V(j)
    for (unsigned int i = 2; i <= N; i++)
    {
        unsigned int mi = i - ColumnHeights[i-1];
        for (unsigned int j = mi; j <= i-1; j++)
            Force[i-1] -= K(j,i) * Force[j-1];
    }

    // ── 第二步: 对角缩放 V̄ = D⁻¹V ──
    for (unsigned int i = 1; i <= N; i++)
        Force[i-1] /= K(i,i);

    // ── 第三步: 回代求解 L^T a = V̄ ──
    // a(i) = V̄(i) - Σ_{j>i} L(i,j) × a(j)
    // 从最后一个方程往前推
    for (unsigned int j = N; j >= 2; j--)
    {
        unsigned int mj = j - ColumnHeights[j-1];
        for (unsigned int i = mj; i <= j-1; i++)
            Force[i-1] -= K(i,j) * Force[j-1];
    }
    // 最终 Force[] 中存储的就是位移向量 a
}
```

### 9.7 CSkylineMatrix::CalculateColumnHeight —— 列高计算

```cpp
void CSkylineMatrix<T_>::CalculateColumnHeight(unsigned int* LocationMatrix, size_t ND)
{
    // 找该单元所有活跃自由度中最小的全局方程号
    // 这决定了 skyline 的"天际线"位置
    unsigned int nfirstrow = INT_MAX;
    for (unsigned int i = 0; i < ND; i++)
        if (LocationMatrix[i] && LocationMatrix[i] < nfirstrow)
            nfirstrow = LocationMatrix[i];

    // 对每个活跃自由度，计算它所在列需要延伸到的高度
    for (unsigned int i = 0; i < ND; i++)
    {
        unsigned int column = LocationMatrix[i];
        if (!column) continue;  // 约束自由度跳过

        unsigned int Height = column - nfirstrow;
        // 取所有单元贡献中的最大值
        if (ColumnHeights_[column-1] < Height)
            ColumnHeights_[column-1] = Height;
    }
    // 列高决定了每列需要存储多少个非零元
    // 列高越大 → 存储越多 → 半带宽越大
}
```

### 9.8 CBar::ElementStress —— 杆单元应力计算

```cpp
void CBar::ElementStress(double* stress, double* Displacement)
{
    CBarMaterial* material_ = dynamic_cast<CBarMaterial*>(ElementMaterial_);

    double DX[3];
    double L2 = 0;
    for (unsigned int i = 0; i < 3; i++)
    {
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];
        L2 += DX[i] * DX[i];  // L² = Σ(dx²)
    }

    // 应力系数向量 S = E/L² × [-dx, -dy, -dz, dx, dy, dz]
    // 物理含义: σ = (E/L) × ε = (E/L) × (ΔL/L) = (E/L²) × (d2-d1)·n̂×L
    double S[6];
    for (unsigned int i = 0; i < 3; i++)
    {
        S[i]   = -DX[i] * material_->E / L2;  // 节点1对应的系数 (负号)
        S[i+3] = -S[i];                         // 节点2对应的系数 (正号)
    }

    // σ = S · d^e (应力 = 系数向量 × 单元节点位移)
    *stress = 0.0;
    for (unsigned int i = 0; i < 6; i++)
    {
        if (LocationMatrix_[i])  // 只累加活跃自由度
            *stress += S[i] * Displacement[LocationMatrix_[i]-1];
            // LocationMatrix_[i]-1 : 方程号从1开始，数组索引从0开始
    }
}
```

### 9.9 CDomain::AssembleForce —— 荷载向量组装

```cpp
bool CDomain::AssembleForce(unsigned int LoadCase)
{
    CLoadCaseData* LoadData = &LoadCases[LoadCase - 1];
    
    clear(Force, NEQ);  // 清零全局力向量（每个工况都需要重新组装）

    // 遍历当前工况的所有集中荷载
    for (unsigned int lnum = 0; lnum < LoadData->nloads; lnum++)
    {
        // 查找荷载对应的全局方程号
        // node[lnum]-1   : 节点数组索引（从0开始）
        // dof[lnum]-1    : 自由度数组索引（从0开始）
        // bcode[...]     : 经过编号后存储的全局方程号
        unsigned int dof = NodeList[LoadData->node[lnum] - 1]
                                   .bcode[LoadData->dof[lnum] - 1];

        if (dof)  // dof>0 表示该自由度是活跃的（未被约束）
            Force[dof - 1] += LoadData->load[lnum];
            // 将荷载值累加到对应的方程位置
    }
    return true;
}
```

### 9.10 CElementGroup::Read —— 单元组数据读入

```cpp
bool CElementGroup::Read(ifstream& Input)
{
    // 读入: 单元类型、单元数量、材料组数量
    Input >> (int&)ElementType_ >> NUME_ >> NUMMAT_;

    // 根据单元类型确定派生类大小（用于多态数组的内存管理）
    CalculateMemberSize();  // 设置 ElementSize_ 和 MaterialSize_

    // ── 读入材料数据 ──
    AllocateMaterials(NUMMAT_);  // 按派生类型分配 (如 new CBarMaterial[NUMMAT_])
    for (unsigned int mset = 0; mset < NUMMAT_; mset++)
        GetMaterial(mset).Read(Input);  // 通过字节偏移访问第 mset 个材料对象

    // ── 读入单元数据 ──
    AllocateElements(NUME_);  // 按派生类型分配 (如 new CBar[NUME_])
    for (unsigned int Ele = 0; Ele < NUME_; Ele++)
    {
        unsigned int N;
        Input >> N;  // 读入单元编号
        // operator[] 通过字节偏移定位到正确的派生类对象
        (*this)[Ele].Read(Input, MaterialList_, NodeList_);
    }
    return true;
}
```

---

## 10. 总结

| 特性             | 描述                           |
| ---------------- | ------------------------------ |
| **求解问题**     | 线性静力学                     |
| **已实现单元**   | 三维杆单元（Bar/Truss）        |
| **刚度矩阵存储** | Skyline 格式（变带宽稀疏存储） |
| **求解方法**     | LDLT 分解 + 回代               |
| **编程范式**     | 面向对象 + 模板                |
| **核心设计模式** | 单例、模板方法、工厂           |
| **跨平台**       | CMake 构建系统                 |

STAPpp 虽然代码量不大（约 1500 行），但完整展示了有限元程序的核心框架：**前处理（数据读入与编号）→ 矩阵组装 → 方程求解 → 后处理（位移与应力输出）**，是学习有限元编程实现的优秀教学代码。
