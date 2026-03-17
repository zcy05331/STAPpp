# STAPpp 项目代码实现详细分析

## 项目概述

STAPpp (Structural Analysis Program Plus Plus) 是一个用C++实现的有限元分析程序，由清华大学航天航空学院计算动力学实验室开发。该程序与STAP90共享相同的输入数据文件格式。

**版本信息**: Release 1.11, November 22, 2017

## 核心架构设计

### 1. 设计模式

#### 1.1 单例模式 (Singleton Pattern)
项目中多个核心类采用单例模式，确保全局只有一个实例：

- **CDomain**: 问题域类，管理整个有限元模型
- **COutputter**: 输出管理类，处理所有输出操作

```cpp
// CDomain单例实现
static CDomain* _instance;
static CDomain* GetInstance() {
    if (!_instance)
        _instance = new CDomain();
    return _instance;
}
```

#### 1.2 工厂模式与多态
使用抽象基类和虚函数实现多态：

- **CElement**: 元素基类
- **CMaterial**: 材料基类
- 派生类如 CBar, CBarMaterial 实现具体功能

### 2. 核心类结构

#### 2.1 CDomain (问题域类)
**职责**: 管理整个有限元分析的数据和流程

**关键成员变量**:
```cpp
char Title[256];                    // 问题标题
unsigned int MODEX;                 // 求解模式 (0:数据检查, 1:执行)
unsigned int NUMNP;                 // 节点总数
CNode* NodeList;                    // 节点列表
unsigned int NUMEG;                 // 单元组数量
CElementGroup* EleGrpList;          // 单元组列表
unsigned int NLCASE;                // 载荷工况数
CLoadCaseData* LoadCases;           // 载荷工况列表
unsigned int NEQ;                   // 方程总数
CSkylineMatrix<double>* StiffnessMatrix;  // 刚度矩阵(天际线存储)
double* Force;                      // 全局力/位移向量
```

**关键方法**:
- `ReadData()`: 读取输入文件
- `CalculateEquationNumber()`: 计算全局方程编号
- `AllocateMatrices()`: 分配矩阵存储空间
- `AssembleStiffnessMatrix()`: 组装全局刚度矩阵
- `AssembleForce()`: 组装全局载荷向量

#### 2.2 CElement (单元基类)
**职责**: 定义所有单元类型的通用接口

**关键成员**:
```cpp
unsigned int NEN_;              // 单元节点数
CNode** nodes_;                 // 单元节点指针数组
CMaterial* ElementMaterial_;    // 单元材料
unsigned int* LocationMatrix_;  // 位置矩阵(全局自由度编号)
unsigned int ND_;              // 位置矩阵维度
```

**纯虚函数**:
```cpp
virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList) = 0;
virtual void Write(COutputter& output) = 0;
virtual void ElementStiffness(double* stiffness) = 0;
virtual void ElementStress(double* stress, double* Displacement) = 0;
```

#### 2.3 CBar (杆单元类)
**职责**: 实现三维杆单元

**单元特性**:
- 节点数: 2 (NEN_ = 2)
- 自由度: 6 (每个节点3个平动自由度)
- 刚度矩阵大小: 21 (上三角存储)

**单元刚度矩阵计算** (ElementStiffness):
```cpp
// 计算杆长度和方向余弦
DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];  // i=0,1,2
L2 = DX[0]^2 + DX[1]^2 + DX[2]^2;
L = sqrt(L2);

// 刚度系数
k = E * Area / (L * L2);

// 刚度矩阵元素 (上三角存储)
K_ij = k * dx_i * dx_j  // 其中dx为方向余弦分量
```

**单元应力计算** (ElementStress):
```cpp
// 应力 = E/L^2 * [-dx, -dy, -dz, dx, dy, dz] * {u}
σ = E/L^2 * Σ(S[i] * Displacement[LocationMatrix[i]-1])
```

### 3. 数据结构

#### 3.1 CSkylineMatrix (天际线矩阵)
**职责**: 高效存储对称稀疏刚度矩阵

**存储策略**:
- 只存储天际线以下的元素
- 利用对称性，只存储上三角或下三角
- 大幅减少内存占用

**关键成员**:
```cpp
T_* data_;                      // 矩阵数据(一维数组)
unsigned int NEQ_;              // 矩阵维度
unsigned int MK_;               // 最大半带宽
unsigned int NWK_;              // 存储空间大小
unsigned int* ColumnHeights_;   // 列高数组
unsigned int* DiagonalAddress_; // 对角元素地址
```

**索引计算** (operator()):
```cpp
// 访问元素(i,j), 编号从1开始
if (j >= i)
    return data_[DiagonalAddress_[j-1] + (j-i) - 1];
else
    return data_[DiagonalAddress_[i-1] + (i-j) - 1];
```

**列高计算** (CalculateColumnHeight):
```cpp
// 找到单元位置矩阵中第一个非零行号
nfirstrow = min(LocationMatrix[i]) where LocationMatrix[i] != 0;

// 更新列高
for each column in LocationMatrix:
    Height = column - nfirstrow;
    if (ColumnHeights[column-1] < Height)
        ColumnHeights[column-1] = Height;
```

**对角地址计算** (CalculateDiagnoalAddress):
```cpp
// 递推公式: M(i+1) = M(i) + H(i) + 1
DiagonalAddress[0] = 1;
for (col = 1; col <= NEQ; col++)
    DiagonalAddress[col] = DiagonalAddress[col-1] + ColumnHeights[col-1] + 1;
```

#### 3.2 CNode (节点类)
**关键成员**:
```cpp
static const unsigned int NDF = 3;  // 每节点自由度数
unsigned int NodeNumber;            // 节点编号
double XYZ[3];                      // 节点坐标
unsigned int bcode[NDF];            // 边界条件码/全局方程号
```

**bcode的双重用途**:
1. 读入时: 0=自由, 1=约束
2. 计算后: 存储全局方程编号 (从1开始)

#### 3.3 CElementGroup (单元组类)
**职责**: 管理相同类型的单元集合

**设计特点**:
- 使用指针算术访问派生类数组
- 动态计算派生类大小 (sizeof)
- 支持多种单元类型扩展

```cpp
// 重载[]运算符，通过指针偏移访问元素
CElement& operator[](unsigned int i) {
    return *(CElement*)((std::size_t)(ElementList_) + i*ElementSize_);
}
```

## 求解流程

### 主程序流程 (main.cpp)

```
1. 读取输入文件 (.dat)
   └─> CDomain::ReadData()
       ├─> ReadNodalPoints()      // 读取节点
       ├─> CalculateEquationNumber()  // 计算方程编号
       ├─> ReadLoadCases()        // 读取载荷
       └─> ReadElements()         // 读取单元

2. 分配矩阵存储
   └─> CDomain::AllocateMatrices()
       ├─> CalculateColumnHeights()    // 计算列高
       ├─> CalculateDiagnoalAddress()  // 计算对角地址
       └─> Allocate()                  // 分配存储空间

3. 组装刚度矩阵
   └─> CDomain::AssembleStiffnessMatrix()
       └─> for each element:
           ├─> ElementStiffness()      // 计算单元刚度
           └─> Assembly()              // 组装到全局矩阵

4. 求解线性方程组
   └─> CLDLTSolver
       ├─> LDLT()                      // LDL^T分解
       └─> for each load case:
           ├─> AssembleForce()         // 组装载荷向量
           ├─> BackSubstitution()      // 回代求解
           ├─> OutputNodalDisplacement()  // 输出位移
           └─> OutputElementStress()   // 输出应力
```

### LDLT求解器 (CLDLTSolver)

#### LDLT分解算法
将对称正定矩阵K分解为: K = L * D * L^T
- L: 单位下三角矩阵
- D: 对角矩阵

**实现细节**:
```cpp
for j = 2 to N:  // 列循环
    mj = j - ColumnHeights[j-1];  // 第j列第一个非零元素行号

    // 计算U_ij (上三角元素)
    for i = mj+1 to j-1:
        mi = i - ColumnHeights[i-1];
        C = Σ(K(r,i) * K(r,j)) for r = max(mi,mj) to i-1;
        K(i,j) -= C;

    // 计算L_rj和D_jj
    for r = mj to j-1:
        L_rj = K(r,j) / K(r,r);
        K(j,j) -= L_rj * K(r,j);
        K(r,j) = L_rj;  // 原位存储
```

#### 回代求解
求解 K*u = F，分三步:
1. 前代: L*v = F
2. 对角: v_bar = D^(-1)*v
3. 回代: L^T*u = v_bar

```cpp
// 前代
for i = 2 to N:
    mi = i - ColumnHeights[i-1];
    Force[i-1] -= Σ(K(j,i) * Force[j-1]) for j = mi to i-1;

// 对角缩放
for i = 1 to N:
    Force[i-1] /= K(i,i);

// 回代
for j = N down to 2:
    mj = j - ColumnHeights[j-1];
    for i = mj to j-1:
        Force[i-1] -= K(i,j) * Force[j-1];
```

## 关键算法实现

### 1. 全局方程编号
```cpp
void CDomain::CalculateEquationNumber() {
    NEQ = 0;
    for each node:
        for each DOF:
            if (bcode[dof] == 0):  // 自由度
                NEQ++;
                bcode[dof] = NEQ;  // 赋予全局方程号
            else:
                bcode[dof] = 0;    // 约束自由度
}
```

### 2. 位置矩阵生成
```cpp
void CElement::GenerateLocationMatrix() {
    unsigned int i = 0;
    for (N = 0; N < NEN_; N++)          // 遍历单元节点
        for (D = 0; D < NDF; D++)       // 遍历节点自由度
            LocationMatrix_[i++] = nodes_[N]->bcode[D];
}
```

### 3. 刚度矩阵组装
```cpp
void CSkylineMatrix::Assembly(double* Matrix, unsigned int* LocationMatrix, size_t ND) {
    for j = 0 to ND-1:
        Lj = LocationMatrix[j];  // 全局方程号
        if (!Lj) continue;       // 跳过约束自由度

        DiagjElement = (j+1)*j/2;  // 单元矩阵中第j列对角元素位置

        for i = 0 to j:
            Li = LocationMatrix[i];
            if (!Li) continue;

            // 组装: K(Li,Lj) += Ke[DiagjElement + j - i]
            (*this)(Li,Lj) += Matrix[DiagjElement + j - i];
}
```

## 内存管理

### 动态内存分配策略
1. **节点数组**: `new CNode[NUMNP]`
2. **单元组数组**: `new CElementGroup[NUMEG]`
3. **单元数组**: 根据类型动态分配 (如 `new CBar[NUME]`)
4. **刚度矩阵**: `new double[NWK]` (天际线存储)
5. **载荷向量**: `new double[NEQ]`

### 析构函数链
```
~CDomain()
  ├─> delete[] NodeList
  ├─> delete[] EleGrpList
  │     └─> ~CElementGroup()
  │           ├─> delete[] ElementList_
  │           │     └─> ~CBar() / ~CElement()
  │           └─> delete[] MaterialList_
  ├─> delete[] LoadCases
  ├─> delete[] Force
  └─> delete StiffnessMatrix
        └─> ~CSkylineMatrix()
              ├─> delete[] ColumnHeights_
              ├─> delete[] DiagonalAddress_
              └─> delete[] data_
```

## 输入输出格式

### 输入文件结构 (.dat)
```
标题行
NUMNP NUMEG NLCASE MODEX
节点数据块
载荷工况数据块
单元组数据块
  ElementType NUME NUMMAT
  材料数据
  单元数据
```

### 输出文件内容 (.out)
1. 问题标题和时间戳
2. 节点信息
3. 方程编号
4. 载荷信息
5. 单元信息
6. 系统数据统计
7. 每个载荷工况的:
   - 节点位移
   - 单元应力
8. 求解时间统计

## 调试功能

### 条件编译宏 _DEBUG_
启用时输出额外调试信息:
- 完整刚度矩阵
- 对角地址数组
- 列高数组
- 位移向量
- 位置矩阵

```cpp
#ifdef _DEBUG_
    Output->PrintStiffnessMatrix();
    Output->PrintDiagonalAddress();
    Output->PrintColumnHeights();
    Output->PrintDisplacement();
#endif
```

## 性能优化

### 1. 天际线存储
- 相比全矩阵存储，节省大量内存
- 对于稀疏矩阵，存储效率显著提高
- 存储量: NWK = Σ(ColumnHeights[i] + 1)

### 2. 原位分解
- LDLT分解直接在原矩阵上进行
- 不需要额外的L和D矩阵存储
- L存储在下三角，D存储在对角线

### 3. 单元组管理
- 相同类型单元集中处理
- 减少虚函数调用开销
- 便于批量操作

### 4. 指针算术
- 使用指针偏移访问数组元素
- 避免虚函数表查找
- 提高访问效率

## 可扩展性设计

### 添加新单元类型的步骤

1. **定义单元类** (继承CElement)
```cpp
class CNewElement : public CElement {
    // 实现纯虚函数
    virtual bool Read(...);
    virtual void Write(...);
    virtual void ElementStiffness(...);
    virtual void ElementStress(...);
};
```

2. **定义材料类** (继承CMaterial)
```cpp
class CNewMaterial : public CMaterial {
    // 添加材料属性
    // 实现Read和Write
};
```

3. **更新ElementTypes枚举**
```cpp
enum ElementTypes {
    ...,
    NewElementType
};
```

4. **更新CElementGroup**
```cpp
void CElementGroup::CalculateMemberSize() {
    case ElementTypes::NewElementType:
        ElementSize_ = sizeof(CNewElement);
        MaterialSize_ = sizeof(CNewMaterial);
        break;
}

void CElementGroup::AllocateElements(size_t size) {
    case ElementTypes::NewElementType:
        ElementList_ = new CNewElement[size];
        break;
}
```

## 代码质量特点

### 优点
1. **清晰的架构**: 单一职责原则，类职责明确
2. **良好的封装**: 使用private/protected/public合理分离接口
3. **可扩展性**: 基于继承和多态的设计便于添加新单元
4. **高效算法**: 天际线存储和LDLT求解器性能优秀
5. **完整注释**: Doxygen风格注释，便于生成文档

### 可改进之处
1. **智能指针**: 使用std::unique_ptr/shared_ptr替代裸指针
2. **异常处理**: 使用C++异常机制替代exit()
3. **STL容器**: 使用std::vector替代手动管理的数组
4. **命名空间**: 避免using namespace std污染全局命名空间
5. **const正确性**: 更多使用const修饰符
6. **移动语义**: 利用C++11移动语义优化性能

## 总结

STAPpp是一个设计良好的有限元分析程序，展示了以下核心技术:

1. **面向对象设计**: 合理使用继承、多态、封装
2. **高效数据结构**: 天际线矩阵存储
3. **数值算法**: LDLT分解和回代求解
4. **软件工程**: 单例模式、工厂模式、模板编程
5. **可扩展架构**: 便于添加新单元类型

该项目为学习有限元方法的程序实现提供了优秀的参考，代码结构清晰，算法实现高效，是理解有限元分析软件开发的良好范例。
