# STAPpp 添加新单元类型完整指南

## 概述

本指南详细说明如何在STAPpp中添加新的单元类型（如梁单元、板单元、壳单元等）。以添加**梁单元(Beam)**为例，展示完整的实现步骤。

---

## 必须修改的文件清单

### 核心文件（必须修改）
1. **src/h/ElementGroup.h** - 添加单元类型枚举（已有定义）
2. **src/h/Beam.h** - 创建梁单元头文件（新建）
3. **src/cpp/Beam.cpp** - 创建梁单元实现文件（新建）
4. **src/h/Material.h** - 添加梁材料类定义
5. **src/cpp/Material.cpp** - 实现梁材料类方法
6. **src/cpp/ElementGroup.cpp** - 更新单元组管理逻辑
7. **src/cpp/Outputter.cpp** - 添加梁单元输出函数
8. **src/h/Outputter.h** - 声明梁单元输出函数
9. **src/CMakeLists.txt** - 添加新源文件到编译系统

### 可选文件
10. **src/h/Node.h** - 如果需要更多自由度（如转角），需修改NDF

---

## 详细实现步骤

### 步骤1: 确认单元类型枚举（已存在）

**文件**: `src/h/ElementGroup.h`

```cpp
enum ElementTypes
{
    UNDEFINED = 0,
    Bar,    // Bar element (已实现)
    Q4,     // 4Q element
    T3,     // 3T element
    H8,     // 8H element
    Beam,   // Beam element ← 我们要实现这个
    Plate,  // Plate element
    Shell   // Shell element
};
```

✅ **已定义，无需修改**

---

### 步骤2: 创建梁单元头文件

**文件**: `src/h/Beam.h` （新建）

```cpp
/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Beam element class (3D Euler-Bernoulli beam)
class CBeam : public CElement
{
public:

    //! Constructor
    CBeam();

    //! Destructor
    ~CBeam();

    //! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //! Write element data to stream
    virtual void Write(COutputter& output);

    //! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //! Calculate element stress (for beam: bending moment, shear force, axial force)
    virtual void ElementStress(double* stress, double* Displacement);

private:
    //! Calculate transformation matrix from local to global coordinates
    void GetTransformationMatrix(double* T);

    //! Calculate local stiffness matrix
    void LocalStiffness(double* Klocal);
};
```

**关键点**:
- 继承自 `CElement`
- 必须实现4个纯虚函数
- 梁单元通常需要坐标转换矩阵

---

### 步骤3: 实现梁单元类

**文件**: `src/cpp/Beam.cpp` （新建）

```cpp
/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*****************************************************************************/

#include "Beam.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//! Constructor
CBeam::CBeam()
{
    NEN_ = 2;    // 梁单元有2个节点
    nodes_ = new CNode*[NEN_];

    // 每个节点6个自由度 (3平动 + 3转动)
    // 如果Node::NDF=3，需要修改为6
    ND_ = 12;    // 2节点 × 6自由度/节点 = 12
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

//! Destructor
CBeam::~CBeam()
{
}

//! Read element data from stream Input
bool CBeam::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;    // Material property set number
    unsigned int N1, N2;  // Node numbers

    Input >> N1 >> N2 >> MSet;

    ElementMaterial_ = dynamic_cast<CBeamMaterial*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];

    return true;
}

//! Write element data to stream
void CBeam::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

//! Calculate element stiffness matrix
void CBeam::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);

    // 计算梁长度
    double DX[3];
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2]);
    double L2 = L * L;
    double L3 = L2 * L;

    // 局部坐标系下的刚度矩阵
    double* Klocal = new double[ND_ * ND_];
    clear(Klocal, ND_ * ND_);

    double E = material->E;
    double A = material->Area;
    double Iy = material->Iy;  // 绕y轴的惯性矩
    double Iz = material->Iz;  // 绕z轴的惯性矩
    double J = material->J;    // 扭转常数
    double G = material->G;    // 剪切模量

    // 轴向刚度
    double EA_L = E * A / L;

    // 弯曲刚度 (绕y轴)
    double EIy_L3 = E * Iy / L3;

    // 弯曲刚度 (绕z轴)
    double EIz_L3 = E * Iz / L3;

    // 扭转刚度
    double GJ_L = G * J / L;

    // 填充局部刚度矩阵 (12×12)
    // 这里简化处理，实际需要完整的12×12矩阵
    // 参考有限元教材中的Euler-Bernoulli梁单元刚度矩阵

    // 轴向自由度 (u1, u7)
    Klocal[0*ND_ + 0] = EA_L;
    Klocal[0*ND_ + 6] = -EA_L;
    Klocal[6*ND_ + 0] = -EA_L;
    Klocal[6*ND_ + 6] = EA_L;

    // 绕z轴弯曲 (v1, θz1, v7, θz7)
    Klocal[1*ND_ + 1] = 12*EIz_L3;
    Klocal[1*ND_ + 5] = 6*EIz_L3*L;
    Klocal[1*ND_ + 7] = -12*EIz_L3;
    Klocal[1*ND_ + 11] = 6*EIz_L3*L;

    Klocal[5*ND_ + 1] = 6*EIz_L3*L;
    Klocal[5*ND_ + 5] = 4*EIz_L3*L2;
    Klocal[5*ND_ + 7] = -6*EIz_L3*L;
    Klocal[5*ND_ + 11] = 2*EIz_L3*L2;

    // ... 其他项类似填充

    // 坐标转换
    double* T = new double[ND_ * ND_];
    GetTransformationMatrix(T);

    // K_global = T^T * K_local * T
    // 这里需要矩阵乘法运算
    // 简化处理，实际需要完整实现

    delete[] Klocal;
    delete[] T;
}

//! Calculate transformation matrix
void CBeam::GetTransformationMatrix(double* T)
{
    // 计算从局部坐标系到全局坐标系的转换矩阵
    // T是12×12的矩阵

    double DX[3];
    for (unsigned int i = 0; i < 3; i++)
        DX[i] = nodes_[1]->XYZ[i] - nodes_[0]->XYZ[i];

    double L = sqrt(DX[0]*DX[0] + DX[1]*DX[1] + DX[2]*DX[2]);

    // 方向余弦
    double cx = DX[0] / L;
    double cy = DX[1] / L;
    double cz = DX[2] / L;

    // 构造转换矩阵（3×3旋转矩阵的块对角形式）
    // 详细实现参考有限元教材
}

//! Calculate element stress
void CBeam::ElementStress(double* stress, double* Displacement)
{
    // 对于梁单元，应力包括：
    // - 轴力 N
    // - 剪力 Vy, Vz
    // - 弯矩 My, Mz
    // - 扭矩 T

    CBeamMaterial* material = dynamic_cast<CBeamMaterial*>(ElementMaterial_);

    // 提取单元位移
    double u[12];
    for (unsigned int i = 0; i < ND_; i++)
    {
        if (LocationMatrix_[i])
            u[i] = Displacement[LocationMatrix_[i]-1];
        else
            u[i] = 0.0;
    }

    // 转换到局部坐标系
    // 计算内力
    // 详细实现...

    stress[0] = 0.0;  // 轴力
    stress[1] = 0.0;  // 剪力Vy
    stress[2] = 0.0;  // 剪力Vz
    stress[3] = 0.0;  // 弯矩My
    stress[4] = 0.0;  // 弯矩Mz
    stress[5] = 0.0;  // 扭矩T
}
```

**注意事项**:
- 梁单元需要6个自由度/节点（3平动+3转角）
- 需要坐标转换矩阵
- 刚度矩阵是12×12
- 应力输出包括内力和弯矩

---

### 步骤4: 添加梁材料类

**文件**: `src/h/Material.h`

在文件末尾添加：

```cpp
//! Material class for beam element
class CBeamMaterial : public CMaterial
{
public:

    double Area;    //!< Cross-sectional area
    double Iy;      //!< Moment of inertia about y-axis
    double Iz;      //!< Moment of inertia about z-axis
    double J;       //!< Torsional constant
    double G;       //!< Shear modulus

public:

    //! Read material data from stream Input
    virtual bool Read(ifstream& Input);

    //! Write material data to stream
    virtual void Write(COutputter& output);
};
```

**文件**: `src/cpp/Material.cpp`

添加实现：

```cpp
//! Read beam material data
bool CBeamMaterial::Read(ifstream& Input)
{
    Input >> nset;    // Material set number
    Input >> E;       // Young's modulus
    Input >> Area;    // Cross-sectional area
    Input >> Iy;      // Moment of inertia about y
    Input >> Iz;      // Moment of inertia about z
    Input >> J;       // Torsional constant
    Input >> G;       // Shear modulus

    return true;
}

//! Write beam material data
void CBeamMaterial::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << Area
           << setw(16) << Iy
           << setw(16) << Iz
           << setw(16) << J
           << setw(16) << G << endl;
}
```

---

### 步骤5: 更新ElementGroup类

**文件**: `src/h/ElementGroup.h`

在开头添加include：

```cpp
#include "Beam.h"  // 添加这一行
```

**文件**: `src/cpp/ElementGroup.cpp`

修改三个函数：

#### 5.1 CalculateMemberSize()

```cpp
void CElementGroup::CalculateMemberSize()
{
    switch (ElementType_)
    {
        case ElementTypes::UNDEFINED:
            std::cerr << "Setting element type to UNDEFINED." << std::endl;
            exit(5);
        case ElementTypes::Bar:
            ElementSize_ = sizeof(CBar);
            MaterialSize_ = sizeof(CBarMaterial);
            break;
        case ElementTypes::Beam:  // 添加这个case
            ElementSize_ = sizeof(CBeam);
            MaterialSize_ = sizeof(CBeamMaterial);
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available." << std::endl;
            exit(5);
            break;
    }
}
```

#### 5.2 AllocateElements()

```cpp
void CElementGroup::AllocateElements(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            ElementList_ = new CBar[size];
            break;
        case ElementTypes::Beam:  // 添加这个case
            ElementList_ = new CBeam[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available." << std::endl;
            exit(5);
    }
}
```

#### 5.3 AllocateMaterials()

```cpp
void CElementGroup::AllocateMaterials(std::size_t size)
{
    switch(ElementType_)
    {
        case ElementTypes::Bar:
            MaterialList_ = new CBarMaterial[size];
            break;
        case ElementTypes::Beam:  // 添加这个case
            MaterialList_ = new CBeamMaterial[size];
            break;
        default:
            std::cerr << "Type " << ElementType_ << " not available." << std::endl;
            exit(5);
    }
}
```

---

### 步骤6: 添加输出函数

**文件**: `src/h/Outputter.h`

在public部分添加声明：

```cpp
//! Output beam element data
void OutputBeamElements(unsigned int EleGrp);
```

**文件**: `src/cpp/Outputter.cpp`

#### 6.1 修改OutputElementInfo()

在switch语句中添加：

```cpp
void COutputter::OutputElementInfo()
{
    // ... 前面的代码 ...

    switch (ElementType)
    {
        case ElementTypes::Bar:
            OutputBarElements(EleGrp);
            break;
        case ElementTypes::Beam:  // 添加这个case
            OutputBeamElements(EleGrp);
            break;
        default:
            *this << ElementType << " has not been implemented yet." << endl;
            break;
    }
}
```

#### 6.2 实现OutputBeamElements()

```cpp
//! Output beam element data
void COutputter::OutputBeamElements(unsigned int EleGrp)
{
    CDomain* FEMData = CDomain::GetInstance();
    CElementGroup& ElementGroup = FEMData->GetEleGrpList()[EleGrp];
    unsigned int NUMMAT = ElementGroup.GetNUMMAT();

    *this << " M A T E R I A L   D E F I N I T I O N" << endl << endl;
    *this << " NUMBER OF DIFFERENT SETS OF MATERIAL" << endl;
    *this << " AND CROSS-SECTIONAL  CONSTANTS  . . . .( NPAR(3) ) . . ="
          << setw(5) << NUMMAT << endl << endl;

    *this << "  SET       YOUNG'S     CROSS-SECTIONAL    MOMENT OF      MOMENT OF      TORSIONAL     SHEAR" << endl
          << " NUMBER     MODULUS          AREA         INERTIA-Y      INERTIA-Z       CONSTANT     MODULUS" << endl
          << "               E              A               Iy             Iz              J            G" << endl;

    *this << setiosflags(ios::scientific) << setprecision(5);

    for (unsigned int mset = 0; mset < NUMMAT; mset++)
    {
        *this << setw(5) << mset+1;
        ElementGroup.GetMaterial(mset).Write(*this);
    }

    *this << endl << endl << " E L E M E N T   I N F O R M A T I O N" << endl;
    *this << " ELEMENT     NODE     NODE       MATERIAL" << endl
          << " NUMBER-N      I        J       SET NUMBER" << endl;

    unsigned int NUME = ElementGroup.GetNUME();

    for (unsigned int Ele = 0; Ele < NUME; Ele++)
    {
        *this << setw(5) << Ele+1;
        ElementGroup[Ele].Write(*this);
    }

    *this << endl;
}
```

#### 6.3 修改OutputElementStress()

在switch语句中添加：

```cpp
void COutputter::OutputElementStress()
{
    // ... 前面的代码 ...

    switch (ElementType)
    {
        case ElementTypes::Bar:
            // ... Bar的处理 ...
            break;

        case ElementTypes::Beam:  // 添加这个case
            *this << "  ELEMENT        AXIAL FORCE      SHEAR-Y      SHEAR-Z      MOMENT-Y     MOMENT-Z      TORQUE" << endl
                  << "  NUMBER" << endl;

            double stress[6];  // 6个内力分量

            for (unsigned int Ele = 0; Ele < NUME; Ele++)
            {
                CElement& Element = EleGrp[Ele];
                Element.ElementStress(stress, Displacement);

                *this << setw(5) << Ele + 1;
                for (int i = 0; i < 6; i++)
                    *this << setw(16) << stress[i];
                *this << endl;
            }

            *this << endl;
            break;

        default:
            cerr << "*** Error *** Element type " << ElementType
                 << " has not been implemented.\n\n";
    }
}
```

---

### 步骤7: 更新CMakeLists.txt

**文件**: `src/CMakeLists.txt`

添加新的源文件：

```cmake
set(SOURCES
    cpp/main.cpp
    cpp/Domain.cpp
    cpp/Node.cpp
    cpp/Bar.cpp
    cpp/Beam.cpp          # 添加这一行
    cpp/Material.cpp
    cpp/ElementGroup.cpp
    cpp/Solver.cpp
    cpp/LoadCaseData.cpp
    cpp/Outputter.cpp
    cpp/Clock.cpp
)
```

---

### 步骤8: 修改节点自由度（如果需要）

**文件**: `src/h/Node.h`

如果梁单元需要转角自由度，需要修改：

```cpp
class CNode
{
public:
    // 从3改为6，包含3个平动和3个转角
    const static unsigned int NDF = 6;  // 修改这里

    // ... 其他代码不变 ...
};
```

**注意**: 修改NDF会影响所有单元类型，需要确保Bar单元也能正确处理。

---

## 输入文件格式

添加梁单元后，输入文件格式示例：

```
Example Beam Structure
4  2  1  1
(NUMNP=4, NUMEG=2, NLCASE=1, MODEX=1)

节点数据:
1  0.0  0.0  0.0  1 1 1 1 1 1
2  5.0  0.0  0.0  0 0 0 0 0 0
3  10.0 0.0  0.0  0 0 0 0 0 0
4  10.0 5.0  0.0  0 0 0 0 0 0

载荷数据:
1  2
4  2  -1000.0
4  3  -500.0

单元组1 (Bar):
1  2  1
(ElementType=1, NUME=2, NUMMAT=1)
1  2.1e11  0.01
(nset=1, E, Area)
1  1  2  1
2  2  3  1

单元组2 (Beam):
5  1  1
(ElementType=5, NUME=1, NUMMAT=1)
1  2.1e11  0.01  8.33e-6  8.33e-6  1.67e-5  8.08e10
(nset=1, E, Area, Iy, Iz, J, G)
1  3  4  1
```

---

## 编译和测试

### 编译步骤

```bash
cd build
cmake ..
make
```

### 测试步骤

1. 创建测试输入文件 `beam_test.dat`
2. 运行程序：
   ```bash
   ./stap++ beam_test
   ```
3. 检查输出文件 `beam_test.out`

---

## 完整修改文件清单总结

| 序号 | 文件路径 | 操作 | 说明 |
|------|---------|------|------|
| 1 | `src/h/Beam.h` | **新建** | 梁单元头文件 |
| 2 | `src/cpp/Beam.cpp` | **新建** | 梁单元实现 |
| 3 | `src/h/Material.h` | **修改** | 添加CBeamMaterial类声明 |
| 4 | `src/cpp/Material.cpp` | **修改** | 实现CBeamMaterial方法 |
| 5 | `src/h/ElementGroup.h` | **修改** | 添加#include "Beam.h" |
| 6 | `src/cpp/ElementGroup.cpp` | **修改** | 3个函数添加Beam case |
| 7 | `src/h/Outputter.h` | **修改** | 声明OutputBeamElements() |
| 8 | `src/cpp/Outputter.cpp` | **修改** | 实现输出函数，修改2个switch |
| 9 | `src/CMakeLists.txt` | **修改** | 添加Beam.cpp到SOURCES |
| 10 | `src/h/Node.h` | **可选修改** | 如需转角自由度，修改NDF |

---

## 其他单元类型

### 添加板单元 (Plate)

类似步骤，关键差异：
- 节点数: 4 (矩形板) 或 3 (三角形板)
- 自由度: 每节点3个 (w, θx, θy)
- 材料参数: E, ν, 厚度t

### 添加壳单元 (Shell)

类似步骤，关键差异：
- 节点数: 4 或 8
- 自由度: 每节点6个 (3平动 + 3转角)
- 材料参数: E, ν, 厚度t
- 需要膜应力和弯曲应力

### 添加实体单元 (H8)

类似步骤，关键差异：
- 节点数: 8 (六面体)
- 自由度: 每节点3个 (u, v, w)
- 材料参数: E, ν
- 需要数值积分 (Gauss积分)

---

## 调试建议

1. **逐步测试**: 先实现Read和Write，确保数据读取正确
2. **简单算例**: 使用手算可验证的简单算例
3. **对比验证**: 与商业软件（ANSYS, ABAQUS）结果对比
4. **单元测试**: 为每个函数编写单元测试
5. **调试输出**: 使用`#ifdef _DEBUG_`输出中间结果

---

## 常见问题

### Q1: 编译错误 "undefined reference to CBeam"
**A**: 检查CMakeLists.txt是否添加了Beam.cpp

### Q2: 运行时崩溃
**A**: 检查LocationMatrix_和ND_的大小是否匹配

### Q3: 刚度矩阵奇异
**A**: 检查边界条件是否充分约束刚体位移

### Q4: 结果不收敛
**A**: 检查单元刚度矩阵是否对称正定

---

## 参考资料

1. **有限元方法基础** - 王勖成
2. **The Finite Element Method** - O.C. Zienkiewicz
3. **STAP90用户手册** - 清华大学
4. **STAPpp源代码** - 现有Bar单元实现

---

## 总结

添加新单元类型的核心步骤：

1. ✅ 创建单元类（继承CElement）
2. ✅ 创建材料类（继承CMaterial）
3. ✅ 更新ElementGroup（3个函数）
4. ✅ 更新Outputter（2个函数）
5. ✅ 更新CMakeLists.txt
6. ✅ 测试验证

遵循这个流程，可以系统地扩展STAPpp支持任意新单元类型。
