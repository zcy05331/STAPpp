/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*****************************************************************************/

#include "Q4.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

CQ4::CQ4()
{
    NEN_ = 4;
    nodes_ = new CNode*[NEN_];

    ND_ = 8;
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    unsigned int N1, N2, N3, N4;

    Input >> N1 >> N2 >> N3 >> N4 >> MSet;

    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    nodes_[2] = &NodeList[N3 - 1];
    nodes_[3] = &NodeList[N4 - 1];

    for (unsigned int i = 0; i < NEN_; i++)
    {
        if (nodes_[i]->bcode[2] != 0)
        {
            cerr << "*** Error *** Q4 element requires the z DOF of every connected node to be constrained." << endl
                 << "    Node number = " << nodes_[i]->NodeNumber << endl;
            return false;
        }
    }

    return true;
}

void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
    {
        LocationMatrix_[i++] = nodes_[N]->bcode[0];
        LocationMatrix_[i++] = nodes_[N]->bcode[1];
    }
}

void CQ4::ConstitutiveMatrix(double D[3][3]) const
{
    clear(&D[0][0], 9);

    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
    double E = material->E;
    double nu = material->nu;

    if (material->analysisType == 0) // plane stress
    {
        double c = E / (1.0 - nu * nu);
        D[0][0] = c;
        D[0][1] = c * nu;
        D[1][0] = c * nu;
        D[1][1] = c;
        D[2][2] = c * (1.0 - nu) / 2.0;
    }
    else // plane strain
    {
        double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D[0][0] = c * (1.0 - nu);
        D[0][1] = c * nu;
        D[1][0] = c * nu;
        D[1][1] = c * (1.0 - nu);
        D[2][2] = c * (1.0 - 2.0 * nu) / 2.0;
    }
}

void CQ4::StrainDisplacementMatrix(double xi, double eta, double B[3][8], double& detJ) const
{
    clear(&B[0][0], 24);

    double dN_dxi[4];
    double dN_deta[4];

    dN_dxi[0] = -0.25 * (1.0 - eta);
    dN_dxi[1] =  0.25 * (1.0 - eta);
    dN_dxi[2] =  0.25 * (1.0 + eta);
    dN_dxi[3] = -0.25 * (1.0 + eta);

    dN_deta[0] = -0.25 * (1.0 - xi);
    dN_deta[1] = -0.25 * (1.0 + xi);
    dN_deta[2] =  0.25 * (1.0 + xi);
    dN_deta[3] =  0.25 * (1.0 - xi);

    double J[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    for (unsigned int i = 0; i < NEN_; i++)
    {
        double x = nodes_[i]->XYZ[0];
        double y = nodes_[i]->XYZ[1];
        J[0][0] += dN_dxi[i] * x;
        J[0][1] += dN_deta[i] * x;
        J[1][0] += dN_dxi[i] * y;
        J[1][1] += dN_deta[i] * y;
    }

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    if (detJ <= 0.0)
    {
        cerr << "*** Error *** Invalid Q4 element geometry: detJ <= 0 at a Gauss point." << endl
             << "    detJ = " << detJ << endl;
        exit(6);
    }

    double invJ[2][2];
    invJ[0][0] =  J[1][1] / detJ;
    invJ[0][1] = -J[0][1] / detJ;
    invJ[1][0] = -J[1][0] / detJ;
    invJ[1][1] =  J[0][0] / detJ;

    for (unsigned int i = 0; i < NEN_; i++)
    {
        double dN_dx = invJ[0][0] * dN_dxi[i] + invJ[0][1] * dN_deta[i];
        double dN_dy = invJ[1][0] * dN_dxi[i] + invJ[1][1] * dN_deta[i];

        unsigned int col = 2 * i;
        B[0][col]     = dN_dx;
        B[1][col + 1] = dN_dy;
        B[2][col]     = dN_dy;
        B[2][col + 1] = dN_dx;
    }
}

void CQ4::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());

    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);

    double K[8][8];
    clear(&K[0][0], 64);

    double D[3][3];
    ConstitutiveMatrix(D);

    const double g = 1.0 / sqrt(3.0);
    const double gps[4][2] = {{-g, -g}, {g, -g}, {g, g}, {-g, g}};

    for (unsigned int gp = 0; gp < 4; gp++)
    {
        double B[3][8];
        double detJ;
        StrainDisplacementMatrix(gps[gp][0], gps[gp][1], B, detJ);

        for (unsigned int i = 0; i < ND_; i++)
        {
            for (unsigned int j = 0; j < ND_; j++)
            {
                double v = 0.0;
                for (unsigned int a = 0; a < 3; a++)
                    for (unsigned int b = 0; b < 3; b++)
                        v += B[a][i] * D[a][b] * B[b][j];

                K[i][j] += v * detJ * material->thickness;
            }
        }
    }

    for (unsigned int j = 0; j < ND_; j++)
        for (unsigned int i = 0; i <= j; i++)
            Matrix[(j + 1) * j / 2 + j - i] = K[i][j];
}

void CQ4::ElementStress(double* stress, double* Displacement)
{
    double D[3][3];
    ConstitutiveMatrix(D);

    double u[8];
    clear(u, 8);
    for (unsigned int i = 0; i < ND_; i++)
    {
        if (LocationMatrix_[i])
            u[i] = Displacement[LocationMatrix_[i] - 1];
    }

    const double g = 1.0 / sqrt(3.0);
    const double gps[4][2] = {{-g, -g}, {g, -g}, {g, g}, {-g, g}};

    for (unsigned int gp = 0; gp < 4; gp++)
    {
        double B[3][8];
        double detJ;
        StrainDisplacementMatrix(gps[gp][0], gps[gp][1], B, detJ);

        double strain[3] = {0.0, 0.0, 0.0};
        for (unsigned int a = 0; a < 3; a++)
            for (unsigned int i = 0; i < ND_; i++)
                strain[a] += B[a][i] * u[i];

        for (unsigned int a = 0; a < 3; a++)
        {
            stress[3 * gp + a] = 0.0;
            for (unsigned int b = 0; b < 3; b++)
                stress[3 * gp + a] += D[a][b] * strain[b];
        }
    }
}
