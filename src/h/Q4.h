/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Four-node bilinear quadrilateral element for 2D plane stress/strain
class CQ4 : public CElement
{
private:
    //! Build the 3x3 constitutive matrix for plane stress or plane strain
    void ConstitutiveMatrix(double D[3][3]) const;

    //! Evaluate shape-function derivatives and Jacobian at natural coordinates
    void StrainDisplacementMatrix(double xi, double eta, double B[3][8], double& detJ) const;

public:
    //! Constructor
    CQ4();

    //! Destructor
    ~CQ4();

    //! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

    //! Write element data to stream
    virtual void Write(COutputter& output);

    //! Generate the Q4 location matrix using only x/y translational DOFs
    virtual void GenerateLocationMatrix();

    //! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

    //! Calculate Gauss-point stresses: 4 * [sigma_x, sigma_y, tau_xy]
    virtual void ElementStress(double* stress, double* Displacement);
};
