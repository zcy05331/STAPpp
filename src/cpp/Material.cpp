/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}


// Read material data for Q4 plane stress/strain element
bool CQ4Material::Read(ifstream& Input)
{
    Input >> nset;
    Input >> E >> nu >> thickness >> analysisType;

    if (E <= 0.0)
    {
        cerr << "*** Error *** Q4 material Young's modulus must be positive." << endl;
        exit(5);
    }

    if (thickness <= 0.0)
    {
        cerr << "*** Error *** Q4 material thickness must be positive." << endl;
        exit(5);
    }

    if (nu <= -1.0 || nu >= 0.5)
    {
        cerr << "*** Error *** Q4 material Poisson's ratio must satisfy -1 < nu < 0.5." << endl;
        exit(5);
    }

    if (analysisType != 0 && analysisType != 1)
    {
        cerr << "*** Error *** Q4 analysis type must be 0 (plane stress) or 1 (plane strain)." << endl;
        exit(5);
    }

    return true;
}

// Write material data for Q4 plane stress/strain element
void CQ4Material::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu
           << setw(16) << thickness
           << setw(10) << analysisType << endl;
}
