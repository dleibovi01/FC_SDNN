
#ifndef TESTINGSUITE_H
#define TESTINGSUITE_H

#include <iostream>
#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include <vector>
#include "VectorField.h"
#include "Node.h"
#include "Patch1D.h"
#include "Mesh.h"
#include "Mesh1D.h"
#include "IC.h"
#include "BC_LA_1D.h"
#include "BC_Euler_1D.h"

#include "BC.h"
#include "printing.h"
#include "SpatDiffScheme.h"
#include "TimeStepScheme.h"
#include "Solver.h"
#include <chrono>
#include "SVW.h"
#include "SpMatrix_csr.h"

#include "TestingMesh.h"
#include "TestingFC.h"
#include "TestingSDNN.h"
#include "TestingSolver.h"


void TestingVector1D();
void TestingNode();
void TestingPatch1D();
void TestingMesh1D();
void TestingSpatDiffScheme();
void TestingTimeStepScheme();
void TestingSolver();
void TestingSVW();
void TestingSDNN();
void TestingEulerSDNN();
void TestingEulerWENO();
void TestingVML();

#endif
