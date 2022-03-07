
#ifndef TESTINGSUITE_H
#define TESTINGSUITE_H

#include <iostream>
#include <complex>
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include <vector>
#include "VectorField1D.h"
#include "Node1D.h"
#include "Patch1D.h"
#include "Mesh.h"
#include "Patch1DUniform.h"
#include "IC.h"
#include "BC.h"
#include "printing.h"
#include "SpatDiffScheme.h"
#include "TimeStepScheme.h"
#include "Solver.h"
#include <chrono>
#include "SVW.h"
#include "SpMatrix_csr.h"


void TestingVector1D();
void TestingNode1D();
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
