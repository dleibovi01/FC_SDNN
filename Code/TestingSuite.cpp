
#include "TestingSuite.h"
#include "mkl.h"

const std::string B1_filename = "../include/B1.txt";
const std::string B2_filename = "../include/B2.txt";
const std::string B3_filename = "../include/B3.txt";
const std::string B4_filename = "../include/B4.txt";
const std::string W1_filename = "../include/W1.txt";
const std::string W2_filename = "../include/W2.txt";
const std::string W3_filename = "../include/W3.txt";
const std::string W4_filename = "../include/W4.txt"; 

void TestingVector1D()
{
    int N = 10;
    double a [N];
    double b [N];
    double c [N];
    for(int i = 0; i < N; i++)
    {
        a[i] = double(i);
        b[i] = double(N) - double(i) - 1.0;
        c[i] = 5.0;
    }
    std::vector<double*> field;
    field.push_back(a);
    field.push_back(b);
    field.push_back(c);



    VectorField1D e1{field, N};
    Print_VectorField1D(e1);
    VectorField1D e2{e1};
    a[0] = 100;
    e1.setField(0, N, b);
    Print_VectorField1D(e1);

    VectorField1D e3 = e1;
    std::cout << "e3" << std::endl;
    Print_VectorField1D(e3);

    e1.setField(1, N, a);
    std::cout << "e3" << std::endl;
    Print_VectorField1D(e3);



    e3 = e2;
    std::cout << "e3" << std::endl;
    Print_VectorField1D(e3);    

    std::cout << "e1" << std::endl;
    Print_VectorField1D(e1);       

    std::cout << "Sum" << std::endl;
    Print_VectorField1D(e1 + e3);   

    std::cout << "Product" << std::endl;
    Print_VectorField1D(e1 * e3);     

    std::cout << "Product scalar" << std::endl;
    Print_VectorField1D(e3*=2.0);     

    std::cout << "Product scalar" << std::endl;
    Print_VectorField1D(0.1*e3);       

/*
    std::cout << "Circshifted" << std::endl;
    e3.circshift(1);
    Print_VectorField1D(e3);  
*/
    std::cout << "Circshifted function" << std::endl;
    VectorField1D e4 = circshift(e3, -1);
    Print_VectorField1D(e4);  
}

void TestingNode1D()
{   
  std::vector<double> a;
    a.push_back(1.0);
    a.push_back(3.5);
    a.push_back(-1.4);  
   std::vector<double> b;
    b.push_back(5.0);
    b.push_back(2.5);
    b.push_back(-1.1);
    Node1D n{0.22, 0, 3, a};
    Print_Node1D(n);
    n.setValue(1, -1.0);
    Print_Node1D(n);
    n.setValues(b);
    Print_Node1D(n);
    Node1D n2{4};
    Print_Node1D(n2);
}

void TestingPatch1D()
{
    double a = 0.0;
    double b = 1.0;
    int N = 10;
    int unknowns = 3;
    bool lb = true;
    bool rb = false; 
    int intrbl = 0; 
    int intrbr = 0;
    Patch1DUniform patch{N, unknowns, a, b, rb, lb, intrbl, intrbr};
    Print_Patch1D(patch);

    double a1 [N];
    double b1 [N];
    double c1 [N];
    for(int i = 0; i < N; i++)
    {
        a1[i] = double(i);
        b1[i] = double(N) - double(i) - 1.0;
        c1[i] = 5.0;
    }
    std::vector<double*> field;
    field.push_back(a1);
    field.push_back(b1);
    field.push_back(c1);


    VectorField1D e1{field, N};
    patch.setField(e1);

    patch.VectorFieldToNodes();
    Print_Patch1D(patch);

    patch.getNode(5)->setValue(2, 1.0);
    Print_Patch1D(patch);

    Patch1DUniform patch2{N, 1, a, b, rb, lb, intrbl, intrbr};
    std::string problem = "LA_Gaussian"; 
    // std::string problem = ""; 

    Print_Patch1D(patch2);
    IC ic{problem, 1};
    // ic(&patch2);
    ic.setIC(&patch2);
    Print_Patch1D(patch2);

}

void TestingMesh1D()
{
    Mesh1DUniform mesh{0, 1, 2, 10, 2, 1, 1, true, false};
    Print_Mesh1D(mesh);
    std::string problem = "step"; 
    IC ic{problem, 1};
    ic(&mesh);
    Print_Mesh1D(mesh);
    mesh.setIntraPatchBC(1);
    Print_Mesh1D(mesh);
    int unknowns = 1;
    BC bc{"", unknowns};
    bc(&mesh, 0.0);
    Print_Mesh1D(mesh);

    Mesh1DUniform mesh2 = mesh;
    std::cout <<"Printing second mesh" << std::endl;
    Print_Mesh1D(mesh2);

}

void TestingSpatDiffScheme()
{
    /*
    FD1_1D diff_scheme;
    int N = 10;
    double a [N];
    double b [N];
    double c [N];
    double h = 0.1;
    for(int i = 0; i < N; i++)
    {
        a[i] = double(i);
        b[i] = double(N) - double(i) - 1.0;
        c[i] = 5.0;
    }
    std::vector<double*> field;
    field.push_back(a);
    field.push_back(b);
    field.push_back(c);
    VectorField1D e1{field, N}; 

    VectorField1D e2 = diff_scheme.diff(e1, h);

    std::cout << "Initial vectorfield" << std::endl;
    Print_VectorField1D(e1);


    std::cout << "Differentiated vectorfield" << std::endl;
    Print_VectorField1D(e2); 
    */   
}

void TestingTimeStepScheme()
{
    /*
    FD1_1D<Mesh1DUniform> diff_scheme;
    int N = 10;
    double a [N];
    double b [N];
    double c [N];
    double h = 0.1;
    double dt = 0.1;
    Fwd_Euler<FD1_1D<Mesh1DUniform> > TS{diff_scheme, dt};
    
    
    for(int i = 0; i < N; i++)
    {
        a[i] = double(i);
        b[i] = double(N) - double(i) - 1.0;
        c[i] = 5.0;
    }
    std::vector<double*> field;
    field.push_back(a);
    field.push_back(b);
    field.push_back(c);
    VectorField1D e1{field, N}; 

    VectorField1D e2 = TS.advance(e1);
global std::string B1_filename = "B1.txt";
    Print_VectorField1D(e1);


    std::cout << "Differentiated vectorfield" << std::endl;
    Print_VectorField1D(e2);    
    */
}

void TestingSDNN()
{  

    // Create a PDE
    double T = 0.7;
    // double dt = 0.00002;
    double dt = 0.001;
    // std::string problem = "smooth_LA";
    std::string problem = "one_wave";
    BC bc{problem, 1};
    IC ic{problem, 1};
    double alpha = 1.0;
    double a = 1.0;
    
    // ANN ann{W1_filename, W2_filename, W3_filename, W4_filename, B1_filename, 
    //     B2_filename, B3_filename, B4_filename, alpha};  
    ANN ann{alpha};  
    LA1D_SDNN pde {ic, bc, T, a, ann};      
    int intrb = 6;
    int npatches = 1;
    int overlap = 24;
    // int N = 250 + (npatches - 1)*overlap;
    int N = 100;


   // Create a time-stepping scheme
    SSPRK_4 TS;
    int stages = TS.getStages();
    

    // Create a mesh
    int unknowns = pde.getPhysUnknowns();
    Mesh1DUniform mesh{0, 1.4, npatches, N, overlap, intrb, unknowns*stages + 3, 
        true, false};
    ic(&mesh);
    //Print_Mesh1D(mesh);   

    VectorField1D v = (mesh.getPatches())[0]->getFlow();
    int C = 27;
    int d = 5;

    double* A = new double [C*d];
    double* Q = new double [d*d];
    double* AQ = new double [C*d];
    double* FAQF = new double [C*d];    
    std::string filename_A = "../include/Ad5C27.txt";
    std::string filename_Q = "../include/Qd5C27.txt";
    read_FC_Data(A, Q, d, C, filename_A, filename_Q);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);

    int stage = 0;
    std::vector<FC_1D<Patch1DUniform> > diff_schemes;
    std::vector<FC_1D<Patch1DUniform> > filters;
    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;
    double delta = 0.1;
    for(int i = 0; i < npatches; i++)
    {
        diff_schemes.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], delta));
        filters.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], alpha0, p_0));
    }
    // Print_Mat(diff_schemes[0].getShiftCoeffs(), N + C, 1);
    for(int i = 0; i < npatches; i++)
    {
        filters.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], alpha0, p));
    }


    // Create a solver
    Solver<Mesh1DUniform, LA1D_SDNN, SSPRK_4, FC_1D<Patch1DUniform>, 
        FC_1D<Patch1DUniform> > slv{mesh, pde, TS, diff_schemes, filters};

    // Run the solver
    double CFL = 1.5;
    bool visc = true;
    bool adaptive = false;
    slv.solve_sdnn(dt, CFL, adaptive, visc);    

    // Print the solution
    Mesh1DUniform mesh1 = slv.getMesh();

    std::cout << "Solution" << std::endl;
    Print_Mesh1D(mesh1);

    std::string result_file = "result.txt";
    Print_Mesh1D(mesh1, unknowns, result_file);
}


void TestingEulerSDNN()
{  

    // Create a PDE
    double T = 0.2;
    // double dt = 0.00002;
    double dt = 0.000005;
    // std::string problem = "smooth_LA";
    std::string problem = "Euler1D_Sod";
    BC bc{problem, 3};
    IC ic{problem, 3};
    double alpha = 1.0;
    double gamma = 7.0/5.0;
    
    // ANN ann{W1_filename, W2_filename, W3_filename, W4_filename, B1_filename, 
    //     B2_filename, B3_filename, B4_filename, alpha};  
    ANN ann{alpha};  
    Euler1D_SDNN pde {ic, bc, T, gamma, ann};      
    int intrb = 6;
    int npatches = 1;
    int overlap = 24;
    // int N = 250 + (npatches - 1)*overlap;
    int N = 100;


   // Create a time-stepping scheme
    SSPRK_4 TS;
    int stages = TS.getStages();
    

    // Create a mesh
    int unknowns = pde.getPhysUnknowns();
    Mesh1DUniform mesh{0, 1.0, npatches, N, overlap, intrb, unknowns*stages + 3, 
        true, true};
    ic(&mesh);
    //Print_Mesh1D(mesh);   

    VectorField1D v = (mesh.getPatches())[0]->getFlow();
    int C = 27;
    int d = 5;

    double* A = new double [C*d];
    double* Q = new double [d*d];
    double* AQ = new double [C*d];
    double* FAQF = new double [C*d];    
    std::string filename_A = "../include/Ad5C27.txt";
    std::string filename_Q = "../include/Qd5C27.txt";
    read_FC_Data(A, Q, d, C, filename_A, filename_Q);
    build_Cont_Mat(A, Q, d, C, AQ, FAQF);

    int stage = 0;
    std::vector<FC_1D<Patch1DUniform> > diff_schemes;
    std::vector<FC_1D<Patch1DUniform> > filters;
    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;
    double delta = 0.1;
    for(int i = 0; i < npatches; i++)
    {
        diff_schemes.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], delta));
        filters.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], alpha0, p_0));
    }
    // Print_Mat(diff_schemes[0].getShiftCoeffs(), N + C, 1);
    for(int i = 0; i < npatches; i++)
    {
        filters.push_back(FC_1D<Patch1DUniform> 
            (N, d, C, mesh.getPatches()[i], alpha0, p));
    }


    // Create a solver
    Solver<Mesh1DUniform, Euler1D_SDNN, SSPRK_4, FC_1D<Patch1DUniform>, 
        FC_1D<Patch1DUniform> > slv{mesh, pde, TS, diff_schemes, filters};

    // Run the solver
    double CFL = 2.0;
    bool visc = true;
    bool adaptive = false;
    slv.solve_sdnn(dt, CFL, adaptive, visc);    

    // Print the solution
//     Mesh1DUniform mesh1 = slv.getMesh();

//     std::cout << "Solution" << std::endl;
//     Print_Mesh1D(mesh1);

//     std::string result_file = "result.txt";
//     Print_Mesh1D(mesh1, unknowns, result_file);
}






void TestingSolver()
{
    /*
    // Create a mesh
    int N = 100;
    Mesh1DUniform mesh{0, 1, 1, N, 2, 1};
    //Print_Mesh1D(mesh);

    // Create a PDE
    double T = 0.005;
    BC bc{"LA_test"};
    IC ic{"LA_test"};
    double a = 1.0;
    LA1D pde{ic, bc, T, a};

    // Create a spatial differentiation scheme
    int d = 5;
    int C = 27;
    std::string filename_A = "Ad5C27.txt";
    std::string filename_Q = "Qd5C27.txt";
    double h = mesh.getH();
    FC_1D<Mesh1DUniform> diff_scheme{N, d, C, &mesh, filename_A, filename_Q};

    // Create a time-stepping scheme
    double dt = 0.000001;
    Fwd_Euler TS{dt}; 
    //Fwd_Euler TS{dt};

    // Create a solver
    Solver<Mesh1DUniform, LA1D, Fwd_Euler, FC_1D<Mesh1DUniform> > 
    slv{mesh, pde, TS, diff_scheme};

    // Run the solver
    slv.solve(dt);

    // Print the solution
    Mesh1DUniform mesh1 = slv.getMesh();
    std::cout << "Solution" << std::endl;
    Print_Mesh1D(mesh1);
    */
}

void TestingSVW()
{
    int N = 20;
    int npatches = 2;
    Mesh1DUniform mesh{0., 1., npatches, N, 2, 1, 1, true, false};
    std::string problem = "step"; 
    IC ic{problem, 1};
    auto patches = mesh.getPatches();
    // int npatches = patches.size();
    std::cout << "npatches = " << npatches << std::endl;
    SVW_mesh svw_m{mesh};
    auto v = svw_m.getSVWs()[0]->getPatchSVWS();
    std::cout << "number of cols = " << v[1]->getCols() << std::endl;
    std::cout << "number of rows = " << v[1]->getRows() << std::endl;
    v[0]->print();
    std::vector<double> sum_vect(N, 0.0);
    double* sum = sum_vect.data();
    v[0]->rowSum(sum);
    Print_Mat(sum, N, 1);

    std::cout << std::endl;
}


int main()
{
    
    //TestingVector1D();
    //TestingSpatDiffScheme();
    //TestingTimeStepScheme();
    //TestingNode1D();
    //TestingPatch1D();
    //TestingMesh1D();
    // TestingSDNN();
    TestingEulerSDNN();
    // TestingSVW();

    
    // auto t1 = std::chrono::high_resolution_clock::now();

    // TestingSDNN();

    // auto t2 = std::chrono::high_resolution_clock::now();

    // // Getting number of milliseconds as an integer. //
    // auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
    // // Getting number of milliseconds as a double. //
    // std::chrono::duration<double, std::milli> ms_double = t2 - t1;

    // std::cout << ms_int.count() << "ms\n" << std::endl;
    // std::cout << ms_double.count() << "ms" << std::endl;        
    


    return 0;
}




