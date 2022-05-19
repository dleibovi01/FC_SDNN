/* 2D Mesh */

#include "Mesh2D.h"
#include "printing.h"

Mesh2DUniform::Mesh2DUniform(double a, double b, double c, double d,
    int npatches_x, int npatches_y, int patchsize_x, int patchsize_y,
    int _overlap, int fringe, int unknowns, bool lb, bool rb, bool db, bool ub)
    : overlap{_overlap}, Npatches_x{npatches_x}, Npatches_y{npatches_y}  
{
    // overlap = _overlap;
    // Npatches_x = npatches_x;
    // Npatches_y = npatches_y;
    int patchsize = patchsize_x * patchsize_y;
    int Nx = npatches_x * patchsize_x - (npatches_x - 1) * overlap;
    int Ny = npatches_y * patchsize_y - (npatches_y - 1) * overlap;
    double hx = (b - a) / (double (Nx - 1));
    double hy = (d - c) / (double (Ny - 1));

    std::vector<int> inner_nodes_onepatch(patchsize);

    std::vector<int> inner_nodes_onlyleftpatch;
    std::vector<int> inner_nodes_onlyrightpatch;
    std::vector<int> inner_nodes_onlybottompatch;
    std::vector<int> inner_nodes_onlytoppatch;

    std::vector<int> inner_nodes_topleftpatch;
    std::vector<int> inner_nodes_toppatch;
    std::vector<int> inner_nodes_toprightpatch;
    std::vector<int> inner_nodes_leftpatch;
    std::vector<int> inner_nodes_rightpatch;
    std::vector<int> inner_nodes_bottomleftpatch;
    std::vector<int> inner_nodes_bottompatch;
    std::vector<int> inner_nodes_bottomrightpatch;
    std::vector<int> inner_nodes_innerpatch;

    // one patch
    for(int i = 0; i < patchsize_x; i++)
    {
        for(int j = 0; j < patchsize_y; j++)
        {
            inner_nodes_onepatch.push_back(patchsize_y * i + j);
        }
    }

    // Only a left and a right patch
    for(int i = 0; i < patchsize_x - fringe; i++)
    {
        for(int j = 0; j < patchsize_y; j++)
        {
            inner_nodes_onlyleftpatch.push_back(patchsize_y * i + j);
        }
    }
    for(int i = fringe; i < patchsize_x; i++)
    {
        for(int j = 0; j < patchsize_y; j++)
        {
            inner_nodes_onlyrightpatch.push_back(patchsize_y * i + j);
        }
    }

    // Only a bottom and a top patch
    for(int i = 0; i < patchsize_x; i++)
    {
        for(int j = 0; j < patchsize_y - fringe; j++)
        {
            inner_nodes_onlybottompatch.push_back(patchsize_y * i + j);
        }
    }
    for(int i = 0; i < patchsize_x; i++)
    {
        for(int j = fringe; j < patchsize_y; j++)
        {
            inner_nodes_onlytoppatch.push_back(patchsize_y * i + j);
        }
    }


    // top left inner nodes
    for(int i = 0; i < patchsize_x - fringe; i++)
    {
        for(int j = fringe; j < patchsize_y; j++)
        {
            inner_nodes_topleftpatch.push_back(patchsize_y * i + j);
        }
    }

    // top inner nodes
    for(int i = 0; i < patchsize_x; i++)
    {
        for(int j = fringe; j < patchsize_y; j++)
        {
            inner_nodes_toppatch.push_back(patchsize_y * i + j);
        }
    }

    // top right inner nodes
    for(int i = fringe; i < patchsize_x; i++)
    {
        for(int j = fringe; j < patchsize_y; j++)
        {
            inner_nodes_toprightpatch.push_back(patchsize_y * i + j);
        }
    }

    // left inner nodes
    for(int i = 0; i < patchsize_x - fringe; i++)
    {
        for(int j = fringe; j < patchsize_y - fringe; j++)
        {
            inner_nodes_leftpatch.push_back(patchsize_y * i + j);
        }
    }

    // middle inner nodes
    for(int i = fringe; i < patchsize_x - fringe; i++)
    {
        for(int j = fringe; j < patchsize_y - fringe; j++)
        {
            inner_nodes_innerpatch.push_back(patchsize_y * i + j);
        }
    }

    // right inner nodes
    for(int i = fringe; i < patchsize_x; i++)
    {
        for(int j = fringe; j < patchsize_y - fringe; j++)
        {
            inner_nodes_rightpatch.push_back(patchsize_y * i + j);
        }
    }

    // bottom left inner nodes
    for(int i = 0; i < patchsize_x - fringe; i++)
    {
        for(int j = 0; j < patchsize_y - fringe; j++)
        {
            inner_nodes_bottomleftpatch.push_back(patchsize_y * i + j);
        }
    }

    // bottom inner nodes
    for(int i = 0; i < patchsize_x; i++)
    {
        for(int j = 0; j < patchsize_y - fringe; j++)
        {
            inner_nodes_bottompatch.push_back(patchsize_y * i + j);
        }
    }

    // bottom right inner nodes
    for(int i = fringe; i < patchsize_x; i++)
    {
        for(int j = 0; j < patchsize_y - fringe; j++)
        {
            inner_nodes_bottomrightpatch.push_back(patchsize_y * i + j);
        }
    }

    bool lb_patch = false;
    bool rb_patch = false;
    bool db_patch = false;
    bool ub_patch = false;

    // std::vector<int> inner_nodes;

    for(int i = 0; i < npatches_x; i++)
    {
        if(i == 0)
            lb_patch = lb;               
        if(i == npatches_x - 1)
            rb_patch = rb;
        for(int j = 0; j < npatches_y; j++)
        {
            if(j == 0)
                db_patch = db;
            if(j == npatches_y - 1)
                ub_patch = ub;
            patches.push_back(new Patch2D(patchsize_x, patchsize_y, unknowns, 
                        a + double(i*(patchsize_x - overlap))*hx, 
                        a + double(i*(patchsize_x - overlap) + 
                        (patchsize_x - 1))*hx,
                        c + double(j*(patchsize_y - overlap))*hy,
                        c + double(j*(patchsize_y - overlap) + 
                        (patchsize_y - 1))*hy,
                        lb_patch, rb_patch, db_patch, ub_patch, fringe));            
            db_patch = false;
            ub_patch = false;
        }
        lb_patch = false;
        rb_patch = false;

        
    }

    // Inner nodes if there is only 1 patch vertically or horizontally ?

    if(npatches_x == 1 && npatches_y == 1)
    {
        patches[0]->setInnerNodes(inner_nodes_onepatch);
    }
    else if(npatches_y == 1 && npatches_x > 1)
    {
        patches[0]->setInnerNodes(inner_nodes_onlyleftpatch);
        patches[1]->setInnerNodes(inner_nodes_onlyrightpatch);
    }
    else if(npatches_x == 1 && npatches_y > 1)
    {
        patches[0]->setInnerNodes(inner_nodes_onlybottompatch);
        patches[1]->setInnerNodes(inner_nodes_onlytoppatch);
    }    
    else
    {
        // Gotta set the inner nodes now
        for(int i = 0; i < npatches_x; i++)
        {
            for(int j = 0; j < npatches_y; j++)
            {
                if(i == 0)
                {
                    if(j == 0)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_bottomleftpatch);
                    else if(j == npatches_y - 1)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_topleftpatch);
                    else
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_leftpatch);                             
                }
                else if(i == npatches_x - 1)
                {
                    if(j == 0)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_bottomrightpatch);    
                    else if(j == npatches_y - 1)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_toprightpatch);
                    else
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_rightpatch);               
                }
                else
                    if(j == 0)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_bottompatch);    
                    else if(j == npatches_y - 1)
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_toppatch);
                    else
                        patches[i*npatches_y + j]->
                            setInnerNodes(inner_nodes_innerpatch);  
            }
        }
    }



}


void Mesh2DUniform::setIntraPatchBC(int unknowns, int stage)
{
    // int N;
    // N = patches.size();
    int fringe = patches[0]->getFringe();

    // // If there is only one patch in the x-direction
    // if(Npatches_x == 1)
    // {
    //     // If there is only one patch in the y-direction
    //     if(Npatches_y == 1)
    //     {
    //         // there's only one patch, so we do nothing.
    //     }
    //     // If there are only two patches in the y-direction
    //     else if(Npatches_y == 2)
    //     {
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);
    //         patches[1]->setInnerBdry(patches[0], 2, fringe, overlap, unknowns,
    //             stage);
    //     }
    //     // else, if there are three or more patches in the y-direction
    //     {
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);           
    //         for(int i = 1; i < Npatches_y - 1; i++)
    //         {
    //             patches[i]->setInnerBdry(patches[i + 1], 8, fringe, overlap,
    //                 unknowns, stage);
    //             patches[i]->setInnerBdry(patches[i - 1], 2, fringe, overlap,
    //                 unknowns, stage);
    //         }
    //         patches[Npatches_y - 1]->setInnerBdry(patches[Npatches_y - 2], 2,
    //             fringe, overlap, unknowns, stage);
    //     }
    // }
    // // If there are only two patches in the x-direction
    // else if(Npatches_x == 2)
    // {
    //     // If there is only one patch in the y-direction
    //     if(Npatches_y == 1)
    //     {
    //         patches[0]->setInnerBdry(patches[1], 6, fringe, overlap, unknowns,
    //             stage);
    //         patches[1]->setInnerBdry(patches[0], 3, fringe, overlap, unknowns,
    //             stage);
    //     }
    //     // If there are only two patches in the y-direction
    //     else if(Npatches_y == 2)
    //     {
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);
    //         patches[0]->setInnerBdry(patches[2], 6, fringe, overlap, unknowns,
    //             stage);
    //         patches[1]->setInnerBdry(patches[0], 2, fringe, overlap, unknowns,
    //             stage);
    //         patches[1]->setInnerBdry(patches[3], 6, fringe, overlap, unknowns,
    //             stage);
    //         patches[2]->setInnerBdry(patches[0], 4, fringe, overlap, unknowns,
    //             stage);
    //         patches[2]->setInnerBdry(patches[3], 8, fringe, overlap, unknowns,
    //             stage);
    //         patches[3]->setInnerBdry(patches[2], 2, fringe, overlap, unknowns,
    //             stage);
    //         patches[3]->setInnerBdry(patches[1], 4, fringe, overlap, unknowns,
    //             stage);
    //         patches[0]->setInnerBdry(patches[3], 9, fringe, overlap, unknowns,
    //             stage);
    //         patches[1]->setInnerBdry(patches[2], 3, fringe, overlap, unknowns,
    //             stage);
    //         patches[2]->setInnerBdry(patches[1], 7, fringe, overlap, unknowns,
    //             stage);
    //         patches[3]->setInnerBdry(patches[0], 1, fringe, overlap, unknowns,
    //             stage);
    //     }
    //     // else, if there are three or more patches in the y-direction
    //     {
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);   
    //         patches[0]->setInnerBdry(patches[Npatches_y], 6, fringe, overlap,
    //             unknowns, stage);
    //         patches[0]->setInnerBdry(patches[Npatches_y + 1], 9, fringe,
    //             overlap, unknowns, stage);    
    //         patches[Npatches_y]->setInnerBdry(patches[Npatches_y + 1], 8,
    //             fringe, overlap, unknowns, stage);   
    //         patches[Npatches_y]->setInnerBdry(patches[0], 4, fringe, overlap,
    //             unknowns, stage);  
    //         patches[Npatches_y]->setInnerBdry(patches[1], 7, fringe, overlap,
    //             unknowns, stage);       
    //         for(int i = 1; i < Npatches_y - 1; i++)
    //         {
    //             patches[i]->setInnerBdry(patches[i - 1], 2, fringe, overlap,
    //                 unknowns, stage);   
    //             patches[i]->setInnerBdry(patches[Npatches_y + i], 6, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[i]->setInnerBdry(patches[i + 1], 8, fringe, overlap,
    //                 unknowns, stage);  
    //             patches[i]->setInnerBdry(patches[Npatches_y + i - 1], 3, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[i]->setInnerBdry(patches[Npatches_y + i + 1], 9, fringe,
    //                 overlap, unknowns, stage);    
    //             patches[Npatches_y + i]->setInnerBdry(
    //                 patches[Npatches_y + i - 1], 2, fringe, overlap, unknowns,
    //                 stage);   
    //             patches[Npatches_y + i]->setInnerBdry(
    //                 patches[Npatches_y + i + 1], 8, fringe, overlap, unknowns,
    //                 stage); 
    //             patches[Npatches_y + i]->setInnerBdry(patches[i], 4, fringe,
    //                 overlap, unknowns, stage);  
    //             patches[Npatches_y + i]->setInnerBdry(patches[i - 1], 1, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[Npatches_y + i]->setInnerBdry(patches[i + 1], 7, fringe,
    //                 overlap, unknowns, stage);     
    //         }
    //         patches[Npatches_y - 1]->setInnerBdry(patches[Npatches_y - 2], 2,
    //             fringe, overlap, unknowns, stage);   
    //         patches[Npatches_y - 1]->setInnerBdry(patches[2*Npatches_y - 1], 6,
    //             fringe, overlap, unknowns, stage); 
    //         patches[Npatches_y - 1]->setInnerBdry(patches[2*Npatches_y - 2], 3,
    //             fringe, overlap, unknowns, stage);     
    //         patches[2*Npatches_y - 1]->setInnerBdry(patches[2*Npatches_y - 2],
    //             2, fringe, overlap, unknowns, stage);   
    //         patches[2*Npatches_y - 1]->setInnerBdry(patches[Npatches_y - 1], 4,
    //             fringe, overlap, unknowns, stage); 
    //         patches[2*Npatches_y - 1]->setInnerBdry(patches[Npatches_y - 2], 1,
    //             fringe, overlap, unknowns, stage);   
    //     }
    // }
    // // else, if there are three or more patches in the x-direction
    // else
    // {
    //     // If there is only one patch in the y-direction
    //     if(Npatches_y == 1)
    //     {
    //         patches[0]->setInnerBdry(patches[1], 6, fringe, overlap, unknowns,
    //             stage);           
    //         for(int i = 1; i < Npatches_x - 1; i++)
    //         {
    //             patches[i]->setInnerBdry(patches[i + 1], 6, fringe, overlap,
    //                 unknowns, stage);
    //             patches[i]->setInnerBdry(patches[i - 1], 4, fringe, overlap,
    //                 unknowns, stage);
    //         }
    //         patches[Npatches_x - 1]->setInnerBdry(patches[Npatches_x - 2], 4,
    //             fringe, overlap, unknowns, stage);
    //     }
    //     // If there are only two patches in the y-direction
    //     else if(Npatches_y == 2)
    //     {
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);   
    //         patches[0]->setInnerBdry(patches[2], 6, fringe, overlap, unknowns,
    //             stage);
    //         patches[0]->setInnerBdry(patches[3], 9, fringe, overlap, unknowns,
    //             stage);    
    //         patches[1]->setInnerBdry(patches[3], 6, fringe, overlap, unknowns,
    //             stage);   
    //         patches[1]->setInnerBdry(patches[0], 4, fringe, overlap, unknowns,
    //             stage);  
    //         patches[1]->setInnerBdry(patches[2], 3, fringe, overlap,
    //             unknowns, stage);       
    //         for(int j = 1; j < Npatches_x - 1; j++)
    //         {
    //             patches[2*j]->setInnerBdry(patches[2*j - 2], 4, fringe, overlap,
    //                 unknowns, stage);   
    //             patches[2*j]->setInnerBdry(patches[2*j - 1], 7, fringe, overlap,
    //                 unknowns, stage); 
    //             patches[2*j]->setInnerBdry(patches[2*j + 1], 8, fringe, overlap,
    //                 unknowns, stage);  
    //             patches[2*j]->setInnerBdry(patches[2*j + 3], 9, fringe, overlap,
    //                 unknowns, stage); 
    //             patches[2*j]->setInnerBdry(patches[2*j + 2], 6, fringe, overlap,
    //                 unknowns, stage);    
    //             patches[2*j + 1]->setInnerBdry(patches[2*j - 2], 7, fringe, 
    //                 overlap, unknowns, stage);   
    //             patches[2*j + 1]->setInnerBdry(patches[2*j - 1], 4, fringe, 
    //                 overlap, unknowns, stage); 
    //             patches[2*j + 1]->setInnerBdry(patches[2*j], 2, fringe, overlap,
    //                 unknowns, stage);  
    //             patches[2*j + 1]->setInnerBdry(patches[2*j + 3], 6, fringe, 
    //                 overlap, unknowns, stage); 
    //             patches[2*j + 1]->setInnerBdry(patches[2*j + 2], 3, fringe, 
    //                 overlap, unknowns, stage);      
    //         }
    //         patches[2*Npatches_x - 2]->setInnerBdry(patches[2*Npatches_x - 1], 
    //             8, fringe, overlap, unknowns, stage);   
    //         patches[2*Npatches_x - 2]->setInnerBdry(patches[2*Npatches_x - 3], 
    //             7, fringe, overlap, unknowns, stage); 
    //         patches[2*Npatches_x - 2]->setInnerBdry(patches[2*Npatches_x - 4], 
    //             4, fringe, overlap, unknowns, stage);     
    //         patches[2*Npatches_x - 1]->setInnerBdry(patches[2*Npatches_x - 2],
    //             2, fringe, overlap, unknowns, stage);   
    //         patches[2*Npatches_x - 1]->setInnerBdry(patches[2*Npatches_x - 3],
    //             4, fringe, overlap, unknowns, stage); 
    //         patches[2*Npatches_x - 1]->setInnerBdry(patches[2*Npatches_x - 4], 
    //             7, fringe, overlap, unknowns, stage);   
    //     }

    //     // else, if there are three or more patches in the y-direction
    //     else
    //     {
    //         // deal with the corner patches first
    //         patches[0]->setInnerBdry(patches[1], 8, fringe, overlap, unknowns,
    //             stage);   
    //         patches[0]->setInnerBdry(patches[Npatches_y], 6, fringe, overlap,
    //             unknowns, stage);
    //         patches[0]->setInnerBdry(patches[Npatches_y + 1], 9, fringe,
    //             overlap, unknowns, stage); 
    //         patches[Npatches_y - 1]->setInnerBdry(patches[Npatches_y - 2], 2,
    //             fringe, overlap, unknowns, stage);   
    //         patches[Npatches_y - 1]->setInnerBdry(patches[2*Npatches_y - 2], 3,
    //             fringe, overlap, unknowns, stage);
    //         patches[Npatches_y - 1]->setInnerBdry(patches[2*Npatches_y - 1], 6,
    //             fringe, overlap, unknowns, stage); 
    //         patches[(Npatches_x - 1)*Npatches_y]->setInnerBdry(
    //             patches[(Npatches_x - 2)*Npatches_y], 4, fringe, overlap,
    //             unknowns, stage);   
    //         patches[(Npatches_x - 1)*Npatches_y]->setInnerBdry(
    //             patches[(Npatches_x - 2)*Npatches_y + 1], 7, fringe, overlap,
    //             unknowns, stage);  
    //         patches[(Npatches_x - 1)*Npatches_y]->setInnerBdry(
    //             patches[(Npatches_x - 1)*Npatches_y + 1], 8, fringe, overlap,
    //             unknowns, stage);  
    //         patches[Npatches_x*Npatches_y - 1]->setInnerBdry(
    //             patches[(Npatches_x - 1)*Npatches_y - 2], 1,
    //             fringe, overlap, unknowns, stage);   
    //         patches[Npatches_x*Npatches_y - 1]->setInnerBdry(
    //             patches[(Npatches_x - 1)*Npatches_y - 1], 4,
    //             fringe, overlap, unknowns, stage);  
    //         patches[Npatches_x*Npatches_y - 1]->setInnerBdry(
    //             patches[Npatches_x*Npatches_y - 2], 2,
    //             fringe, overlap, unknowns, stage);  

    //         // Patches on the left edge but not the corner
    //         for(int i = 1; i < Npatches_y - 1; i++)
    //         {
    //             patches[i]->setInnerBdry(patches[i - 1], 2, fringe, overlap,
    //                 unknowns, stage);   
    //             patches[i]->setInnerBdry(patches[Npatches_y + i - 1], 3, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[i]->setInnerBdry(patches[Npatches_y + i], 6, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[i]->setInnerBdry(patches[Npatches_y + i + 1], 9, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[i]->setInnerBdry(patches[i + 1], 8, fringe, overlap,
    //                 unknowns, stage);  
    //         }

    //         // Patches on the right edge but not the corner
    //         for(int i = 1; i < Npatches_y - 1; i++)
    //         {
    //             patches[(Npatches_x - 1)*Npatches_y + i]->setInnerBdry(
    //                 patches[(Npatches_x - 1)*Npatches_y + i - 1], 2, fringe,
    //                 overlap, unknowns, stage);   
    //             patches[(Npatches_x - 1)*Npatches_y + i]->setInnerBdry(
    //                 patches[(Npatches_x - 2)*Npatches_y + i - 1], 1, fringe,
    //                 overlap, unknowns, stage);  
    //             patches[(Npatches_x - 1)*Npatches_y + i]->setInnerBdry(
    //                 patches[(Npatches_x - 2)*Npatches_y + i], 4, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[(Npatches_x - 1)*Npatches_y + i + 1]->setInnerBdry(
    //                 patches[(Npatches_x - 2)*Npatches_y + i + 1], 7, fringe,
    //                 overlap, unknowns, stage); 
    //             patches[(Npatches_x - 1)*Npatches_y + i]->setInnerBdry(
    //                 patches[(Npatches_x - 1)*Npatches_y + i + 1], 8, fringe,
    //                 overlap, unknowns, stage);   
    //         }

    //         // Patches on the top edge but not corner
    //         for(int j = 1; j < Npatches_x - 1; j++)
    //         {
    //             patches[Npatches_y*j]->setInnerBdry(patches[Npatches_y*(j - 1)],
    //                 4, fringe, overlap, unknowns, stage);   
    //             patches[Npatches_y*j]->setInnerBdry(
    //                 patches[Npatches_y*(j - 1) + 1], 7, fringe, overlap,
    //                 unknowns, stage);   
    //             patches[Npatches_y*j]->setInnerBdry(patches[Npatches_y*j + 1],
    //                 8, fringe, overlap, unknowns, stage); 
    //             patches[Npatches_y*j]->setInnerBdry(
    //                 patches[Npatches_y*(j + 1) + 1], 9, fringe, overlap,
    //                 unknowns, stage); 
    //             patches[Npatches_y*j]->setInnerBdry(patches[Npatches_y*(j + 1)],
    //                 6, fringe, overlap, unknowns, stage);   
    //         }   

    //         // Patches on the bottom edge but not corner
    //         for(int j = 1; j < Npatches_x - 1; j++)
    //         {
    //             patches[Npatches_y*(j + 1) - 1]->setInnerBdry(
    //                 patches[Npatches_y*j - 2], 1, fringe, overlap, unknowns,
    //                 stage);   
    //             patches[Npatches_y*(j + 1) - 1]->setInnerBdry(
    //                 patches[Npatches_y*j - 1], 4, fringe, overlap, unknowns,
    //                 stage);   
    //             patches[Npatches_y*(j + 1) - 1]->setInnerBdry(
    //                 patches[Npatches_y*(j + 1) - 2], 2, fringe, overlap,
    //                 unknowns, stage); 
    //             patches[Npatches_y*(j + 1) - 1]->setInnerBdry(
    //                 patches[Npatches_y*(j + 2) - 2], 3, fringe, overlap,
    //                 unknowns, stage);  
    //             patches[Npatches_y*(j + 1) - 1]->setInnerBdry(
    //                 patches[Npatches_y*(j + 2) - 1], 6, fringe, overlap,
    //                 unknowns, stage);   
    //         }            

    //         // All the rest of the inner patches
    //         for(int j = 1; j < Npatches_x - 1; j++)
    //         {
    //             for(int i = 1; i < Npatches_y - 1; i++)
    //             {
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j - 1) + i - 1], 1, fringe, overlap,
    //                     unknowns, stage);   
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j - 1) + i], 4, fringe, overlap,
    //                     unknowns, stage);  
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j - 1) + i + 1], 7, fringe, overlap,
    //                     unknowns, stage); 
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*j + i - 1], 2, fringe, overlap,
    //                     unknowns, stage); 
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*j + i + 1], 8, fringe, overlap,
    //                     unknowns, stage); 
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j + 1) + i - 1], 3, fringe, overlap,
    //                     unknowns, stage);  
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j + 1) + i], 6, fringe, overlap,
    //                     unknowns, stage); 
    //                 patches[Npatches_y*j + i]->setInnerBdry(
    //                     patches[Npatches_y*(j + 1) + i + 1], 9, fringe, overlap,
    //                     unknowns, stage);  
    //             }
    //         }


    //     }

    // }

    for(int j = 0; j < Npatches_x; j++)
    {
        for(int i = 0; i < Npatches_y; i++)
        {
            // std::cout << std::endl;
            // std::cout << "Setting patch " << j*Npatches_y + i << std::endl;

            // Print_Patch2D(*patches[Npatches_y*j + i]);
            if (j > 0) 
            {
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j - 1) + i], 4, fringe, overlap,
                    unknowns, stage);  
            }
            if (i > 0) 
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*j + i - 1], 2, fringe, overlap,
                    unknowns, stage); 
            if (i < Npatches_y - 1)
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*j + i + 1], 8, fringe, overlap,
                    unknowns, stage); 
            if (j < Npatches_x - 1) 
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j + 1) + i], 6, fringe, overlap,
                    unknowns, stage); 
            if ((i > 0) && (j > 0))
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j - 1) + i - 1], 1, fringe, overlap,
                    unknowns, stage);  
            if ((i < Npatches_y - 1) && (j > 0))
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j - 1) + i + 1], 7, fringe, overlap,
                    unknowns, stage); 
            if ((i > 0) && (j < Npatches_x - 1))
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j + 1) + i - 1], 3, fringe, overlap,
                    unknowns, stage); 
            if ((i < Npatches_y - 1) && (j < Npatches_x - 1))
                patches[Npatches_y*j + i]->setInnerBdry(
                    patches[Npatches_y*(j + 1) + i + 1], 9, fringe, overlap,
                    unknowns, stage);

            // std::cout << std::endl;
            // patches[Npatches_y*j + i]->VectorFieldToNodes();
            // std::cout << "New patch " << j*Npatches_y + i << std::endl;
            // Print_Patch2D(*patches[Npatches_y*j + i]);
        }
    }

    
}