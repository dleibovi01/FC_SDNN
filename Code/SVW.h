// Smooth viscosity window object

#ifndef SVW_H
#define SVW_H

#include <vector>
#include "mkl_spblas.h"
#include "printing.h"
#include "SpMatrix_csr.h"
#include <cmath>


const double tolerance = std::pow(10.0, -15.0);
constexpr double window_cutoff = 9.0; // 9.0
const double pi = std::acos(-1);


class SVW {

std::vector<int> patch_ids;
std::vector<SpMatrix_csr*> patch_svws;

public:

    ~SVW()
    {
        // std::cout << "In SVW destructor" << std::endl;
        while(patch_svws.size() > 0)
        {
            delete patch_svws.back();
            patch_svws.pop_back();
        }
    }

    SVW() {};
    SVW(const SVW &svw)
    {
        std::cout << "copy constructor " << std::endl;
        patch_ids = svw.patch_ids;
    }

    SVW & operator= (const SVW &svw) 
    {
        return *this;
    }

    const std::vector<int> & getPatchIds() const {return patch_ids;}
    const std::vector<SpMatrix_csr*> & getPatchSVWS() const {return patch_svws;}

    template<typename Patch>
    void compute_svw_data(Patch * patch_origin, Patch * patch_target, 
        int target_id)
    {
        int rows = patch_origin->getNnodes();
        int cols = patch_target->getNnodes();
        auto nodes_origin = patch_origin->getNodes();
        auto nodes_target = patch_target->getNodes();
        std::vector<int> inner_nodes = patch_target->getInnerNodes();
        std::vector<double> values;
        std::vector<int> rows_start;
        std::vector<int> rows_end;
        std::vector<int> col_indx;
        double value;
        double h = patch_target->getH();
        double cutoff_h = window_cutoff*h;
        auto node_origin = nodes_origin[0];
        auto node_target = nodes_target[0];
        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ONE;
        sparse_status_t status;
        bool new_row;
        int last_index = 1;
        rows_start.push_back(last_index);

        for(int i = 0; i < rows; i++)
        {
            new_row = false;
            node_origin = nodes_origin[i];
            for(int j = 0; j < inner_nodes.size(); j++)
            {
                node_target = nodes_target[inner_nodes[j]];
                value = visc_window(*node_origin, *node_target, cutoff_h); // jkl
                if(value > tolerance)
                {
                    values.push_back(value);
                    col_indx.push_back(inner_nodes[j] + 1);
                    last_index++;
                    new_row = true;
                }
            }
            rows_end.push_back(last_index);
            if(i < rows - 1)
            {
                rows_start.push_back(last_index);
            }             
        }   
        if(values.size() > tolerance)
        {
            patch_ids.push_back(target_id);
            patch_svws.push_back(new SpMatrix_csr{indexing, rows, cols,
                rows_start.data(), rows_end.data(), col_indx.data(),
                values.data()});
            patch_svws.back()->setMVHint();
            patch_svws.back()->optimize();
            
        }
    }


    void setScaledSVWS()
    {
        int N = patch_svws.size();
        int Nnodes = patch_svws[0]->getRows();
        std::vector<double> zeros(Nnodes, 0.0);
        double sum_total[Nnodes];
        double sum_temp[Nnodes];
        for(int i = 0; i < Nnodes; i++)
        {
            sum_total[i] = 0.0;
            sum_temp[i] = 0.0;
        }
        for(int i = 0; i < N; i++)
        {
            patch_svws[i]->rowSum(sum_temp);
            vdAdd(Nnodes, sum_temp, sum_total, sum_total);
        }
        for(int i = 0; i < Nnodes; i++)
        {
            sum_total[i] = 1.0 / sum_total[i];
        }
        for(int i = 0; i < N; i++)
        {
            patch_svws[i]->dotMV(sum_total);
        }
    }

private:

    template<typename Node>
    double visc_window(const Node & node1, const Node & node2, double cutoff)
    {
        double distance = node1.getDist(node2);
        // const double pi = std::acos(-1);
        double H;
        if (distance < cutoff)
        {
            H = std::cos(pi*distance / 2.0 / cutoff);
            return H*H;
        }
        else
        {
            return 0.0;
        }
    }   
    


};

class SVW_mesh{

std::vector<SVW*> svws;

public:

    SVW_mesh(){};
    
    template<typename Mesh>
    SVW_mesh(const Mesh &mesh)
    {
        auto patches = mesh.getPatches();
        int npatches = patches.size();
        for(int i = 0; i < npatches; i++)
        {
            svws.push_back(new SVW);
        }
        setSVWs(mesh);
    }

    ~SVW_mesh()
    {
        while(svws.size() > 0)
        {
            delete svws.back();
            svws.pop_back();
        }        
    }

    const std::vector<SVW*> & getSVWs() const {return svws;}

private: 

    template<typename Mesh>
    void setSVWs(const Mesh &mesh)
    {
        auto patches = mesh.getPatches();
        auto patch_origin = patches[0];
        auto patch_target = patches[0];
        int npatches = patches.size();
        // Compute the unscaled windows
        for(int i = 0; i < npatches; i++)
        {
            for(int j = 0; j < npatches; j++)
            {
                svws[i]->compute_svw_data(patches[i], patches[j], j);
            }
        }
        // Scale the windows
        for(int i = 0; i < npatches; i++)
        {
            svws[i]->setScaledSVWS();
        }
    }    
};


#endif