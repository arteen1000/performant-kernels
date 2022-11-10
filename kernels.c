#include "kernels.h"

// Ni,Nj -- Dimensions of the In/Out matricies
void compute_transpose(int Ni, int Nj,
                    const char In[Ni][Nj], char Out[Nj][Ni]) {
  for(int j = 0; j < Nj; ++j) {
    for(int i = 0; i < Ni; ++i) {
      move_one_value(i,j,j,i,Ni,Nj,In,Out);
    }
  }
}


// Ni,Nj,Nk -- Dimensions of the output matrix
// S -- width/length/height of the stencil
void compute_stencil(int Ni, int Nj, int Nk, int S, 
            const float In[Ni+S][Nj+S][Nk+S], float Out[Ni][Nj][Nk], 
            const float Stencil[S][S][S]) {
  for(int k = 0; k < Nk; ++k) {
    for(int j = 0; j < Nj; ++j) { 
      for(int i = 0; i < Ni; ++i) { 
        set_to_zero(i,j,k,Ni,Nj,Nk,Out);
        for(int z = 0; z < S; ++z) {
          for(int y = 0; y < S; ++y) {
            for(int x = 0; x < S; ++x) {
              macc_element(&In[i+x][j+y][k+z],&Out[i][j][k],&Stencil[x][y][z]);
        } } }
  } } }
}




