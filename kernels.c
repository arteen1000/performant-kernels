#include "kernels.h"

// *******
// BEFORE
//********

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



// *******
// AFTER
// *******

// #include "kernels.h"

// void compute_transpose(int Ni, int Nj, const char In[restrict 10000][10000], char Out[restrict 10000][10000])
//  {
//   const int B = 125;  // 100 50 125 200 250
//   int i, j;
//   for ( j = 0 ; j < Nj ; j+=B)     // go through all rows of In
//   {
//     for (i = 0 ; i < Ni ; i+=B)   // go through all cols of In
//     {
//       for (int j1 = j; j1 < j + B; ++j1)
//       {
//         for (int i1 = i; i1 < i + B ; ++i1)
//         {
//           Out[j1][i1] = In[i1][j1];
//         }
//       }
//     }
//   }
// }

// void compute_transpose(int Ni, int Nj, const char In[restrict 10000][10000], char Out[restrict 10000][10000])
//  {
//   const int B1 = 250;  // 100 50 125 200 250 (500 is ideal!)    // optimizes blocks for L2 (256K)
//   const int B2 = 125; // (optimizes blocks for L1)

//   int i, j;
//   for (j = 0 ; j < Ni ; j+=B1)
//   {
//     for (i = 0 ; i < Nj ; i+=B1)   // go through all cols of In
//     {
//       for (int j1 = j; j1 < j + B1; j1+=B2)
//       {
//         for (int i1 = i; i1 < i + B1; i1+=B2)
//         {
//           for (int j2 = j1; j2 < j1 + B2; j2++)
//           {
//             for(int i2 = i1; i2 < i1 + B2; i2++)
//             {
//                   Out[j2][i2] = In[i2][j2];
//             }
//           }
//         }
//       }
//     }
//   }
// }

void compute_transpose(int Ni, int Nj, const char In[restrict 10000][10000], char Out[restrict 10000][10000])
 {
  const int B1 = 500;  // 100 50 125 200 250 (500 is ideal!)    // optimizes blocks for L2 (256K)
  const int B2 = 125; // (optimizes blocks for L1)


  for (int j = 0 ; j < Ni ; j+=B1)
  {
    for (int i = 0 ; i < Nj ; i+=B1)   // go through all cols of In
    {
      for (int j1 = j; j1 < j + B1; j1+=B2)
      {
        for (int i1 = i; i1 < i + B1; i1+=B2)
        {
          for (int j2 = j1; j2 < j1 + B2; j2++)
          {
            for(int i2 = i1; i2 < i1 + B2; i2++)
            {
                  Out[i2][j2] = In[j2][i2];
            }
          }
        }
      }
    }
  }
}

// Ni,Nj,Nk -- Dimensions of the output matrix
// S -- width/length/height of the stencil
// void compute_stencil(int Ni, int Nj, int Nk, int S, const float In[Ni+S][Nj+S][Nk+S], float Out[Ni][Nj][Nk], const float Stencil[S][S][S]) 
// {
//   for(int k = 0; k < Nk; ++k) {
//     for(int j = 0; j < Nj; ++j) { 
//       for(int i = 0; i < Ni; ++i) { 
//         set_to_zero(i,j,k,Ni,Nj,Nk,Out);
//         array[i][j][k]=0;
//         for(int z = 0; z < S; ++z) {
//           for(int y = 0; y < S; ++y) {
//             for(int x = 0; x < S; ++x) {
//               macc_element(&In[i+x][j+y][k+z],&Out[i][j][k],&Stencil[x][y][z]);
//               *Out += (*In) * (*Stencil);
//         } } }
//   } } }
// }

// void compute_stencil(int Ni, int Nj, int Nk, int S, const float In[Ni+S][Nj+S][Nk+S], float Out[Ni][Nj][Nk], const float Stencil[ S][S][S]) 
// {
//   float sum = 0;
//   for(int i = 0; i < Ni; ++i) {
//     for(int j = 0; j < Nj; ++j) { 
//       for(int k = 0; k < Nk; ++k) { 
//         // Out[i][j][k]=0;
//         sum = 0;
//         for(int x = 0; x < S; ++x) {

//         for(int y = 0; y < S; ++y) {
          
//         for(int z = 0; z < S; ++z) {
//               sum += In[i+x][j+y][k+z] * Stencil[x][y][z];
//         } 
//         } 
//         }
//         Out[i][j][k] = sum;
//   } } }
// }


void compute_stencil(int Ni, int Nj, int Nk, int S, const float In[restrict Ni+S][Nj+S][Nk+S], float Out[restrict Ni][Nj][Nk], const float Stencil[restrict S][S][S]) 
{
  if (S == 8)
  {
    compute_stencil1(In, Out, Stencil);
  }
  else if (S == 20)
  {
    compute_stencil2(In, Out, Stencil);
  }
  else
  {
    compute_stencil3(In, Out, Stencil);
}
}

// void compute_stencil1(const float In[restrict 136][136][136], float Out[restrict 128][128][128], const float Stencil[restrict 8][8][8]) 
// {
//   const int N = 128;
//   const int S = 8;

//   for(int i = 0; i < N; ++i) {
//     for(int j = 0; j < N; ++j) { 
//       for(int k = 0; k < N; ++k) { 
//         Out[i][j][k] = 0; } } }
    
//   for (int i = 0 ; i < N ; i++){
//     for (int x = 0 ; x < S ; x++){
//       for (int j = 0; j < N ; j++){
//         for (int y = 0 ; y < S ; y++){
//           for (int k = 0 ; k < N ; k++){
//             for (int z = 0 ; z < S ; z++){
//               Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];
//             }
//           }
//         }
//       }
//     }
//   }
// }

void compute_stencil1(const float In[restrict 136][136][136], float Out[restrict 128][128][128], const float Stencil[restrict 8][8][8]) 
{
  const int N = 128;
  const int S = 8;

  const int B1 = 16;
  const int B2 = 4;

    for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) { 
      for(int k = 0; k < N; ++k) { 
        Out[i][j][k] = 0; } } }
    // Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];

  for (int i = 0 ; i < N ; i+=B1){
    for (int x = 0 ; x < S ; x+=B2){
      for (int j = 0; j < N ; j+=B1){
        for (int y = 0 ; y < S ; y+=B2){
          for (int k = 0 ; k < N ; k+=B1){
            for (int z = 0 ; z < S ; z+=B2){
              for (int i1 = i ; i1 < i + B1 ; i1++){
                for (int x1 = x; x1 < x + B2 ; x1++){
                  for (int j1 = j ; j1 < j + B1 ; j1++){
                    for (int y1 = y ; y1 < y + B2 ; y1++){
                      for (int k1 = k; k1 < k + B1 ; k1++){
                        for (int z1 = z ; z1 < z + B2 ; z1++){
                          Out[i1][j1][k1] += In[i1+x1][j1+y1][k1+z1] * Stencil[x1][y1][z1];
                          // Out[i1][j1][k1] += In[i1+x1][j1+y1][k1+z1] * Stencil[x1][y1][z1];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

// void compute_stencil1(const float In[restrict 136][136][136], float Out[restrict 128][128][128], const float Stencil[restrict 8][8][8]) 
// {
//   const int N = 128;
//   const int S = 8;
//   const int B1 = 8;
//   const int BB = 16;

//   const int B2 = 4;

//     for(int i = 0; i < N; ++i) {
//     for(int j = 0; j < N; ++j) { 
//       for(int k = 0; k < N; ++k) { 
//         Out[i][j][k] = 0; } } }
//     // Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];

// for (int ii = 0; ii < N ; ii+=BB){
//   for (int jj = 0; jj < N ; jj += BB){
//   for (int kk = 0 ; kk < N ; kk += BB){
//   for (int i = ii ; i < ii+BB ; i+=B1){
//     for (int x = 0 ; x < S ; x+=B2 ){
//       for (int j = jj; j < jj+BB ; j+=B1){
//         for (int y = 0 ; y < S ; y+=B2){
//           for (int k = kk ; k < kk + BB ; k+=B1){
//             for (int z = 0 ; z < S ; z+=B2){
//               for (int i1 = i ; i1 < i + B1 ; i1++){
//                 for (int x1 = x; x1 < x + B2 ; x1++){
//                   for (int j1 = j ; j1 < j + B1 ; j1++){
//                     for (int y1 = y ; y1 < y + B2 ; y1++){
//                       for (int k1 = k; k1 < k + B1 ; k1++){
//                         for (int z1 = z ; z1 < z + B2 ; z1++){
//                           Out[i1][j1][k1] += In[i1+x1][j1+y1][k1+z1] * Stencil[x1][y1][z1];
//                         }
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }}}}
// }



// void compute_stencil2(const float In[restrict 84][84][84], float Out[restrict 64][64][64], const float Stencil[restrict 20][20][20]) 
// {
//   float sum = 0;
//   const int N = 64;
//   const int S = 20;

//     for(int i = 0; i < N; ++i) {
//     for(int j = 0; j < N; ++j) { 
//       for(int k = 0; k < N; ++k) { 
//         Out[i][j][k] = 0; } } }
    
//   for (int i = 0 ; i < N ; i++){
//     for (int x = 0 ; x < S ; x++){
//       for (int j = 0; j < N ; j++){
//         for (int y = 0 ; y < S ; y++){
//           for (int k = 0 ; k < N ; k++){
//             for (int z = 0 ; z < S ; z++){
//               Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];
//             }
//           }
//         }
//       }
//     }
//   }
// }

// void compute_stencil2(const float In[restrict 84][84][84], float Out[restrict 64][64][64], const float Stencil[restrict 20][20][20]) 
// {
//   const int N = 64;
//   const int S = 20;
//   const int B1 = 8;
//   const int B2 = 4;

//     for(int i = 0; i < N; ++i) {
//     for(int j = 0; j < N; ++j) { 
//       for(int k = 0; k < N; ++k) { 
//         Out[i][j][k] = 0; } } }
//     // Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];

//   for (int i = 0 ; i < N ; i+=B1){
//     for (int x = 0 ; x < S ; x+=B2 ){
//       for (int j = 0; j < N ; j+=B1){
//         for (int y = 0 ; y < S ; y+=B2){
//           for (int k = 0 ; k < N ; k+=B1){
//             for (int z = 0 ; z < S ; z+=B2){
//               for (int i1 = i ; i1 < i + B1 ; i1++){
//                 for (int x1 = x; x1 < x + B2 ; x1++){
//                   for (int j1 = j ; j1 < j + B1 ; j1++){
//                     for (int y1 = y ; y1 < y + B2 ; y1++){
//                       for (int k1 = k; k1 < k + B1 ; k1++){
//                         for (int z1 = z ; z1 < z + B2 ; z1++){
//                           Out[i1][j1][k1] += In[i1+x1][j1+y1][k1+z1] * Stencil[x1][y1][z1];
//                         }
//                       }
//                     }
//                   }
//                 }
//               }
//             }
//           }
//         }
//       }
//     }
//   }
// }

void compute_stencil2(const float In[restrict 84][84][84], float Out[restrict 64][64][64], const float Stencil[restrict 20][20][20]) 
{
  const int N = 64;
  const int S = 20;
  const int B1 = 8;
  const int BB = 16;

  const int B2 = 4;

    for(int i = 0; i < N; ++i) {
    for(int j = 0; j < N; ++j) { 
      for(int k = 0; k < N; ++k) { 
        Out[i][j][k] = 0; } } }
    // Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];

for (int ii = 0; ii < N ; ii+=BB){
  for (int jj = 0; jj < N ; jj += BB){
  for (int kk = 0 ; kk < N ; kk += BB){
  for (int i = ii ; i < ii+BB ; i+=B1){
    for (int x = 0 ; x < S ; x+=B2 ){
      for (int j = jj; j < jj+BB ; j+=B1){
        for (int y = 0 ; y < S ; y+=B2){
          for (int k = kk ; k < kk + BB ; k+=B1){
            for (int z = 0 ; z < S ; z+=B2){
              for (int i1 = i ; i1 < i + B1 ; i1++){
                for (int x1 = x; x1 < x + B2 ; x1++){
                  for (int j1 = j ; j1 < j + B1 ; j1++){
                    for (int y1 = y ; y1 < y + B2 ; y1++){
                      for (int k1 = k; k1 < k + B1 ; k1++){
                        for (int z1 = z ; z1 < z + B2 ; z1++){
                          Out[i1][j1][k1] += In[i1+x1][j1+y1][k1+z1] * Stencil[x1][y1][z1];
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }}}}
}

// void compute_stencil3(const float In[restrict 6][6][1048578], float Out[restrict 4][4][1048576], const float Stencil[restrict 2][2][2]) 
// {
//   const int Nij = 4;
//   const int Nk = 1048576;
//   const int S = 2;
//   float sum = 0;
//   for(int i = 0; i < Nij; ++i) {
//     for(int j = 0; j < Nij; ++j){ 
//        for(int k = 0; k < Nk; ++k){ 
//         sum=0;
//         for(int x = 0; x < S; ++x) {
//         for(int y = 0; y < S; ++y) {
//         for(int z = 0; z < S; ++z) {
//               sum += In[i+x][j+y][k+z] * Stencil[x][y][z];
//         } 
//         } 
//         }
//         Out[i][j][k]=sum;
//   } } }

// }

// void compute_stencil3(const float In[restrict 6][6][1048578], float Out[restrict 4][4][1048576], const float Stencil[restrict 2][2][2]) 
// {
//   const int Nij = 4;
//   const int Nk = 1048576;
//   const int S = 2;
//   const int B = 2000;
//   float sum = 0;

// for(int k = 0; k < 1048000; k+=B){ 
//   for(int i = 0; i < Nij; ++i) {
//     for(int j = 0; j < Nij; ++j){ 
//        for(int k1 = k; k1 < k+B ; k1++){ 
          
//         sum=0;
//         // Out[i][j][k1] = 0;
//         for(int x = 0; x < S; ++x) {

//         for(int y = 0; y < S; ++y) {
          
//         for(int z = 0; z < S; ++z) {
//               sum += In[i+x][j+y][k1+z] * Stencil[x][y][z];
//         }
//         } 
//         }
//         Out[i][j][k1] = sum;
          
//   } } }
// }


//   for(int i = 0; i < Nij; ++i) {
//     for(int j = 0; j < Nij; ++j){ 
//        for(int k = 1048000; k < Nk; k++){ 
          
//         sum=0;
//         // Out[i][j][k1] = 0;
//         for(int x = 0; x < S; ++x) {

//         for(int y = 0; y < S; ++y) {
          
//         for(int z = 0; z < S; ++z) {
//               sum += In[i+x][j+y][k+z] * Stencil[x][y][z];
//         } 
//         } 
//         }
//         Out[i][j][k] = sum;
       
//   } } }


// }

void compute_stencil3(const float In[restrict 6][6][1048578], float Out[restrict 4][4][1048576], const float Stencil[restrict 2][2][2]) 
{
  const int Nij = 4;
  const int Nk = 1048576;
  const int S = 2;
  const int B1 = 8000;
  const int B2 = 2000;
  float sum = 0;

for(int k = 0; k < 1048000; k+=B1){ 
  for(int i = 0; i < Nij; ++i){
   for (int k1 = k ; k1 < k + B1; k1+=B2){
    for(int j = 0; j < Nij; ++j){ 
       for(int k2 = k1; k2 < k1+B2 ; k2++){ 
          
        sum = 0;
        // Out[i][j][k1] = 0;
        for(int x = 0; x < S; ++x) {

        for(int y = 0; y < S; ++y) {
          
        for(int z = 0; z < S; ++z) {
              sum += In[i+x][j+y][k2+z] * Stencil[x][y][z];
        }
        } 
        }
        Out[i][j][k2] = sum;
          
  } } }}
}


  for(int i = 0; i < Nij; ++i) {
    for(int j = 0; j < Nij; ++j){ 
       for(int k = 1048000; k < Nk; k++){ 
          
        sum=0;
        // Out[i][j][k1] = 0;
        for(int x = 0; x < S; ++x) {

        for(int y = 0; y < S; ++y) {
          
        for(int z = 0; z < S; ++z) {
              sum += In[i+x][j+y][k+z] * Stencil[x][y][z];
        } 
        } 
        }
        Out[i][j][k] = sum;
       
  } } }


}



// void compute_stencil3(const float In[restrict 6][6][1048578], float Out[restrict 4][4][1048576], const float Stencil[restrict 2][2][2]) 
// {
//   const int Nij = 4;
//   const int Nk = 1048576;
//   const int S = 2;
//   const int B = 500;
//   float sum = 0;

//   for(int i = 0; i < Nij; ++i) {
//     for(int j = 0; j < Nij; ++j){ 
//        for(int k = 0; k < 1048000; k+=B){ 
//           for(int k1 = k; k1 < k+B; k1++){
//         sum=0;
//         // Out[i][j][k1] = 0;
//         for(int x = 0; x < S; ++x) {

//         for(int y = 0; y < S; ++y) {
          
//         for(int z = 0; z < S; ++z) {
//               sum += In[i+x][j+y][k1+z] * Stencil[x][y][z];
//         } 
//         } 
//         }
//         Out[i][j][k1] = sum;
//        }
//   } } }

//   for(int i = 0; i < Nij; ++i) {
//     for(int j = 0; j < Nij; ++j){ 
//        for(int k = 1048000; k < Nk; k++){ 
//         Out[i][j][k]=0;
//         for(int x = 0; x < S; ++x) {

//         for(int y = 0; y < S; ++y) {
          
//         for(int z = 0; z < S; ++z) {
//               Out[i][j][k] += In[i+x][j+y][k+z] * Stencil[x][y][z];
//         } 
//         } 
//         }
//   } } }

// }


/*
32K L1 cache
Stencil: 
1: 8 x 8 (with 4 byte floats) = 64 * 4 = 256 bytes for the stencil, row major ordering S[0][0] - S[0][7] .. S[1][0] - S[1][7] .. S[7][0] - S[7][7]
2: 20 x 20 (with 4 byte floats) = 400 * 4 = 1600 bytes for the stencil

Cache lines are 64 bytes
*/
// look for reuse and EXPLOIT IT (tiling)






