#include <stdio.h>
#include "skinGhostAvg_mod_decl.h"



void print_arr(double* arr, int length){
    for(int i=0; i<length; i++){
      //printf("%d\n", i);
      printf("%f,",arr[i]);
      //printf("%p\n", &arr[i]);
    }
    printf("\n");
}

void fill_ones(double* arr, int length){
    for(int i=0; i<length; i++){
      arr[i] = 1.0;
    }
}

void fill_ascending(double* arr, int length){
    for(int i=0; i<length; i++){
      arr[i] = i+1.0;
    }
}

int main(void){
    printf("%s:%s:%d \n", __FILE__, __FUNCTION__, __LINE__);
    double phi[8];
    double phiGhost[8];
    int l_phi = sizeof(phi)/sizeof(phi[0]);
    //printf("lphi = %d\n",l_phi );
    //print_arr(phi,l_phi);
    fill_ones(phi,l_phi);
    printf("phi interior, original modal coeffs: ");
    print_arr(phi,l_phi);

    fill_ascending(phiGhost,l_phi);
    printf("phi ghost , original modal coeffs: ");
    print_arr(phiGhost,l_phi);

    skin_ghost_avg_lower_3x_Ser_p1(phi, phiGhost);
    //printf("test function\n");
    //test_func(phi, length);
    printf("Final modal coeffs (lower): ");
    print_arr(phi,l_phi);

    printf("\n Now test upper boundary\n");
    fill_ones(phi,l_phi);
    printf("phi interior, original modal coeffs: ");
    print_arr(phi,l_phi);

    fill_ascending(phiGhost,l_phi);
    printf("phi ghost , original modal coeffs: ");
    print_arr(phiGhost,l_phi);

    skin_ghost_avg_upper_3x_Ser_p1(phi, phiGhost);
    //printf("test function\n");
    //test_func(phi, length);
    printf("Final modal coeffs (upper): ");
    print_arr(phi,l_phi);

    return 0;
}

