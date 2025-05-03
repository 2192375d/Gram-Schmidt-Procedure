#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dot_product(int n, double v[n], double w[n]){
    double sum = 0.0;
    for(int i = 0; i < n; i++){
        sum += v[i] * w[i];
    }

    return sum;
}

double norm(int n, double v[n]){
    return sqrt(dot_product(n, v, v));
}

void vector_addition(int n, double v[n], double w[n], double result[n]){
    for(int i = 0; i < n; i++){
        result[i] = v[i] + w[i];
    }
}

void scalar_multiplication(int n, double v[n], double c, double result[n]){
    for(int i = 0; i < n; i++){
        result[i] = c * v[i];
    }
}

void assign_equal(int n, double source[n], double result[n]){
    for(int i = 0; i < n; i++){
        result[i] = source[i];
    }
}

void print_vector(int n, double v[n]){
    printf("(%f", v[0]);
    for(int i = 1; i < n; i++){
        printf(", %f", v[i]);
    }
    printf(")\n");
}

void Gram_Schmidt_Procedure(int m, int n, double LI_list[m][n], double orthonormal_list[m][n]){
    
    //Deal with e_1
    scalar_multiplication(n, LI_list[0], 1/norm(n, LI_list[0]), orthonormal_list[0]);

    //Deal with e_2,...,e_m
    for(int i = 1; i < m; i++){
        double v[n]; // v[i] = v_i  - <v_i, e_1>e_1-...-<v_i, e_{i-1}>e_{i-1}
        assign_equal(n, LI_list[i], v);

        //Create v
        double temp[n];

        for(int j = 0; j < i; j++){
            
            scalar_multiplication(n, orthonormal_list[j], -1, temp);
            scalar_multiplication(n, temp, dot_product(n, LI_list[i], orthonormal_list[j]), temp);
            vector_addition(n, v, temp, v);
        }

        scalar_multiplication(n, v, 1/norm(n, v), orthonormal_list[i]);
    }
}

void test1(){
    double LI_list[3][3] = {
        {1.0, 0.0, 1.0},
        {0.0, 2.0, 1.0},
        {3.0, 4.0, 2.0}
    };

    double orthonormal_list[3][3];
    Gram_Schmidt_Procedure(3,  3, LI_list, orthonormal_list);

    for(int i = 0; i < 3; i++){
        print_vector(3, orthonormal_list[i]);
    }
}

// void test2(){
//     double LI_list[2][3] = {
//         {1.0, 0.0, 1.0},
//         {0.0, 2.0, 1.0}
//     };

//     double orthonormal_list[2][3];
//     Gram_Schmidt_Procedure(2,  3, LI_list, orthonormal_list);

//     for(int i = 0; i < 2; i++){
//         print_vector(3, orthonormal_list[i]);
//     }
// }

int main(){
    test1();

    return 0;
}
