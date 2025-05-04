#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double dot_product(int n, double v[n], double w[n]){
    /*
     * The function returns the dot product of two vectors
     *
     * int n: Size of both vectors
     * double v[n]: first vector
     * double w[n]: second vector
    */

    double sum = 0.0;
    for(int i = 0; i < n; i++){
        sum += v[i] * w[i];
    }

    return sum;
}

double norm(int n, double v[n]){
    /*
     * The function returns the norm(magnitude) of a vector
     *
     * int n: Size of the vectors
     * double v[n]: the vector to find the norm for
    */

    return sqrt(dot_product(n, v, v));
}

void vector_addition(int n, double v[n], double w[n], double result[n]){
    /*
     * The function takes two vectors as input, and pass their sum to the 
     * result vector
     *
     * int n: Size of the vectors
     * double v[n]: the first vector to add
     * double w[n]: the second vector to add
     * double result[n]: the vector to pass the sum to
     * 
     * ex:  a[2] = {1, 1}; b[2] = {2, 2}; result[2];
     *      vector_addition(2, a, b, result);
     *      printf("(%f, %f)", result[0], result[1]);
     * 
     * The output will be: (3, 3)
    */

    for(int i = 0; i < n; i++){
        result[i] = v[i] + w[i];
    }
}

void scalar_multiplication(int n, double v[n], double c, double result[n]){
    /*
     * The function takes a vectors and a scalar as input, and pass the result 
     * of scalar multiplication to the result vector
     *
     * int n: Size of the vectors
     * double v[n]: the first vector to add
     * double c: the scalar to multiply to
     * double result[n]: the vector to pass the result of scalar multiplication 
     * to
     * 
     * ex:  a[2] = {1, 2}; c = 3; result[2];
     *      scalar_multiplication(2, a, c, result);
     *      printf("(%f, %f)", result[0], result[1]);
     * 
     * The output will be: (3, 6)
    */

    for(int i = 0; i < n; i++){
        result[i] = c * v[i];
    }
}

void assign_equal(int n, double source[n], double result[n]){
    /*
     * The function copies one array into the other array
     * 
     * int n: size of the arrays
     * double source[n]: the array to copy from
     * double result[n]: the array to copy to
    */

    for(int i = 0; i < n; i++){
        result[i] = source[i];
    }
}

void print_vector(int n, double v[n]){
    /*
     * The function prints individual vectors in form of (a_1, a_2, ..., a_n)
     * 
     * int n: size of the vector
     * double v[n]: the vector to print
     * 
     * ex:  a[4] = {1, 2, 5, 4};
     *      print_vector(4, a);
     * 
     * The output will be: (1, 2, 5, 4)
     * 
     * (note that it skips the line after the end)
    */

    printf("(%f", v[0]);
    for(int i = 1; i < n; i++){
        printf(", %f", v[i]);
    }
    printf(")\n");
}

void Gram_Schmidt_Procedure(int m, int n, double LI_list[m][n], double orthonormal_list[m][n]){
    /*
     * The function takes an input linearly independent list of vectors, and applies
     * Gram-Schmidt Procedure to get an orthonormal list with the same span as the
     * input linearly independent list.
     * (What is orthonormal, and how exactly does that work? Check what is in the 
     * readme.md file)
     * 
     * int m: the number of vectors in the linearly independent list
     * int n: the number of components in each vectors in the linearly independent
     * list (ex: (8, 5, 1) has 3 components)
     * double LI_list[m][n]: The input linearly independent list (assumed to be
     * linearly independent as a precondition/resonable assumption)
     * orthonormal_list[m][n]: the orthonormal list obtained from this function, 
     * which will have the same span as the linearly independent list.
    */
    
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
