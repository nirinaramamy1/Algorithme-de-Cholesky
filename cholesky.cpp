/*
*   RAHENINTSOA
*   Maminirina
*   L3 MISA 2020-2021
*   nirinaramamy1@gmail.com
*   +261 34 23 944 07
*/

/*---------------------------------
Pour compiler le programme :

g++ maminirina.cpp -o maminirina
./maminirina

---------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

typedef vector<vector<float>> MatrixF;
typedef vector<float> VectorF;


void createMatrix(MatrixF& A, int dimension);
void createvector(VectorF& X, int dimension);
void getData(MatrixF& A, VectorF& B);
void showMatrix(MatrixF& matrix);
void showVector(VectorF& vector);
void showSystem(MatrixF A, VectorF B);
void solveSystemTriangularSuperior(MatrixF& A, VectorF& X, VectorF& B);
void solveSystemTriangularInferior(MatrixF& A, VectorF& X, VectorF& B);
int getDimensionFromData();
float sumInTriangularSuperior(MatrixF& A, VectorF& X, int i);
float sumInTriangularInferior(MatrixF& A, VectorF& X, int i);
float sumSquare(MatrixF& B, int i);
float sumProduct(MatrixF& B, int i, int j);
MatrixF transposeMatrix(MatrixF matrix);
MatrixF cholesky(MatrixF& matrix);


int main() {
    MatrixF A;
    VectorF X;
    VectorF b;
    
    MatrixF B;
    MatrixF TB;
    VectorF Y;

    getData(A, b);
    createMatrix(B, getDimensionFromData());
    createMatrix(TB, getDimensionFromData());
    createvector(Y, getDimensionFromData());
    createvector(X, getDimensionFromData());


    cout << "Système d'équation : " << endl << endl;
    showSystem(A, b);
    cout << endl << endl;
    cout << "La matrice A : " << endl << endl;
    showMatrix(A);

    cout << endl;
    cout << "Le vecteur b : " << endl << endl;
    showVector(b);

    cout << endl << endl;
    cout << "On cherche une matrice B telle que B est une matrice triangulaire dérivant de A " << endl;
    cout << "en utilisant l'Algorithme de Cholesky. Alors B est : " << endl << endl;
    B = cholesky(A);
    showMatrix(B);
 
    // Résolution du système B.Y = b
    cout << endl << endl;
    cout << "Puis on résoud deux systèmes triangulaires tels que : " << endl;
    cout << "B.Y = b" << endl << "On cherche Y = (X1,X2,X3,X4) " << endl << endl;
    showSystem(B, b);
    solveSystemTriangularInferior(B, Y, b);

    cout << endl << endl;
    cout << "Après résolution on obtient Y : " << endl << endl;
    showVector(Y);

    // Résolution du système TB.X = Y
    cout << endl << endl;
    cout << "et TB.X = Y, où TB est la matrice transposée de B.";
    TB = transposeMatrix(B);
    cout << endl << "On cherche X = (X1,X2,X3,X4) " << endl << endl;
    showSystem(TB, Y);
    solveSystemTriangularSuperior(TB, X, Y);

    cout << endl << endl;
    cout << "Après résolution on obtient X : " << endl << endl;
    showVector(X);
    cout << endl << endl;

    cout << "Alors la solution du système d'équation A.X = b : " << endl << endl;
    showSystem(A, b);
    cout << endl;

    cout << " est X : " << endl << endl;
    showVector(X);

    return EXIT_SUCCESS;
}


MatrixF cholesky(MatrixF& A) {
    MatrixF B;
    createMatrix(B, getDimensionFromData());
    
    for(int i = 0; i < B.size(); i++) {
        for(int j = 0; j < B.size(); j++) {
            if(j < i) {
                B[i][j] = (1.0/B[j][j])*(A[i][j] - sumProduct(B, i, j));
            } else if(j == i) {
                B[i][i] = sqrt(A[i][i] - sumSquare(B, i));
            } else {
                B[i][j] = 0;
            }
        }
    }
    return B;
}

MatrixF transposeMatrix(MatrixF matrix) {
    MatrixF T;
    createMatrix(T, getDimensionFromData());
    float temp{0};
    for(int i = 0; i < matrix.size(); i++) {
        for(int j = 0; j < matrix.size(); j++) {
            T[j][i] = matrix[i][j];
        }
    }
    return T;
}

float sumProduct(MatrixF& B, int i, int j) {
    float s{0};
    for(int k = 0; k < j; k++) {
        s += B[i][k]*B[j][k];
    }
    return s;
}

float sumSquare(MatrixF& B, int i) {
    float s{0};
    for(int k = 0; k < i; k++) {
        s += B[i][k]*B[i][k];
    }
    return s;
}

void createMatrix(MatrixF& A, int dimension) {
    VectorF X;
    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < dimension; j++) {
            X.push_back(0);
        }
        A.push_back(X);
        X.clear();
    }
}
void createvector(VectorF& X, int dimension) {
    for(int i = 0; i < dimension; i++) {
        X.push_back(0);
    }
}

void showSystem(MatrixF A, VectorF B) {
    for(int i = 0; i < A.size(); i++) {
        cout << '|';
        for(int j = 0; j < A.size(); j++) {
            if(j == A.size()-1) {
                cout << A[i][j] << '|';
            } else {
                cout << A[i][j] << setw(10);
            }
        }
        cout << setw(5) << "|X" << i+1 << '|' <<  setw(3) << '=' << setw(5) << B[i] << endl;
    }
}

void showMatrix(MatrixF& matrix) {
    for(int i = 0; i < matrix.size(); i++) {
        for(int j = 0; j < matrix.size(); j++) {
            cout << setw(10) << matrix[i][j];
        }
        cout << endl;
    }
}

void showVector(VectorF& vector) {
    for(int i = 0; i < vector.size(); i++) {
        cout << setw(10) << vector[i] << endl;
    }
}

int getDimensionFromData() {
    ifstream inFile{"data.txt"};
    string data{""};
    int N{0};
    if(inFile.is_open()) {
        inFile >> data;
        N = stoi(data);
    }
    inFile.close();
    return N;
}

void solveSystemTriangularSuperior(MatrixF& A, VectorF& X, VectorF& B) {
    for(int i = A.size()-1; i >= 0; i--) {
        X[i] = (1.0/(A[i][i]))*(B[i]-sumInTriangularSuperior(A,X,i));
    }
}
float sumInTriangularSuperior(MatrixF& A, VectorF& X, int i) {
    float s{0};
    for(int j = i+1; j < A.size(); j++) {
        s += A[i][j]*(X[j]);
    }
    return s;
}

void solveSystemTriangularInferior(MatrixF& A, VectorF& X, VectorF& B) {
    for(int i = 0; i < A.size(); i++) {
        X[i] = (1.0/(A[i][i]))*(B[i]-sumInTriangularInferior(A,X,i));
    }
}

float sumInTriangularInferior(MatrixF& A, VectorF& X, int i) {
    float s{0};
    for(int j = A.size()-1; j >= 0; j--) {
        s += A[i][j]*(X[j]);
    }
    return s;
}

void getData(MatrixF& A, VectorF& B) {
    //Récupérer les données
    ifstream inFile{"data.txt"};
    string data{""};
    VectorF X;
    int N{0};

    if(inFile.is_open()) {
        //Récupérer la dimension
        inFile >> data;
        N = stoi(data);

        //Récupérer la matrice
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                inFile >> data;
                X.push_back(stof(data));
            }
            A.push_back(X);
            X.clear();
        }

        //Récupérer le vecteur
        for(int i = 0; i < N; i++) {
            inFile >> data;
            B.push_back(stof(data));
        }
    }
    inFile.close();
}