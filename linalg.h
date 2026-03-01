#pragma once

#include <vector>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <initializer_list>

//NOTE: Would be good to add checks for invert and dimensionality

// Vector class

class Vector
{
private:
    
    std::vector<double> data;

public:

    // Constructors
    explicit Vector(size_t n);
    Vector(std::initializer_list<double> list);

    //Get size method
    size_t size() const;

    //Get elements methods
    double& operator[](size_t i);
    const double& operator[](size_t i) const;

    //Utility method
    void wipe();

    // Arithmetic operators
    Vector operator+(const Vector& other) const;
    Vector operator-(const Vector& other) const;
    Vector operator*(double scalar) const; // for dividing entire vec by constant just use inverse

    // Dot and cross product
    double dot(const Vector& other) const;
    Vector cross(const Vector& other) const;

    // Element wise operations
    Vector ewiseDiv(const Vector& other) const; // not dividing all elements by constant but element-wise division (not standard but practical)
    Vector ewiseMult(const Vector& other) const; // element-wise multiplication (gets rid of remainder )
    // if later find setting vectors to higher numbers define ^ but this is usually used for bitwise XOR

    bool operator==(const Vector& other) const; // returns element-wise equality bool
    double norm() const;

};


// Matrix class

class Matrix
{
private:
    
    std::vector<double> data; //store as 1D vec then use index formula where index = row*columns + column

    size_t columns;
    size_t rows;

public:

    //Constructors
    explicit Matrix(size_t r, size_t c);
    Matrix(size_t r, size_t c,std::initializer_list<double> list);

    //Get dimensions 
    size_t numRows() const;
    size_t numCols() const; 

    //Get elements methods (use index formula in c++ file!) 
    double& operator()(size_t r, size_t c);
    const double& operator()(size_t r, size_t c) const;

    //Make I matrix (using static allows call without requiring existing object)
    static Matrix identity(size_t n);

    //Arithemtic operations
    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(double scalar) const;
    
    //Matrix specific properties
    Vector eig() const;
    Matrix multMat(const Matrix& other) const;
    Matrix outerProd(const Matrix& other) const;
    double det() const; // only for 2 and 3 square mat
    Matrix inv() const; // only for 2 and 3 square mat

};

Vector triDiagSolve(Matrix A,Vector b);

// MAYBE ADD GENERAL TENSOR CLASS LATER + EXPAND det,inv AND eig to n dimensions?
