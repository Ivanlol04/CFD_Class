#include "linalg.h"

// Constructor 
Vector::Vector(size_t n) : data(n) {}

// Initializer list constructor
Vector::Vector(std::initializer_list<double> list): data(list) {}

// Size
size_t Vector::size() const
{
    return data.size();
}

// Element access 1
double& Vector::operator[](size_t i)
{
    return data[i];
}

// Element access 2 (ref)
const double& Vector::operator[](size_t i) const
{
    return data[i];
}

// Arithmetic operators
Vector Vector::operator+(const Vector& other) const
{

    if (size() != other.size()){
        throw std::invalid_argument("Dimension mismatch in + operator");
    }

    Vector result(size());

    for (size_t i = 0; i < size(); i++)
    {
        result[i] = data[i]+ other[i];
    }

    return result;
};
    
Vector Vector::operator-(const Vector& other) const
{

    if (size() != other.size())
    {
        throw std::invalid_argument("Dimension mismatch in - operator");
    };
    
    Vector result(size());

    for (size_t i = 0; i < size(); i++)
    {
        result[i] = data[i] - other[i];
    }
    
    return result; 

};

Vector Vector::operator*(double scalar) const//
{

    Vector result(size());

    for (size_t i = 0; i < size(); i++)
    {
        result[i] = data[i] * scalar; 
    }
    
    return result;

}; 
    
double Vector::dot(const Vector& other) const
{

    if (size() != other.size())
    {
        throw std::invalid_argument("Dimension mismatch in dot product");
    }

    double result = 0.0;

    for (size_t i = 0; i < size(); i++)
    {
        result += data[i] * other[i];
    }

    return result;    

};

Vector Vector::cross(const Vector& other) const
{

    if (other.size()!=3 || size()!=3)
    {
        throw std::invalid_argument("Dimension error in cross product");
    }
    
    Vector result(3);

    result[0] = data[1]*other[2] - data[2]*other[1];
    result[1] = data[2]*other[0] - data[0]*other[2];
    result[2] = data[0]*other[1] - data[1]*other[0];

    return result;

};

Vector Vector::ewiseDiv(const Vector& other) const // not dividing all elements by constant but element-wise division (not standard but practical)
{

    if (size()!=other.size())
    {
        throw std::invalid_argument("Dimension mismatch");
    }
    
    Vector result(size());

    for (size_t i = 0; i < size(); i++)
    {
        result[i] = data[i]/other[i];
    }

    return result;

}; 

Vector Vector::ewiseMult(const Vector& other) const // element-wise multiplication (gets rid of remainder )
{

    if (size()!=other.size())
    {
        throw std::invalid_argument("Dimension mismatch");
    }

    Vector result(size());

    for (size_t i = 0; i < size(); i++)
    {
        result[i] = data[i]*other[i];
    }

    return result;
    
}; 

bool Vector::operator==(const Vector& other) const // returns element-wise equality bool, maybe change to accept additional tolerance parameter epsilon
{

    if (size()!=other.size())
    {
        throw std::invalid_argument("Dimension mismatch");
    }

    for (size_t i = 0; i < size(); i++)
    {
        if (data[i] != other[i])
        {
            return false;
        }
        
    }
    
    return true;
}; 

double Vector::norm() const
{

    double result = 0.0;

    for (size_t i = 0; i < size(); i++)
    {
    
        result += data[i] * data[i];

    }
    
    result = std::sqrt(result);

    return result;

};

void Vector::wipe()
{

    for (size_t i = 0; i < size(); i++)
    {
        data[i] = 0;
    }

};

//Constructors
Matrix::Matrix(size_t r, size_t c): data(r*c) {rows = r; columns = c;}

Matrix::Matrix(size_t r, size_t c,std::initializer_list<double> list): data(list)  {rows = r; columns = c;}

//Get dimensions 
size_t Matrix::numRows() const
{
    return rows;
};

size_t Matrix::numCols() const
{
    return columns;
}; 

//Get elements methods (use index formula in c++ file!) 
    
double& Matrix::operator()(size_t r, size_t c)
{

    return data[r*columns + c];

};
    
const double& Matrix::operator()(size_t r, size_t c) const
{

    return data[r*columns + c];

};
    
//Make I matrix (using static allows call without requiring existing object)    
Matrix Matrix::identity(size_t n)
{
    if (n<1)
    {
        throw std::invalid_argument("Identity matrix must be defined by positive nonzero integer");
    }

    Matrix I(n,n);

    for (size_t i = 0; i < n; i++)
    {
        for (size_t j = 0; j < n; j++)
        {
            if (i==j)
            {
                I(i,j) = 1;
            }
        }
    }
    
    return I;

    // Example of how to call the static method to create a 3×3 identity matrix in a different script
    //Matrix I = Matrix::identity(3);
};

//Arithemtic operations
Matrix Matrix::operator+(const Matrix& other) const
{

    if (columns != other.columns || rows != other.rows)
    {
        throw std::invalid_argument("Matricies have to have the same dimensions in arithmetic operations");
    }

    Matrix result(rows,columns);
    size_t index = 0;

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            index = i*columns + j;
            result.data[index] = data[index] + other.data[index];
        }
        
    }
    
    return result;

};

Matrix Matrix::operator-(const Matrix& other) const
{

    if (columns != other.columns || rows != other.rows)
    {
        throw std::invalid_argument("Matricies have to have the same dimensions in arithmetic operations");
    }

    Matrix result(rows,columns);
    size_t index = 0;

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            index = i*columns + j;
            result.data[index] = data[index] - other.data[index];
        }
        
    }
    
    return result;
    
};

Matrix Matrix::operator*(double scalar) const
{

    Matrix result(rows,columns);
    size_t index = 0;

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            index = i*columns + j;
            result.data[index] = data[index] *scalar;
        }
        
    }
    
    return result;

};

//Matrix specific properties
Matrix Matrix::multMat(const Matrix& other) const
{
    // Standard matrix multiplication: (rows x columns) * (other.rows x other.columns)
    if (columns != other.rows)
    {
        throw std::invalid_argument("Matrix multiplication dimension mismatch");
    }

    Matrix result(rows, other.columns);
    double sum = 0.0;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < other.columns; ++j)
        {
            for (size_t k = 0; k < columns; ++k)
            {
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;

            sum = 0.0;
        }
    }

    return result;
};

Matrix Matrix::outerProd(const Matrix& other) const
{
    // Kronecker product
    // Result size: (rows * other.rows) x (columns * other.columns)
    size_t r_res = rows * other.rows;
    size_t c_res = columns * other.columns;

    Matrix result(r_res, c_res);

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < columns; ++j)
        {
            double a = (*this)(i, j);
            for (size_t k = 0; k < other.rows; ++k)
            {
                for (size_t l = 0; l < other.columns; ++l)
                {
                    result(i * other.rows + k, j * other.columns + l) = a * other(k, l);
                }
            }
        }
    }

    return result;
};

Vector Matrix::eig() const
{
    // Simplified so only returns eigenvalues for a 2x2 matrix
    if (rows != 2 || columns != 2)
    {
        throw std::invalid_argument("Eigenvalues currently implemented only for 2x2 matrices");
    }

    double a = (*this)(0,0);
    double b = (*this)(0,1);
    double c = (*this)(1,0);
    double d = (*this)(1,1);

    double trace = a + d;
    double det = a*d - b*c;
    double discriminant = trace*trace - 4*det;

    // For simplicity, return real parts even if discriminant < 0 (could use std::complex, but Vector holds double)
    // Here we return two real eigenvalues (if discriminant >= 0) or the real parts if complex.
    double sqrt_disc = std::sqrt(std::abs(discriminant));
    Vector ev(2);
    ev[0] = (trace + sqrt_disc) / 2.0;
    ev[1] = (trace - sqrt_disc) / 2.0;

    return ev;
};

double Matrix::det() const
{
    if (rows != columns)
    {
        throw std::invalid_argument("Determinant defined only for square matrices");
    }

    if (rows == 2)
    {
        // 2x2 determinant: ad - bc
        return (*this)(0,0) * (*this)(1,1) - (*this)(0,1) * (*this)(1,0);
    }
    else if (rows == 3)
    {
        // 3x3 determinant using Sarrus' rule / cofactor expansion
        double a = (*this)(0,0), b = (*this)(0,1), c = (*this)(0,2);
        double d = (*this)(1,0), e = (*this)(1,1), f = (*this)(1,2);
        double g = (*this)(2,0), h = (*this)(2,1), i = (*this)(2,2);

        return a * (e*i - f*h) - b * (d*i - f*g) + c * (d*h - e*g);
    }
    else
    {
        throw std::invalid_argument("Determinant implemented only for 2x2 and 3x3 matrices");
    }
}

Matrix Matrix::inv() const
{
    if (rows != columns)
    {
        throw std::invalid_argument("Inverse defined only for square matrices");
    }

    double d = det();                      // take advantage of det()
    if (std::fabs(d) < 1e-12)               // tolerance for singularity
    {
        throw std::runtime_error("Matrix is singular (zero determinant) - cannot invert");
    }

    if (rows == 2)
    {
        // Inverse of 2x2: (1/det) * [ d  -b; -c  a ]
        Matrix inv(rows, columns);
        inv(0,0) =   (*this)(1,1) / d;
        inv(0,1) = - (*this)(0,1) / d;
        inv(1,0) = - (*this)(1,0) / d;
        inv(1,1) =   (*this)(0,0) / d;
        return inv;
    }
    else if (rows == 3)
    {
        // Inverse of 3x3 via adjugate matrix (transpose of cofactor matrix)
        double a = (*this)(0,0), b = (*this)(0,1), c = (*this)(0,2);
        double d = (*this)(1,0), e = (*this)(1,1), f = (*this)(1,2);
        double g = (*this)(2,0), h = (*this)(2,1), i = (*this)(2,2);

        // Compute cofactors
        double cof00 =  (e*i - f*h);
        double cof01 = -(d*i - f*g);
        double cof02 =  (d*h - e*g);
        double cof10 = -(b*i - c*h);
        double cof11 =  (a*i - c*g);
        double cof12 = -(a*h - b*g);
        double cof20 =  (b*f - c*e);
        double cof21 = -(a*f - c*d);
        double cof22 =  (a*e - b*d);

        // Adjugate is transpose of cofactor matrix
        Matrix inv(rows, columns);
        inv(0,0) =  cof00 / d;   inv(0,1) =  cof10 / d;   inv(0,2) =  cof20 / d;
        inv(1,0) =  cof01 / d;   inv(1,1) =  cof11 / d;   inv(1,2) =  cof21 / d;
        inv(2,0) =  cof02 / d;   inv(2,1) =  cof12 / d;   inv(2,2) =  cof22 / d;

        return inv;
    }
    else
    {
        throw std::invalid_argument("Inverse implemented only for 2x2 and 3x3 matrices");
    }
}


//Helper methods
//Use Thomas algorithm - only applicable to tridiagonal A matricies
Vector triDiagSolve(Matrix A,Vector r)
{

    // might need to check for zeros along main diagonal or although this should not be an issue for current HW

    size_t n = r.size();

    Vector x(n); //holds solutions

    Vector c(n); //new coefficient for forward sub
    Vector d(n); //new r after forward sub

    for (size_t i = 0; i < n; i++)
    {
        // c' = c/b for i = 1 and c' = c/(b-ac'_{i-1})
        if (i==0)
        {
            c[i] = A(i,i+1)/A(i,i);

            d[i] = r[i]/A(i,i);
        }
        else if (i<n-1)
        {
            c[i] = A(i,i+1)/(A(i,i) - A(i,i-1)*c[i-1]);
            d[i] = (r[i] - A(i,i-1)*d[i-1])/(A(i,i) - A(i,i-1)*c[i-1]);
        }
        else
        {        
        //d' = (r - a*d'_{i-1}) / (b - a*c'_{i-1})
        d[i] = (r[i] - A(i,i-1)*d[i-1])/(A(i,i) - A(i,i-1)*c[i-1]);        
        }
    }
    
    x[n-1] = d[n-1];

    for (size_t j = n-2; j > 0 ; j--)
    {
        x[j] = d[j] - c[j]*x[j+1];
    }

    x[0] = d[0] - c[0]*x[1];
    
    return x;

};

