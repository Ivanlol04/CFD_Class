#include "include/linalg.h"
#include <iostream>
#include <numbers>
#include <fstream>
#include <string>
#include <iomanip>

// Analytical function
double analiticF(double x,double t)
{
    double f = std::sin(x) * std::exp(-t);

    return f;
};

//evaluates IC for specific grid size
Vector IC(double delX,size_t numPoints)
{

    double x0 = 0.0;

    // numPoints refers to the number of nodes
    Vector IC(numPoints);

    for (size_t i = 0; i < numPoints; i++)
    {
        IC[i] = std::sin(x0+i*delX);
    }
    
    return IC;

};

//part a)
double fEuler(Vector u,double delX, double delT,size_t index,size_t numPoints)
{
    double res = 0.0;

    if (index == 0 || index == numPoints-1)
    {
        u[index] = 0.0; //hardcoded BCs   
    }
    else
    {
        res = u[index] + ((delT)/(delX *delX))*(u[index+1]-u[index]*2+u[index-1]);
    }

    return res;

};

Vector bEuler(Vector b,Matrix A)
{


    size_t n = b.size();

    Vector u(n);

    u = triDiagSolve(A,b);

    return u;


};

Vector CrankNicholson(Vector b, Matrix A)
{

    size_t n = b.size();

    Vector u(n);

    u = triDiagSolve(A,b);

    return u;

};

double L1norm(Vector eps)
{

    double norm = 0.0;

    int N = eps.size();

    for (size_t i = 0; i < N; i++)
    {
     
        norm += std::fabs(eps[i]); 

    }
    
    norm = norm/(N);

    return norm;

};

double L2norm(Vector eps)
{

    double norm = 0.0;

    int N = eps.size();

    for (size_t i = 0; i < N; i++)
    {
     
        norm += eps[i]*eps[i]; 

    }

    norm = norm/(N);

    norm = std::sqrt(norm);

    return norm;
};

Matrix initializeAMatBEuler(size_t n, double delT, double delX)
{

    Matrix A(n,n);

    double a = -(1.0/(delX*delX));
    double b = 1.0/delT + 2.0/(delX*delX);
    double c = -(1.0/(delX*delX));

    //Initialize A matrix
    for (size_t i = 0; i < n; i++)
    {
        
        if (i==0)
        {
            A(i,i) = b;
            A(i,i+1) = c;
        }
        else if (i==n-1)
        {
            A(i,i-1) = a;
            A(i,i) = b;
        }
        else
        {
            A(i,i-1) = a;
            A(i,i) = b;
            A(i,i+1) = c;
        }
    }

    return A;

};

Matrix initializeAMatCN(size_t n, double delT, double delX)
{


    Matrix A(n,n);

    double a = -(1.0/(2.0*delX*delX));
    double b = 1.0/delT + 1.0/(delX*delX);
    double c = -(1.0/(2.0*delX*delX));

    //Initialize A matrix
    for (size_t i = 0; i < n; i++)
    {
        
        if (i==0)
        {
            A(i,i) = b;
            A(i,i+1) = c;
        }
        else if (i==n-1)
        {
            A(i,i-1) = a;
            A(i,i) = b;
        }
        else
        {
            A(i,i-1) = a;
            A(i,i) = b;
            A(i,i+1) = c;
        }
    }

    return A;

};

Vector initializeUvecCN(Vector u, double delX, double delT)
{

    size_t n = u.size();

    Vector U(n);

    for (size_t i = 0; i < n; i++)
    {
        if (i==0)
        {

            U[i] = u[i]/delT + 0.5*((u[i+1]-2*u[i])/(delX*delX));

        }
        else if (i==n-1)
        {

            U[i] = u[i]/delT + 0.5*((u[i-1]-2*u[i])/(delX*delX));

        }
        else
        {

            U[i] = u[i]/delT + 0.5*((u[i+1]-2*u[i]+u[i-1])/(delX*delX));

        }
        
        
    }
    

    return U;

};


int main(){

    double bc1 = 0.0;
    double bc2 = 0.0;
    
    constexpr double pi = 3.14159265358979323846; // using older compiler   
    double t = 0.0;

    //part b)
    double gamma1 = 1.0/2.0; //gamma = delT/delX^2
    double gamma2 = 1.0/6.0;    
    double gamma3 = 1.0/100.0;

    std::vector<double> gammas = {gamma1,gamma2,gamma3};

    double delX = (2.0*pi)/64.0;

    double delT1 = gamma1*delX*delX;
    double delT2 = gamma2*delX*delX;
    double delT3 = gamma3*delX*delX;

    double L1normTimeT1 = 0.0;
    double L2normTimeT1 = 0.0;

    std::vector<double> L1normTimesTimeT1;
    std::vector<double> L2normTimesTimeT1;

    std::vector<double> delTs = {delT1,delT2,delT3};

    size_t numPoints = (2*pi)/delX + 1;

    Vector eps(numPoints);
    
    Vector u(numPoints);

    Matrix A(numPoints,numPoints);

    Vector uLast(numPoints);

    Vector uAnalyticFE1(numPoints); // store three analytic solutions as depend on final t so how runs are made influences each
    Vector uAnalyticFE2(numPoints);
    Vector uAnalyticFE3(numPoints);
    Vector uAnalyticBE1(numPoints); 
    Vector uAnalyticBE2(numPoints);
    Vector uAnalyticBE3(numPoints);
    Vector uAnalyticCN1(numPoints); 
    Vector uAnalyticCN2(numPoints);
    Vector uAnalyticCN3(numPoints);
    std::vector<Vector> hist;

    // Use pointers so writes through analyticSols[i] reach the named Vectors above.
    std::vector<Vector*> analyticSols = {&uAnalyticFE1,&uAnalyticFE2,&uAnalyticFE3,
                                         &uAnalyticBE1,&uAnalyticBE2,&uAnalyticBE3,
                                         &uAnalyticCN1,&uAnalyticCN2,&uAnalyticCN3};
    
    bool flagSteadyState = 0;

    size_t iter = 0;
    size_t maxIter = 500000;

    u = IC(delX,numPoints); //IC function takes delX and numPoints as first and second input (same for u1,u2,u3 only delT changes)

    Vector uAt1(numPoints); // used to find L1 and L2 norms for part b) at time close to or equal t=1 second

    double tCrit = 0.0; //stores the value of the time at which u for the norm calculations is written

    double tReachedByAllPartb = 1.00238; //not sure why but this value seems to be great to evaluate around as all three different gammas pass through with least deviation

    //For forward euler
    for (size_t i = 0; i < 3; i++)
    {    

        while (iter<maxIter && flagSteadyState == 0)
        {
            
            uLast = u;

            for (size_t j = 0; j < numPoints; j++)
            {

                u[j] = fEuler(uLast,delX,delTs[i],j,numPoints);

            }

            t += delTs[i];
            
            if (tReachedByAllPartb - 0.000005<=t && t<=tReachedByAllPartb + 0.000005)
            {
                uAt1 = u;
                std::cout<<std::setprecision(15)<<"Part b) FE for gamma"<<gammas[i]<<" written at t = "<<t<<" ,delT = "<<delTs[i]<<std::endl;
                tCrit = t;
            }

            iter += 1;

            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        u = IC(delX,numPoints);
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < numPoints; k++)
        {
        
            (*analyticSols[i])[k] = analiticF(k*delX,tCrit);
        
        }
        
        eps = *analyticSols[i] - uAt1;
        L1normTimesTimeT1.push_back(L1norm(eps));
        L2normTimesTimeT1.push_back(L2norm(eps));

        t = 0.0;

    }    
 
    //For backward euler
    for (size_t i = 0; i < 3; i++)
    {    
        while (iter<maxIter && flagSteadyState == 0)
        {
            
            uLast = u;

            A = initializeAMatBEuler(numPoints,delTs[i],delX);

            Vector rhs = u * (1.0/delTs[i]);
            u = bEuler(rhs,A);

            t += delTs[i];

            if (tReachedByAllPartb - 0.000005<=t && t<=tReachedByAllPartb + 0.000005)
            {
                uAt1 = u;
                std::cout<<std::setprecision(15)<<"Part b) BE for gamma"<<gammas[i]<<" written at t = "<<t<<" ,delT = "<<delTs[i]<<std::endl;
                tCrit = t;
            }
            
            iter += 1;

            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        u = IC(delX,numPoints);
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < numPoints; k++)
        {
        
            (*analyticSols[i+3])[k] = analiticF(k*delX,tCrit);
        
        }
        eps = *analyticSols[i+3] - uAt1;
        L1normTimesTimeT1.push_back(L1norm(eps));
        L2normTimesTimeT1.push_back(L2norm(eps));

        t = 0.0;

    }   

    //For CN
    Vector U(numPoints);

    for (size_t i = 0; i < 3; i++)
    {    
        while (iter<maxIter && flagSteadyState == 0)
        {
            
            uLast = u;

            A = initializeAMatCN(numPoints,delTs[i],delX);

            U = initializeUvecCN(u,delX,delTs[i]);

            u = bEuler(U,A);

            t += delTs[i];
            
            if (tReachedByAllPartb - 0.000005<=t && t<=tReachedByAllPartb + 0.000005)
            {
                uAt1 = u;
                std::cout<<std::setprecision(15)<<"Part b) CN for gamma"<<gammas[i]<<" written at t = "<<t<<" ,delT = "<<delTs[i]<<std::endl;
                tCrit = t;
            }

            iter += 1;

            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        u = IC(delX,numPoints);
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < numPoints; k++)
        {
        
            (*analyticSols[i+6])[k] = analiticF(k*delX,tCrit);
        
        }
        eps = *analyticSols[i+6] - uAt1;
        L1normTimesTimeT1.push_back(L1norm(eps));
        L2normTimesTimeT1.push_back(L2norm(eps));

        t = 0.0;

    }   

    //part c)

    double delT = 0.0001;

    double delX1 = (2.0*pi)/64.0;
    double delX2 = (2.0*pi)/128.0;
    double delX3 = (2.0*pi)/256.0;

    std::vector<double> delXs = {delX1,delX2,delX3};

    size_t numPoints1 = 2*pi/delX1 + 1; //maybe do values in calc so to not carry any error when dividing and multiplying by pi
    size_t numPoints2 = 2*pi/delX2 + 1;
    size_t numPoints3 = 2*pi/delX3 + 1; 

    std::vector<size_t> allNumPoints = {numPoints1,numPoints2,numPoints3};

    Vector uAnalyticFE4(numPoints1);
    Vector uAnalyticFE5(numPoints2);
    Vector uAnalyticFE6(numPoints3);
    Vector uAnalyticBE4(numPoints1);
    Vector uAnalyticBE5(numPoints2);
    Vector uAnalyticBE6(numPoints3);
    Vector uAnalyticCN4(numPoints1);
    Vector uAnalyticCN5(numPoints2);
    Vector uAnalyticCN6(numPoints3);

    Vector uAt1np1(numPoints1);
    Vector uAt1np2(numPoints2);
    Vector uAt1np3(numPoints3);

    Vector epsAtnp1(numPoints1);
    Vector epsAtnp2(numPoints2);
    Vector epsAtnp3(numPoints3);

    std::vector<Vector*> uAt1diffNp = {&uAt1np1,&uAt1np2,&uAt1np3};

    std::vector<Vector*> epsAtdiffNp = {&epsAtnp1,&epsAtnp2,&epsAtnp3};

    std::vector<Vector*> analyticSols2 = {&uAnalyticFE4,&uAnalyticFE5,&uAnalyticFE6,
                                          &uAnalyticBE4,&uAnalyticBE5,&uAnalyticBE6,
                                          &uAnalyticCN4,&uAnalyticCN5,&uAnalyticCN6};    
    //For forward euler
    for (size_t i = 0; i < 3; i++)
    {    
        Vector u(allNumPoints[i]);
        u = IC(delXs[i], allNumPoints[i]);  
        Vector uLast(allNumPoints[i]);

        while (iter<maxIter && flagSteadyState == 0)
        {

            uLast = u;

            for (size_t j = 0; j < allNumPoints[i]; j++)
            {

                u[j] = fEuler(uLast,delXs[i],delT,j,allNumPoints[i]);

            }

            iter += 1;

            t += delT;

            if (iter == 10000)
            {
                *uAt1diffNp[i] = u;
                std::cout<<std::setprecision(15)<<"Part c) FE for delX = "<<delXs[i]<<" written at t = "<<t<<" ,delT = "<<delT<<std::endl;
            }
            
            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < allNumPoints[i]; k++)
        {
            (*analyticSols2[i])[k] = analiticF(k*delXs[i],1.0);
        }
        *epsAtdiffNp[i] = *analyticSols2[i] - *uAt1diffNp[i];
        L1normTimesTimeT1.push_back(L1norm(*epsAtdiffNp[i]));
        L2normTimesTimeT1.push_back(L2norm(*epsAtdiffNp[i]));
        
        t = 0.0;

    }    
 
    //For backward euler
    for (size_t i = 0; i < 3; i++)
    {    
        Vector u(allNumPoints[i]);
        u = IC(delXs[i], allNumPoints[i]);  
        Vector uLast(allNumPoints[i]);

        while (iter<maxIter && flagSteadyState == 0)
        {
            
            uLast = u;

            A = initializeAMatBEuler(allNumPoints[i],delT,delXs[i]);

            Vector rhs = u * (1.0/delT);
            u = bEuler(rhs,A);

            iter += 1;
            
            t += delT;

            if (iter == 10000)
            {
                *uAt1diffNp[i] = u;
                std::cout<<std::setprecision(15)<<"Part c) BE for delX = "<<delXs[i]<<" written at t = "<<t<<" ,delT = "<<delT<<std::endl;
            }

            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < allNumPoints[i]; k++)
        {
            (*analyticSols2[i+3])[k] = analiticF(k*delXs[i],1.0);
        }
        *epsAtdiffNp[i] = *analyticSols2[i+3] - *uAt1diffNp[i];
        L1normTimesTimeT1.push_back(L1norm(*epsAtdiffNp[i]));
        L2normTimesTimeT1.push_back(L2norm(*epsAtdiffNp[i]));
        
        t = 0.0;

    }   

    //For CN
    for (size_t i = 0; i < 3; i++)
    {    
        Vector u(allNumPoints[i]);
        u = IC(delXs[i], allNumPoints[i]); 
        Vector U(allNumPoints[i]);
        U = initializeUvecCN(u,delXs[i],delT);
        Vector uLast(allNumPoints[i]);

        while (iter<maxIter && flagSteadyState == 0)
        {
            
            uLast = u;

            A = initializeAMatCN(allNumPoints[i],delT,delXs[i]);

            U = initializeUvecCN(u,delXs[i],delT);

            u = bEuler(U,A);

            iter += 1;

            t += delT;

            if (iter == 10000)
            {
                *uAt1diffNp[i] = u;
                std::cout<<std::setprecision(15)<<"Part c) CN for delX = "<<delXs[i]<<" written at t = "<<t<<" ,delT = "<<delT<<std::endl;
            }
            
            if (uLast == u)
            {
                flagSteadyState = 1;
            }

        }
        // Store and reset u for next run
        hist.push_back(u);
        u.wipe();
        u = IC(delXs[i],allNumPoints[i]);
        iter = 0;
        flagSteadyState = 0;
        for (size_t k = 0; k < allNumPoints[i]; k++)
        {
            (*analyticSols2[i+6])[k] = analiticF(k*delXs[i],1.0);
        }
        *epsAtdiffNp[i] = *analyticSols2[i+6] - *uAt1diffNp[i];
        L1normTimesTimeT1.push_back(L1norm(*epsAtdiffNp[i]));
        L2normTimesTimeT1.push_back(L2norm(*epsAtdiffNp[i]));
        
        t = 0.0;

    }   

    // Use hist to find norms alongside uAnalyticX and Lnorm functions

    
    std::vector<double> L1norms(18, 0.0);
    std::vector<double> L2norms(18, 0.0);

    // Part b) — fixed size, direct subtraction is safe
    for (size_t j = 0; j < 9; ++j)
    {
        Vector eps = hist[j] - (*analyticSols[j]);
        L1norms[j] = L1norm(eps);
        L2norms[j] = L2norm(eps);
    }

    for (size_t j = 0; j < 9; ++j)
    {
        Vector eps = hist[j+9] - (*analyticSols2[j]);
        L1norms[j+9] = L1norm(eps);
        L2norms[j+9] = L2norm(eps);
    }


    auto vectorIsInvalid = [](const Vector& v) -> bool
    {
        // Only flag NaN to make sure no memory leaks or wrong runs
        for (size_t k = 0; k < v.size(); ++k)
            if (std::isnan(v[k])) return true;
        return false;
    };

    // CSV export 
    // Run labels: 18 entries matching the order solutions were pushed into hist.
    // Part b): scheme × gamma, Part c): scheme × grid
    const std::vector<std::string> runLabels = {
        // Part b) — Forward Euler
        "b_FE_g0.5", "b_FE_g0.167", "b_FE_g0.01",
        // Part b) — Backward Euler
        "b_BE_g0.5", "b_BE_g0.167", "b_BE_g0.01",
        // Part b) — Crank-Nicolson
        "b_CN_g0.5", "b_CN_g0.167", "b_CN_g0.01",
        // Part c) — Forward Euler
        "c_FE_dx64", "c_FE_dx128", "c_FE_dx256",
        // Part c) — Backward Euler
        "c_BE_dx64", "c_BE_dx128", "c_BE_dx256",
        // Part c) — Crank-Nicolson
        "c_CN_dx64", "c_CN_dx128", "c_CN_dx256"
    };

    // delX value for each of the 18 runs, in the same order as runLabels.
    // Part b) all use delX = 2pi/16; part c) cycles {2pi/64, 2pi/128, 2pi/256} per scheme.
    const std::vector<double> runDelX = {
        // Part b) — all three schemes, three gammas each
        delX, delX, delX,   // FE
        delX, delX, delX,   // BE
        delX, delX, delX,   // CN
        // Part c) — all three schemes, three grids each
        delX1, delX2, delX3,  // FE
        delX1, delX2, delX3,  // BE
        delX1, delX2, delX3   // CN
    };

    // delT value for each of the 18 runs, in the same order as runLabels.
    // Part b) delT = gamma*delX^2 varies per run; part c) all use fixed delT = 0.001.
    const std::vector<double> runDelT = {
        // Part b) — delT = gamma * delX^2, three gammas per scheme
        delT1, delT2, delT3,  // FE
        delT1, delT2, delT3,  // BE
        delT1, delT2, delT3,  // CN
        // Part c) — fixed delT = 0.001 for all grids and schemes
        delT, delT, delT,   // FE
        delT, delT, delT,   // BE
        delT, delT, delT    // CN
    };

    // 1. Norm summary — one row per run, skipped if the solution is invalid
    {
        std::ofstream normFile("out_norms.csv");
        normFile << "run,delX,delT,L1_norm_final,L2_norm_final,L1_norm_t1,L2_norm_t1\n";
        for (size_t i = 0; i < 18; ++i)
        {
            const std::string& part = (i < 9) ? "PART B" : "PART C";
            if (vectorIsInvalid(hist[i]))
            {
                std::cout << "Error in " << part << " run " << runLabels[i]
                          << ": solution contains NaN skipping norm row.\n";
                continue;
            }
            normFile << runLabels[i]           << ","
                     << runDelX[i]             << ","
                     << runDelT[i]             << ","
                     << L1norms[i]             << ","
                     << L2norms[i]             << ","
                     << L1normTimesTimeT1[i]   << ","
                     << L2normTimesTimeT1[i]   << "\n";
        }
    }

    // 2. Per-run solution CSVs — x, u_numeric, u_analytic, pointwise_error
    //    Part b) runs share delX_b; part c) runs use delXs[j % 3]
    for (size_t i = 0; i < 18; ++i)
    {
        const std::string& part = (i < 9) ? "PART B" : "PART C";

        // Validate numeric solution
        if (vectorIsInvalid(hist[i]))
        {
            std::cout << "Error in " << part << " run " << runLabels[i]
                      << ": solution contains NaN skipping CSV write.\n";
            continue;
        }

        // Validate analytic solution for this run
        const Vector& analytic = (i < 9) ? (*analyticSols[i])
                                          : (*analyticSols2[i - 9]);
        if (vectorIsInvalid(analytic))
        {
            std::cout << "Error in " << part << " run " << runLabels[i]
                      << ": analytic solution contains NaN skipping CSV write.\n";
            continue;
        }

        std::ofstream sol("out_" + runLabels[i] + ".csv");
        sol << "x,u_numeric,u_analytic,error\n";

        if (i < 9)
        {
            // Part b) — fixed delX
            for (size_t k = 0; k < numPoints; ++k)
            {
                double x   = k * delX;
                double err = hist[i][k] - (*analyticSols[i])[k];
                sol << x << "," << hist[i][k] << ","
                    << (*analyticSols[i])[k] << "," << err << "\n";
            }
        }
        else
        {
            // Part c) — variable delX; index within the 9-element block
            size_t j      = i - 9;
            double dxHere = delXs[j % 3];
            size_t npHere = (*analyticSols2[j]).size();
            for (size_t k = 0; k < npHere; ++k)
            {
                double x   = k * dxHere;
                double err = hist[i][k] - (*analyticSols2[j])[k];
                sol << x << "," << hist[i][k] << ","
                    << (*analyticSols2[j])[k] << "," << err << "\n";
            }
        }
    }

    std::cout << "Wrote out_norms.csv and 18 solution CSVs.\n";
    
    return 0;
}

