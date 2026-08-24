#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <functional>

struct DLUResult {
    std::vector<double> D, L, U;
};

double computeSpectralR(std::vector<double> C, int n, int iterations, double tol) {
    std::vector<double> v(n, 1.0);
    double lambda = 0.0;

    for (int i = 0; i < iterations; ++i) {
        std::vector<double> w(n, 0.0);
        for (int r = 0; r < n; ++r) {
            for (int j = 0; j < n; ++j) {
                w[r] += C[r*n + j] * v[j];
            }
        }

        double sum_sq = 0.0;
        for (int j = 0; j < n; ++j) {
            sum_sq += w[j] * w[j];
        }
        double norm = std::sqrt(sum_sq);

        double rayleigh = 0.0;
        
        for (int j = 0; j < n; ++j) {
            rayleigh += v[j] * w[j];   // OLD v, before normalization update
        }
        // (v is already unit norm from the previous iteration, so no division needed)

        for (int j = 0; j < n; ++j) {
            v[j] = w[j] / norm;   // now update v for next round
        }

std::cout << "iter " << i << " | norm: " << norm << " | rayleigh: " << rayleigh << "\n";

        double diff = std::abs(norm - lambda);
        lambda = norm;

        if (diff < tol) break;
    }
    return lambda;
}

std::vector<double> forwardSub(std::vector<double> dl, 
    std::vector<double> rhs, int n) {
    std::vector<double> x(n);

    /*
    dl*x=rhs
    x[0] = rhs[0] / dl[0];
    x[1] = (rhs[1] - dl[3]*x[0]) / dl[4];
    x[2] = (rhs[2] - dl[6]*x[0] - dl[7]*x[1]) / dl[8];
    
    sum += dl[i*n+k]*x[k]
    x[i] = (rhs[i] - sum) / dl[i*n + i]
    */

    for (int i = 0; i < n; ++i) {
    double sum = 0.0;
    for (int k = 0; k < i; ++k) {
        sum += dl[i*n + k] * x[k];
    }
    x[i] = (rhs[i] - sum) / dl[i*n + i];
    }
    return x;
}

std::vector<double> computeCj(const DLUResult& dlu, int n) {
    std::vector<double> l = dlu.L;
    std::vector<double> u = dlu.U;
    std::vector<double> d = dlu.D;

    std::vector<double> Cj(n*n, 0);

    std::vector<double> lu;
    lu.reserve(n*n);
    std::transform(l.begin(), l.end(), u.begin(), 
    std::back_inserter(lu), std::plus<double>());

    for (int i = 0; i<n; ++i) {
        for (int j = 0; j<n; ++j) {
            Cj[i*n + j] = ((-1.0)*(lu[i*n + j]))/(d[i*n + i]);
        }
    }
    return Cj;
}

std::vector<double> computeCgs(const DLUResult& dlu, int n) {
    std::vector<double> l = dlu.L;
    std::vector<double> u = dlu.U;
    std::vector<double> d = dlu.D;

    std::vector<double> Cgs(n*n, 0);

    std::vector<double> dl;
    dl.reserve(n*n);
    std::transform(d.begin(), d.end(), l.begin(), 
    std::back_inserter(dl), std::minus<double>());

    for (int i = 0; i<n; ++i) {
        std::vector<double> rhs(n, 0);
        for (int j = 0; j<n; ++j) {
            rhs[j] = u[j*n + i];
        }
        std::vector<double> m = forwardSub(dl, rhs, n);

        for (int j = 0; j<n; ++j) {
            Cgs[j*n + i] = m[j];
        }
    }
    return Cgs;
}

DLUResult splitDLU(const std::vector<double>& A, int n) {
    std::vector<double> D(A.size());
    std::vector<double> L(A.size());
    std::vector<double> U(A.size());
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i==j) {
                D[i*n + i] = A[i*n + i];
            }
            if (i>j) {
                L[i*n + j] = A[i*n + j];
            }
            if (i<j) {
                U[i*n + j] = A[i*n + j];
            }
        }
    }
    return {D, L, U};
};

int main() {
    std::vector<double> A = {5.0, 1.0, 2.0,
                        1.0, 3.0, -1.0,
                        1.0, 1.0, 2.5};
    int n = sqrt(A.size());

    DLUResult result = splitDLU(A, n);
    std::vector<double> Cj = computeCj(result, n);
    std::vector<double> Cgs = computeCgs(result, n);

    for (int i = 0; i < Cj.size(); ++i) {
    std::cout << "Cj[" << i << "] = " << Cj[i] << "\n";}
    for (int i = 0; i < Cgs.size(); ++i) {
        std::cout << "Cgs[" << i << "] = " << Cgs[i] << "\n";}

    double rhoJ = computeSpectralR(Cj, n, 100, 1e-8);
    double rhoGS = computeSpectralR(Cgs, n, 100, 1e-8);

    std::cout << "Spectral radius (Jacobi):       " << rhoJ << "\n";
    std::cout << "Spectral radius (Gauss-Seidel): " << rhoGS << "\n";
    return 0;
}
