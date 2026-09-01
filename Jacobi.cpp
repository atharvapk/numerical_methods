#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>

struct jacobiResult {
    std::vector<double> x;
    int iterations;
    bool converged;
};

bool isDiagonallyDominant(const std::vector<double>& A, int n){
    //testing if the component with the largest absolute value is greater than the sum of other components' absolute value in the same row
    for (int i=0; i<n; ++i) {
        double diag = std::abs(A[i*n + i]);
        double off = 0.0;
        for (int j=0; j<n; ++j) {
            if (j!=i){
                off += std::abs(A[i*n + j]);
            }
        }
        if (off >= diag) {
            return false;
        }
    }
    return true;
}

double computeResidualNorm(const std::vector<double>& A,
                              const std::vector<double>& b,
                              const std::vector<double>& x) {
    int n = b.size();
    double max_r = 0.0;

    for (int i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (int j = 0; j < n; ++j) { 
            Ax_i += A[i*n + j] * x[j];
        }

        double r_i = std::abs(Ax_i - b[i]);
        if (r_i > max_r) max_r = r_i;
    }
    return max_r;
}

jacobiResult solveJacobi(const std::vector<double>& A, 
                            const std::vector<double>& b, 
                            std::vector<double> x, 
                            const int iterations, double tol){
    int n = b.size();
    std::vector<double> x_next(n);
    std::ofstream csv_log("jacobi_convergence.csv");
    csv_log << "iteration,delta_norm\n";
    
    /*
    b = Ax
    error = b - Ax
    xnext = error/diag
    so for each row find the sum of the residual across columns - i.e. sum residuals calculated at j, evaluated per i
    */

    for (int k=0; k<iterations; ++k) {
        double delta_norm = 0;
        for (int i=0; i<n; ++i) {
            double sum = 0;
            for (int j=0; j<n; ++j) {
                if (i!=j) {sum += (A[i*n + j])*(x[j]);}
            }
            x_next[i] = (b[i] - sum)/(A[i*n + i]);
            delta_norm += std::abs(x_next[i] - x[i]);
        }

        csv_log << (k+1) << "," << delta_norm << "\n";

        if (delta_norm <= tol) {
            return {x_next, k+1, true};
        }

        x = x_next;
        std::cout << k + 1 << ": " << x_next[0] << " & " << x_next[1] << " & " << x_next[2] 
        /* << " & " << computeResidualNorm(A, b, x_next)% */ << std::endl;
    }
    return {x, iterations, false};
}

int main() {
    std::vector<double> A = {5.0, 1.0, 2.0,
                            1.0, 3.0, -1.0,
                            1.0, 1.0, 2.5};
    std::vector<double> b = {8.0, 3.0, 4.5};
    std::vector<double> x0 = {1.1, 0.9, 1.00};

    bool dominant = isDiagonallyDominant(A, 3);
    std::cout << "Diagonally dominant?: " << (dominant ? "yes" : "no, cannot proceed") << "\n";

    if (dominant == true){
        jacobiResult result = solveJacobi(A, b, x0, 50, 1e-6);
        double residual = computeResidualNorm(A, b, result.x);

        for (int i = 0; i < result.x.size(); ++i) {
        std::cout << "x[" << i << "] = " << result.x[i] << "\n";}
        std::cout << "Residual norm: " << residual << "\n";
        std::cout << "Iterations: " << result.iterations << "\n";
    }
    return 0;
}