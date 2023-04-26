// Evgeny Bobkunov
// CS-03
// e.bobkunov@innopolis.university

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>

using namespace std;

#ifdef WIN32
#define GNUPLOT_NAME "D:\\gnuplot\\bin\\gnuplot -persist"
#else
#define GNUPLOT_NAME "gnuplot -persist"
#endif



// Computes the least square approximation for a given data set
vector <double> leastSquareApproximation(int m, vector<pair<double, double>>& data, int n) {
    // Construct the matrix A and vector b
    vector<vector<double>> A(m, vector<double>(n + 1));
    vector<double> b(m);
    for (int i = 0; i < m; i++) {
        double t = data[i].first;
        double bValue = data[i].second;
        for (int j = 0; j <= n; j++) {
            A[i][j] = pow(t, j);
        }
        b[i] = bValue;
    }

    cout << "A:" << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j <= n; j++) {
            cout << fixed << setprecision(4) << (A[i][j] == -0.0 ? 0.0 : A[i][j]) << " ";
        }
        cout << endl;
    }

    // Compute A_T * A and A_T * b
    vector<vector<double>> AT_A(n + 1, vector<double>(n + 1));
    vector<double> AT_b(n + 1);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            for (int k = 0; k < m; k++) {
                AT_A[i][j] += A[k][i] * A[k][j];
            }
        }
        for (int k = 0; k < m; k++) {
            AT_b[i] += A[k][i] * b[k];
        }
    }

    cout << "A_T*A:" << endl;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            cout << fixed << setprecision(4) << (AT_A[i][j] == -0.0 ? 0.0 : AT_A[i][j]) << " ";
        }
        cout << endl;
    }

    // Compute (A_T * A)^-1
    vector<vector<double>> AT_A_inv(n + 1, vector<double>(n + 1));
    for (int i = 0; i <= n; i++) {
        AT_A_inv[i][i] = 1;
    }
    for (int i = 0; i <= n; i++) {
        double pivot = AT_A[i][i];
        for (int j = 0; j <= n; j++) {
            AT_A[i][j] /= pivot;
            AT_A_inv[i][j] /= pivot;
        }
        for (int j = 0; j <= n; j++) {
            if (i != j) {
                double factor = AT_A[j][i];
                for (int k = 0; k <= n; k++) {
                    AT_A[j][k] -= factor * AT_A[i][k];
                    AT_A_inv[j][k] -= factor * AT_A_inv[i][k];
                }
            }
        }
    }

    cout << "(A_T*A)^-1:" << endl;
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            cout << fixed << setprecision(4) << (AT_A_inv[i][j] == -0.0 ? 0.0 : AT_A_inv[i][j]) << " ";
        }

        cout << endl;
    }

    cout << "A_T*b:" << endl;
    for (int i = 0; i <= n; i++) {
        cout << fixed << setprecision(4) << (AT_b[i] == -0.0 ? 0.0 : AT_b[i]) << endl;
    }


    // Compute x
    vector<double> x(n + 1);
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n; j++) {
            x[i] += AT_A_inv[i][j] * AT_b[j];
        }
    }


    cout << "x~:" << endl;
    for (int i = 0; i <= n; i++) {
        cout << fixed << setprecision(4) << (x[i] == -0.0 ? 0.0 : x[i]) << endl;
    }
    cout << endl;
    return x;
}

int main() {

#ifdef WIN32
    FILE* pipe = _popen(GNUPLOT_NAME, "w");
#else
    FILE* pipe = popen(GNUPLOT_NAME, "w");
#endif

    // Read input
    int m, n;
    cin >> m;
    vector<pair<double, double>> data(m);
    for (int i = 0; i < m; i++) {
        double t, b;
        cin >> t >> b;
        data[i] = make_pair(t, b);
    }
    cin >> n;
    // Compute the least square approximation
    vector<double> h = leastSquareApproximation(m, data, n);

    // Plot the least square approximation
    fprintf(pipe, "plot [-20 : 20] [-20 : 20]  %lf*x**3 + %lf*x**2 + %lf*x**1 + %lf*x**0 , '-' using 1:2 with points\n",  h[3], h[2], h[1], h[0]);
    for (int i = 0; i < m; i++) {
        fprintf(pipe, "%f\t%f\n", data[i].second, data[i].first);
    }
    fprintf(pipe, "e\n");
    fflush(pipe);

#ifdef WIN32
    _pclose(pipe);
#else
    pclose(pipe);
#endif
    return 0;
}
