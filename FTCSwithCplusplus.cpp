#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>

#define _USE_MATH_DEFINES

#include <ctime>
#include <algorithm>
#include <string>
#include <vector>
#include <limits>


using namespace std;

const double M_PI = M_PI;

vector<double> setMeshGrid(double x0, double xN, double dx) {
    int numPoints = static_cast<int>((xN-x0)/dx) + 1;
    vector<double> x(numPoints);
    double stepSize;

    for (int i = 0; i < numPoints; i++) {
        x[i] = x0 + i * dx;
    }

    return x;
}




double calculateDiffusionCriterion(double alpha, double xStep, double tStep) {
    double diffusionCriterion = (alpha * tStep)/pow(xStep,2);
    if (diffusionCriterion > 0.5){
        cout << "Unstable FTCS Solution with Diffusion Coefficient = " << diffusionCriterion << endl;
    }
    else if (diffusionCriterion <= 0.5){
            cout << "Stable FTCS Solution with Diffusion Coefficient = " << diffusionCriterion << endl;
    }

    return diffusionCriterion;
}


void applyInitalConditions(vector<double>& u, const vector<double>& x){
    for (size_t i = 0; i < u.size(); ++i) {
        u[i] = 100.0 * sin(M_PI * x[i] / 2.0);
    }
}


void applyBoundaryConditions(vector<double>& u, const vector<double>& x){
   u[0] = 0;
   u[u.size() - 1] = 0;
}


void solver(double alpha, double dx, double dt, double tStop) {
    vector<double> x = setMeshGrid(0,2,dx);

    int meshSize = x.size();
    int timeSize = static_cast<int>(tStop / dt);

    vector<double> uKcurrent(meshSize);
    vector<double> uKnext(meshSize);

    applyInitalConditions(uKcurrent, x);
    applyBoundaryConditions(uKcurrent, x);
    double diffCrit = calculateDiffusionCriterion(alpha, dx, dt);

    for (int i = 0; i < timeSize; ++i) {
        for (int j = 1; j < meshSize - 1; ++j) {
            uKnext[j] = uKcurrent[j] + diffCrit * (uKcurrent[j+1] - 2 * uKcurrent[j] + uKcurrent[j-1]);
        }
        applyBoundaryConditions(uKnext, x);
        uKcurrent = uKnext;
    }

   string saveFileName = "d05solutionT0";
    ofstream outfile("X:/Projects/cfdTraining/1DHeatEqnSolverFTCS/Solutions/" + saveFileName + ".csv");
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file for writing." << endl;
        return;
    }

    outfile << "x,u\n"; //Add header to .csv output file
    for (int i = 0; i < meshSize; ++i) {
        outfile << x[i] << "," << uKcurrent[i] << "\n";
    }
    outfile.close();
    cout << "Final solution at t = " << tStop << " saved to " + saveFileName + ".csv" << endl;
}


int main() {
    cout << "Current Working Directory: " << filesystem::current_path() << endl;
    double k = 0.13;
    double c = 0.11;
    double rho = 7.8;
    double alpha = k/(c * rho);
    double dx = 0.1;
    double diffCrit = 0.5;
    double endTime = 0;
    double dt = diffCrit * pow(dx,2)/alpha;
    solver(alpha, dx, dt, endTime);


    return 0;
}
