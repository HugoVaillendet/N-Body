#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>

constexpr size_t N = 2;
constexpr double G = 1.0;
constexpr double epsilon = 1e-5;

constexpr double t_0 = 0.0;
constexpr double t_f = 100.0;
constexpr double dt = 0.001;

constexpr size_t M = static_cast<size_t>((t_f - t_0) / dt);

std::vector<double> t(M);
std::vector<double> E(M);

struct Bodies {

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> vx;
    std::vector<double> vy;
    std::vector<double> vz;
    std::vector<double> ax;
    std::vector<double> ay;
    std::vector<double> az;
    std::vector<double> m;

    void init() {
        this->x.resize(N);
        this->y.resize(N);
        this->z.resize(N);
        this->vx.resize(N);
        this->vy.resize(N);
        this->vz.resize(N);
        this->ax.resize(N);
        this->ay.resize(N);
        this->az.resize(N);
        this->m.resize(N);
    }
};

void acceleration(Bodies& b) {

    for (int i = 0; i < N; i++) {
        b.ax[i] = 0.0;
        b.ay[i] = 0.0;
        b.az[i] = 0.0;
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i==j) continue;

            double dx = b.x[j] - b.x[i];
            double dy = b.y[j] - b.y[i];
            double dz = b.z[j] - b.z[i];

            //Distance squared between i and j.
            double r2 = dx*dx + dy*dy + dz*dz + epsilon*epsilon;

            double inverse_r3 = 1.0 / (r2 * std::sqrt(r2));

            //Evaluate force once.
            double F =(G * b.m[j]) * inverse_r3;

            //Apply acceleration on the 3 components.
            b.ax[i] += F * dx;
            b.ay[i] += F * dy;
            b.az[i] += F * dz;
        }
    }
}

void leapfrog_integrator(Bodies& b) {

    for (int i = 0; i < N; i++) {
        //First speed half step.
        b.vx[i] = b.vx[i] + 0.5 * dt * b.ax[i];
        b.vy[i] = b.vy[i] + 0.5 * dt * b.ay[i];
        b.vz[i] = b.vz[i] + 0.5 * dt * b.az[i];
    }

    for (int i = 0; i < N; i++) {
        //Position recalculation at half time step
        b.x[i] = b.x[i] + dt * b.vx[i];
        b.y[i] = b.y[i] + dt * b.vy[i];
        b.z[i] = b.z[i] + dt * b.vz[i];
    }

    //Acceleration recalculation
    acceleration(b);

    for (int i = 0; i < N; i++) {
        //Second speed half step.
        b.vx[i] = b.vx[i] + 0.5 * dt * b.ax[i];
        b.vy[i] = b.vy[i] + 0.5 * dt * b.ay[i];
        b.vz[i] = b.vz[i] + 0.5 * dt * b.az[i];
    }
}

double energy_evaluation(Bodies& b) {

    double T = 0.0;
    double V = 0.0;

    for (int i = 0; i < N; i++) {

        T += 0.5 * b.m[i] * (b.vx[i]*b.vx[i] + b.vy[i]*b.vy[i] + b.vz[i]*b.vz[i]);

        for(int j = i + 1; j < N; j++) {

            double dx = b.x[j] - b.x[i];
            double dy = b.y[j] - b.y[i];
            double dz = b.z[j] - b.z[i];

            //Distance squared between i and j.
            double r2 = dx*dx + dy*dy + dz*dz + epsilon*epsilon;

            //-= here since the potential is negative.
            V -= G * (b.m[i] * b.m[j]) / std::sqrt(r2);
        }
    }

    return (T + V);
}

void initialize(Bodies& b) {

    b.init();

    /*{
        b.x[0] = -1.0;
        b.y[0] = 0.0;
        b.z[0] = 0.0;

        b.vx[0] = 0.412103;
        b.vy[0] = 0.283384;
        b.vz[0] = 0.0;

        b.ax[0] = 0.0;
        b.ay[0] = 0.0;
        b.az[0] = 0.0;

        b.x[1] = 1.0;
        b.y[1] = 0.0;
        b.z[1] = 0.0;

        b.vx[1] = 0.412103;
        b.vy[1] = 0.283384;
        b.vz[1] = 0.0;

        b.ax[1] = 0.0;
        b.ay[1] = 0.0;
        b.az[1] = 0.0;

        b.x[2] = 0.0;
        b.y[2] = 0.0;
        b.z[2] = 0.0;

        b.vx[2] = -2.0 * 0.412103;
        b.vy[2] = -2.0 * 0.283384;
        b.vz[2] = 0.0;

        b.ax[2] = 0.0;
        b.ay[2] = 0.0;
        b.az[2] = 0.0;

        b.m[0] = 10.0;
        b.m[1] = 10.0;
        b.m[2] = 10.0;
    }*/

    {

        b.m[0] = 10.0;
        b.m[1] = 10.0;

        b.x[0] = -50.0;
        b.y[0] = 0.0;
        b.z[0] = 0.0;

        b.vx[0] = 0.0;
        b.vy[0] = -0.2236;
        b.vz[0] = 0.0;

        b.ax[0] = 0.0;
        b.ay[0] = 0.0;
        b.az[0] = 0.0;

        b.x[1] = 50.0;
        b.y[1] = 0.0;
        b.z[1] = 0.0;

        b.vx[1] = 0.0;
        b.vy[1] = 0.2236;
        b.vz[1] = 0.0;

        b.ax[1] = 0.0;
        b.ay[1] = 0.0;
        b.az[1] = 0.0;
    }

    for (int i = 0; i < M; i ++) {
        t[i] = dt * static_cast<double>(i);
    }

    std::cout << "Initalization successfull !\n";
}

void simulation(Bodies& b, std::vector<double>& t, std::vector<double>& E, 
    const std::string& path1) {

    initialize(b);

    std::ofstream Energy(path1);

    for (int i = 0; i < t.size(); i++) {
        leapfrog_integrator(b);
        E[i] = energy_evaluation(b);
        Energy << E[i] << '\n';
    }

    std::cout << "Simulation successful !\n";
}

int main() {
    
    Bodies b;

    simulation(b, t, E, "energy.dat");

    return 0;
}