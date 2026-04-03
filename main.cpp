#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include "json.hpp"

using json = nlohmann::json;

constexpr double F_KM_M = 1e3;
constexpr double G = 6.6743e-11;

constexpr double t_0 = 0.0;
constexpr double t_f = 20.0 * 31557600.0;
constexpr double dt = 800.0;

constexpr size_t M = static_cast<size_t>((t_f - t_0) / dt);

size_t N = 0;

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
};

void load_config(const std::string& path, Bodies& b) {

    std::ifstream config(path);
    json data = json::parse(config);

    N = data["bodies"].size();
    b.x.resize(N);
    b.y.resize(N);
    b.z.resize(N);
    b.vx.resize(N);
    b.vy.resize(N);
    b.vz.resize(N);
    b.ax.resize(N);
    b.ay.resize(N);
    b.az.resize(N);
    b.m.resize(N);
    
    for (auto& body : data["bodies"]) {
        int id = body["id"];

        b.x[id]  = body["position"]["x"].get<double>() * F_KM_M;
        b.y[id]  = body["position"]["y"].get<double>() * F_KM_M;
        b.z[id]  = body["position"]["z"].get<double>() * F_KM_M;

        b.vx[id] = body["velocity"]["vx"].get<double>() * F_KM_M;
        b.vy[id] = body["velocity"]["vy"].get<double>() * F_KM_M;
        b.vz[id] = body["velocity"]["vz"].get<double>() * F_KM_M;

        b.m[id] = body["mass"].get<double>();
    }

}

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
            double r2 = dx*dx + dy*dy + dz*dz;

            double inverse_r3 = (r2 * std::sqrt(r2));

            //Evaluate acceleration per unit of distance.
            double A =(G * b.m[j]) / inverse_r3;

            //Apply acceleration on the 3 components.
            b.ax[i] += A * dx;
            b.ay[i] += A * dy;
            b.az[i] += A * dz;
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
            double r2 = dx*dx + dy*dy + dz*dz;

            //-= here since the potential is negative.
            V -= G * (b.m[i] * b.m[j]) / std::sqrt(r2);
        }
    }

    return (T + V);
}

void initialize(Bodies& b, const std::string& path) {


    load_config(path, b);
    acceleration(b);

    for (int i = 0; i < M; i ++) {
        t[i] = dt * static_cast<double>(i);
    }

    std::cout << "Initalization successfull !\n";

}

void simulation(Bodies& b, std::vector<double>& t, std::vector<double>& E, 
                const std::string& path1) {

    std::ofstream Energy(path1);
    Energy.precision(std::numeric_limits<double>::max_digits10);
    Energy << std::scientific;
    
    //For each timestep we apply leapfrog integration and record energy.
    for (int i = 0; i < t.size(); i++) {
        leapfrog_integrator(b);
        E[i] = energy_evaluation(b);
        Energy << t[i] << " " << E[i] << '\n';
    }

    Energy.close();

    std::cout << "Simulation successful !\n";
}

std::vector<double> RG_evaluate(Bodies& b) {

    std::vector<double> RG(N-1);
    
    for (int i = 1; i < N; i++) {
        double dx = b.x[i] - b.x[0];
        double dy = b.y[i] - b.y[0];
        double dz = b.z[i] - b.z[0];

        double r2 = dx*dx + dy*dy + dz*dz;

        RG[i-1] = std::sqrt(r2);
    }
    return RG;
}

void export_data(Bodies& b, const std::string& path_in, const std::string& path_out) {

    std::vector<double> RG = RG_evaluate(b);

    std::ifstream config(path_in);
    json data = json::parse(config);

    std::ofstream File(path_out + ".dat");
    File.precision(std::numeric_limits<double>::max_digits10);
    File << std::scientific;

    for (auto& body : data["bodies"]) { 
        size_t idx = static_cast<size_t>(body["id"].get<int>());
        File << body["name"] << '\n';
        File << "Positions: ";
        File << b.x[idx] << " " << b.y[idx] << " " << b.z[idx] << '\n';
        File << "Velocities: ";
        File << b.vx[idx] << " " << b.vy[idx] << " " << b.vz[idx] << '\n';
        if (idx > 0) {
            File << "RG: " << RG[idx-1] << '\n';
        }

        File << '\n';
    }

    File.close();
}

int main() {
    
    Bodies b;
    initialize(b, "solar_system.json");
    export_data(b, "solar_system.json", "begin_sys");

    simulation(b, t, E, "energy.dat");
    export_data(b, "solar_system.json", "end_sys");

    return 0;
}