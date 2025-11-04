#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <omp.h>
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <nlohmann/json.hpp>
#include <yaml-cpp/yaml.h>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <stb_image.h>
#include <vtk-9.0/vtkRenderer.h>
#include <vtk-9.0/vtkRenderWindow.h>
#include <vtk-9.0/vtkRenderWindowInteractor.h>
#include <vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0/vtkSTLWriter.h>

using json = nlohmann::json;

namespace CoAnQi {
	namespace Physics {
		const double PI = 3.141592653589793;
		const double c = 3.0e8;
		const double G = 6.67430e-11;

		double Omega_g = 7.3e-16;
		double Mbh = 8.15e36;
		double dg = 2.55e20;

		double v_SCm = 0.99 * c;
		double rho_A = 1e-23;
		double rho_sw = 8e-21;
		double v_sw = 5e5;
		double QA = 1e-10;
		double Qs = 0.0;
		double kappa = 0.0005;
		double alpha = 0.001;
		double gamma = 0.00005;
		double delta_sw = 0.01;
		double epsilon_sw = 0.001;
		double delta_def = 0.01;
		double HSCm = 1.0;
		double UUA = 1.0;
		double eta = 1e-22;
		double k1 = 1.5, k2 = 1.2, k3 = 1.8, k4 = 2.0;
		double beta_i = 0.6;
		double rho_v = 6e-27;
		double C_concentration = 1.0;
		double f_feedback = 0.1;
		const double num_strings = 1e9;
		double Ts00 = 1.27e3 + 1.11e7;
		std::vector<std::vector<double>> g_mu_nu = {
		{1.0, 0.0, 0.0, 0.0},
		{0.0, -1.0, 0.0, 0.0},
		{0.0, 0.0, -1.0, 0.0},
		{0.0, 0.0, 0.0, -1.0}
		};

		struct CelestialBody {
			std::string name;
			double Ms;
			double Rs;
			double Rb;
			double Ts_surface;
			double omega_s;
			double Bs_avg;
			double SCm_density;
			double QUA;
			double Pcore;
			double PSCm;
			double omega_c;
		};

		double step_function(double r, double Rb) {
			return (r > Rb) ? 1.0 : 0.0;
		}

		double compute_Ereact(double t, double rho_SCm, double v_SCm, double rho_A, double kappa) {
			if (rho_A <= 0.0) throw std::runtime_error("Invalid rho_A value");
			return (rho_SCm * v_SCm * v_SCm / rho_A) * std::exp(-kappa * t);
		}

		double compute_mu_s(double t, double Bs, double omega_c, double Rs, double SCm_contrib = 1e3) {
			double Bs_t = Bs + 0.4 * std::sin(omega_c * t) + SCm_contrib;
			return Bs_t * std::pow(Rs, 3);
		}

		double compute_grad_Ms_r(double Ms, double Rs) {
			if (Rs <= 0.0) throw std::runtime_error("Invalid Rs value");
			return G * Ms / (Rs * Rs);
		}

		double compute_Bj(double t, double omega_c, double SCm_contrib = 1e3) {
			return 1e-3 + 0.4 * std::sin(omega_c * t) + SCm_contrib;
		}

		double compute_omega_s_t(double t, double omega_s, double omega_c) {
			return omega_s - 0.4e-6 * std::sin(omega_c * t);
		}

		double compute_mu_j(double t, double omega_c, double Rs, double SCm_contrib = 1e3) {
			double Bj = compute_Bj(t, omega_c, SCm_contrib);
			return Bj * std::pow(Rs, 3);
		}

		double compute_Ug1(const CelestialBody& body, double r, double t, double tn, double alpha, double delta_def, double k1) {
			if (r <= 0.0) throw std::runtime_error("Invalid r value");
			double mu_s = compute_mu_s(t, body.Bs_avg, body.omega_c, body.Rs);
			double grad_Ms_r = compute_grad_Ms_r(body.Ms, body.Rs);
			double defect = 1.0 + delta_def * std::sin(0.001 * t);
			return k1 * mu_s * grad_Ms_r * std::exp(-alpha * t) * std::cos(PI * tn) * defect;
		}

		double compute_Ug2(const CelestialBody& body, double r, double t, double tn, double k2, double QA, double delta_sw, double v_sw, double HSCm, double rho_A, double kappa) {
			if (r <= 0.0) throw std::runtime_error("Invalid r value");
			double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
			double S = step_function(r, body.Rb);
			double wind_mod = 1.0 + delta_sw * v_sw;
			return k2 * (QA + body.QUA) * body.Ms / (r * r) * S * wind_mod * HSCm * Ereact;
		}

		double compute_Ug3(const CelestialBody& body, double r, double t, double tn, double theta, double rho_A, double kappa, double k3) {
			double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
			double omega_s_t = compute_omega_s_t(t, body.omega_s, body.omega_c);
			double Bj = compute_Bj(t, body.omega_c);
			return k3 * Bj * std::cos(omega_s_t * t * PI) * body.Pcore * Ereact;
		}

		double compute_Ug4(double t, double tn, double rho_v, double C_concentration, double Mbh, double dg, double alpha, double f_feedback, double k4) {
			if (dg <= 0.0) throw std::runtime_error("Invalid dg value");
			double decay = std::exp(-alpha * t);
			double cycle = std::cos(PI * tn);
			return k4 * rho_v * C_concentration * Mbh / dg * decay * cycle * (1 + f_feedback);
		}

		double compute_Um(const CelestialBody& body, double t, double tn, double rj, double gamma, double rho_A, double kappa, double num_strings, double phi_hat) {
			if (rj <= 0.0) throw std::runtime_error("Invalid rj value");
			double Ereact = compute_Ereact(t, body.SCm_density, v_SCm, rho_A, kappa);
			double mu_j = compute_mu_j(t, body.omega_c, body.Rs);
			double decay = 1.0 - std::exp(-gamma * t * std::cos(PI * tn));
			double single = mu_j / rj * decay * phi_hat;
			return single * num_strings * body.PSCm * Ereact;
		}

		double compute_Ubi(double Ugi, double beta_i, double Omega_g, double Mbh, double dg, double epsilon_sw, double rho_sw, double UUA, double tn) {
			if (dg <= 0.0) throw std::runtime_error("Invalid dg value");
			double wind_mod = 1.0 + epsilon_sw * rho_sw;
			return -beta_i * Ugi * Omega_g * Mbh / dg * wind_mod * UUA * std::cos(PI * tn);
		}

		std::vector<std::vector<double>> compute_A_mu_nu(double tn, double eta, double Ts00) {
			std::vector<std::vector<double>> A = g_mu_nu;
			double mod = eta * Ts00 * std::cos(PI * tn);
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < 4; ++j) {
					A[i][j] += mod;
				}
			}
			return A;
		}

		double compute_FU(const CelestialBody& body, double r, double t, double tn, double theta) {
			try {
				double Ug1 = compute_Ug1(body, r, t, tn, alpha, delta_def, k1);
				double Ug2 = compute_Ug2(body, r, t, tn, k2, QA, delta_sw, v_sw, HSCm, rho_A, kappa);
				double Ug3 = compute_Ug3(body, r, t, tn, theta, rho_A, kappa, k3);
				double Ug4 = compute_Ug4(t, tn, rho_v, C_concentration, Mbh, dg, alpha, f_feedback, k4);
				double sum_Ugi = Ug1 + Ug2 + Ug3 + Ug4;

				double Ubi1 = compute_Ubi(Ug1, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
				double Ubi2 = compute_Ubi(Ug2, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
				double Ubi3 = compute_Ubi(Ug3, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
				double Ubi4 = compute_Ubi(Ug4, beta_i, Omega_g, Mbh, dg, epsilon_sw, rho_sw, UUA, tn);
				double sum_Ubi = Ubi1 + Ubi2 + Ubi3 + Ubi4;

				double Um = compute_Um(body, t, tn, body.Rb, gamma, rho_A, kappa, num_strings);

				auto A = compute_A_mu_nu(tn, eta, Ts00);
				double A_scalar = A[0][0] + A[1][1] + A[2][2] + A[3][3];

				return sum_Ugi + sum_Ubi + Um + A_scalar;
			}
			catch (const std::exception& e) {
				std::cerr << "Error in compute_FU for " << body.name << ": " << e.what() << std::endl;
				return 0.0;
			}
		}

		void output_json_params(const CelestialBody& body) {
			std::cout << "{" << std::endl;
			std::cout << " \"name\": \"" << body.name << "\"," << std::endl;
			std::cout << " \"SCm_density\": " << body.SCm_density << "," << std::endl;
			std::cout << " \"UA\": " << body.QUA << "," << std::endl;
			std::cout << " \"Qs\": " << Qs << std::endl;
			std::cout << "}" << std::endl;
		}

		std::vector<CelestialBody> load_bodies(const std::string& filename) {
			// Full implementation as before
			std::vector<CelestialBody> bodies;
			// ... load from JSON/YAML/CSV
			return bodies;
		}
	}

	namespace MUGE {
		struct ResonanceParams {
			// Full as before
		};

		struct MUGESystem {
			// Full as before
		};

		// All compute functions as before
		std::vector<MUGESystem> load_muge_systems(const std::string& filename) {
			// Full implementation
			std::vector<MUGESystem> systems;
			// ... load
			return systems;
		}
	}

	namespace Fluid {
		class FluidSolver {
			// Full class as before
		};
	}

	namespace Testing {
		void run_unit_tests() {
			// Full tests as before
		}
	}

	namespace Graphics3D {
		struct MeshData {
			std::vector<glm::vec3> vertices;
			std::vector<glm::vec3> normals;
			std::vector<glm::vec2> texCoords;
			std::vector<unsigned int> indices;
		};

		bool loadOBJ(const std::string& path, MeshData& mesh) {
			// Full as before
			return true;
		}

		void exportOBJ(const std::string& path, const MeshData& mesh) {
			// Full as before
		}

		GLuint loadTexture(const std::string& path) {
			// Full as before
			return 0;
		}

		class Shader {
			// Full as before
		};

		class Camera {
			// Full as before
		};

		struct SimulationEntity {
			glm::vec3 position;
			glm::vec3 velocity;
			MeshData model;
			void update(float dt) {
				position += velocity * dt;
			}
		};

		void renderMultiViewports(const std::vector<Camera>& cameras, const std::vector<SimulationEntity>& entities) {
			// Full as before
		}

		struct KeyPosition {
			glm::vec3 position;
			double timeStamp;
		};

		struct KeyRotation {
			glm::quat orientation;
			double timeStamp;
		};

		struct KeyScale {
			glm::vec3 scale;
			double timeStamp;
		};

		struct BoneInfo {
			int id;
			glm::mat4 offset;
		};

		class Bone {
			// Full as before
		};

		MeshData generateProceduralLandscape(int width, int height, float scale) {
			// Full as before
			MeshData mesh;
			return mesh;
		}

		MeshData extrudeMesh(const MeshData& base, float height) {
			// Full as before
			MeshData extruded;
			return extruded;
		}

		MeshData booleanUnion(const MeshData& mesh1, const MeshData& mesh2) {
			// Full as before
			MeshData unionMesh;
			return unionMesh;
		}

		void renderLaTeX(const std::string& latexCode, float x, float y) {
			// Full as before
		}
	}

	namespace Plugins {
		class SIMPlugin {
			// Full as before
		};
	}

	namespace Utils {
		void simulate_quasar_jet(double initial_velocity, const std::string& output_file = "") {
			// Full as before
		}

		void print_summary_stats(const std::vector<double>& values, const std::string& name) {
			// Full as before
		}

		std::vector<Graphics3D::SimulationEntity> populate_simulation_entities(const std::vector<MUGE::MUGESystem>& muge_systems) {
			// Full as before
			std::vector<Graphics3D::SimulationEntity> entities;
			return entities;
		}

		void initOpenGL(GLFWwindow** window) {
			// Full as before
		}
	}
}

int main(int argc, char** argv) {
	std::string input_file_bodies = "bodies.json";
	std::string input_file_muge = "muge.json";
	std::string output_file = "output.txt";
	std::string plugin_path = "plugin.dll";

	std::vector<double> fu_values, compressed_values, resonance_values;

	try {
		std::vector<CoAnQi::Physics::CelestialBody> bodies = CoAnQi::Physics::load_bodies(input_file_bodies);
		if (bodies.empty()) {
			CoAnQi::Physics::CelestialBody sun = { "Sun", 1.989e30, 6.96e8, 1.496e13, 5778.0, 2.5e-6, 1e-4, 1e15, 1e-11, 1.0, 1.0, 2 * CoAnQi::Physics::PI / (11.0 * 365.25 * 24 * 3600) };
			// ... add default bodies
			bodies = { sun /* , other bodies */ };
		}

		double r = 1e13;
		double t = 0.0;
		double tn = t;
		double theta = 0.0;

		for (const auto& body : bodies) {
			r = body.Rb;
			double FU = CoAnQi::Physics::compute_FU(body, r, t, tn, theta);
			fu_values.push_back(FU);
			std::cout << "Unified Field Strength (FU) for " << body.name << " at t=" << t << ", r=" << r << ": " << FU << " (normalized units)" << std::endl;
			// ... print individual components
			CoAnQi::Physics::output_json_params(body);
		}

		CoAnQi::Utils::print_summary_stats(fu_values, "FU");

		std::string velocity_csv = output_file.empty() ? "" : output_file + "_velocity.csv";
		CoAnQi::Utils::simulate_quasar_jet(CoAnQi::Physics::v_SCm, velocity_csv);

		std::vector<CoAnQi::MUGE::MUGESystem> muge_systems = CoAnQi::MUGE::load_muge_systems(input_file_muge);
		if (muge_systems.empty()) {
			// ... add default systems
		}

		for (const auto& sys : muge_systems) {
			double compressed_g = CoAnQi::MUGE::compute_compressed_MUGE(sys);
			compressed_values.push_back(compressed_g);
			double resonance_g = CoAnQi::MUGE::compute_resonance_MUGE(sys, CoAnQi::MUGE::ResonanceParams{});
			resonance_values.push_back(resonance_g);
			std::cout << "Compressed MUGE g for " << sys.name << ": " << compressed_g << " m/s2" << std::endl;
			std::cout << "Resonance MUGE g for " << sys.name << ": " << resonance_g << " m/s2" << std::endl;
		}

		CoAnQi::Utils::print_summary_stats(compressed_values, "Compressed MUGE");
		CoAnQi::Utils::print_summary_stats(resonance_values, "Resonance MUGE");

		// Output to file
		if (!output_file.empty()) {
			std::ofstream out(output_file);
			// ... write values
		}

		// 3D Visualization
		GLFWwindow* window;
		CoAnQi::Utils::initOpenGL(&window);
		std::vector<CoAnQi::Graphics3D::SimulationEntity> entities = CoAnQi::

}
}