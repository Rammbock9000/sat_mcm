#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <memory>
#include <cctype>
#include <algorithm>

#include <scm.h>

#ifdef USE_CADICAL
#include <scm_cadical.h>
#endif

#ifdef USE_Z3
#include <scm_z3.h>
#endif

int main(int argc, char** argv) {
	std::unique_ptr<scm> solver;
	int C = 42;
	int timeout = 300;
	bool quiet = true;
	int word_size = 0; // means "calculate necessary word size automatically"
	std::string solver_name = "no_solver";
#ifdef USE_Z3
	solver_name = "z3";
#endif
#ifdef USE_CADICAL
	solver_name = "cadical";
#endif
	if (argc == 1) {
		std::cout << "Please call satscm like this: ./satscm <constant> <solver_name> <timeout> <quiet>" << std::endl;
		std::cout << "  e.g. './satscm 14709 cadical 300 1' for constant 14709 with 300 sec time budget and without debug outputs using the cadical backend" << std::endl;
		return 0;
	}
	if (argc > 1) {
		std::string s(argv[1]);
		try {
			C = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 2) {
		std::string s(argv[2]);
		std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){return std::tolower(c);});
		solver_name = s;
	}
	if (argc > 3) {
		std::string s(argv[3]);
		try {
			timeout = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 4) {
		std::string s(argv[4]);
		try {
			quiet = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	std::cout << "Starting OSCM for C = " << C << " and timeout = " << timeout << " seconds with solver " << solver_name << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	if (solver_name == "cadical") {
#ifdef USE_CADICAL
		solver = std::make_unique<scm_cadical>(C, timeout, quiet, word_size);
#else
		throw std::runtime_error("Link CaDiCaL lib to use CaDiCaL backend");
#endif
	}
	else if (solver_name == "z3") {
#ifdef USE_Z3
		solver = std::make_unique<scm_z3>(C, timeout, quiet, word_size);
#else
		throw std::runtime_error("Link Z3 lib to use Z3 backend");
#endif
	}
	else
		throw std::runtime_error("unknown solver name '"+solver_name+"'");
	solver->solve();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	std::cout << "Finished solving after " << elapsed_time << " seconds" << std::endl;
	solver->print_solution();
	return 0;
}