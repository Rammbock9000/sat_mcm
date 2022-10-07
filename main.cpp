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
	std::vector<int> C;
	int timeout = 300;
	bool quiet = true;
	bool allow_negative_numbers = false;
	std::string solver_name = "no_solver";
	int threads = 1;
#ifdef USE_Z3
	solver_name = "z3";
#endif
#ifdef USE_CADICAL
	solver_name = "cadical";
#endif
	if (argc == 1) {
		std::cout << "Please call satscm like this: ./satscm <constant(s)> <solver_name> <timeout> <threads> <quiet>" << std::endl;
		std::cout << "If specifying multiple constants (for MCM instead of SCM), make sure to give them as a colon-separated list (duplicates, negative numbers and even numbers are getting pre-processed)" << std::endl;
		std::cout << "  e.g. './satscm 14709 z3 300 2 1' for constant 14709 with 300 sec time budget, 2 allowed CPU threads and without debug outputs using the Z3 backend" << std::endl;
		std::cout << "  e.g. './satscm 1:3:11:3:22 cadical 42 1 0' for constants 1, 3, 11, 3 and 22 (which gets transformed to constants 3 and 11) with 42 sec time budget, 1 allowed CPU thread and with debug outputs using the CaDiCaL backend" << std::endl;
		return 0;
	}
	if (argc > 1) {
		std::string s(argv[1]);
		try {
			std::stringstream c_str(s);
			std::string buff;
			while(std::getline(c_str, buff, ':')) {
				C.emplace_back(std::stoi(buff));
			}
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer(s)" << std::endl;
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
			threads = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 5) {
		std::string s(argv[5]);
		try {
			quiet = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 6) {
		std::string s(argv[6]);
		try {
			allow_negative_numbers = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	std::cout << "Starting OSCM for constant" << (C.size()>1?"s\n":" ");
	for (auto &c : C) {
		std::cout << (C.size()>1?"  ":"") << c << std::endl;
	}
	std::cout << "and " << timeout << " seconds timeout with solver " << solver_name << " and " << threads << " allowed threads" << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	if (solver_name == "cadical") {
#ifdef USE_CADICAL
		solver = std::make_unique<scm_cadical>(C, timeout, quiet, allow_negative_numbers);
#else
		throw std::runtime_error("Link CaDiCaL lib to use CaDiCaL backend");
#endif
	}
	else if (solver_name == "z3") {
#ifdef USE_Z3
		solver = std::make_unique<scm_z3>(C, timeout, quiet, threads, allow_negative_numbers);
#else
		throw std::runtime_error("Link Z3 lib to use Z3 backend");
#endif
	}
	else
		throw std::runtime_error("unknown solver name '"+solver_name+"'");
	solver->solve();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	std::cerr << "Finished solving after " << elapsed_time << " seconds" << std::endl;
	solver->print_solution();
	return 0;
}