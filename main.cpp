#include <iostream>
#include <scm_cadical.h>
#include <sstream>
#include <string>
#include <chrono>

int main(int argc, char** argv) {
	int C = 42;
	int timeout = 300;
	bool quiet = true;
	if (argc == 1) {
		std::cout << "Please call satscm like this: ./satscm <constant> <timeout> <quiet>" << std::endl;
		std::cout << "  e.g. ./satscm 699829 300 1 for constant 699829 with 300 sec time budget and without debug outputs" << std::endl;
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
		try {
			timeout = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 3) {
		std::string s(argv[3]);
		try {
			quiet = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	std::cout << "Starting OSCM for C = " << C << " and timeout = " << timeout << " seconds" << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	scm_cadical s(C, timeout, quiet);
	s.solve();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	std::cout << "Finished solving after " << elapsed_time << " seconds" << std::endl;
	s.print_solution();
	return 0;
}