#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <memory>
#include <cctype>
#include <algorithm>

#include <mcm.h>

#ifdef USE_CADICAL
#include <mcm_cadical.h>
#endif

#ifdef USE_Z3
#include <mcm_z3.h>
#endif

#ifdef USE_SYRUP
#include <mcm_syrup.h>
#endif

int main(int argc, char** argv) {
	std::unique_ptr<mcm> solver;
	std::vector<std::vector<int>> C;
	int timeout = 300;
	mcm::verbosity_mode verbosity = mcm::verbosity_mode::normal_mode;
	bool allow_negative_numbers = false;
	std::string solver_name = "no_solver";
	int threads = 1;
	bool also_minimize_full_adders = false;
	bool allow_node_output_shift = false;
	bool write_cnf = false;
	bool enumerate_all = false;
	int allow_coefficient_sign_inversion = 0;
	int min_num_adders = -1;
#ifdef USE_Z3
	solver_name = "z3";
#endif
#ifdef USE_SYRUP
	solver_name = "syrup";
#endif
#ifdef USE_CADICAL
	solver_name = "cadical";
#endif
	if (argc == 1) {
		std::cout << "Please call satmcm like this: ./satmcm <constant(s)> <solver name> <timeout> <threads> <quiet> <minimize full adders> <allow post adder right shfits> <allow negative coefficients> <write cnf files> <allow coefficient sign inversion> <min num adders> <enumerate all>" << std::endl;
		std::cout << "  => constant(s): <int:int:int>;<...>: colon-separated list of integers that should be computed" << std::endl;
		std::cout << "  => solver name: <string>: cadical, z3, syrup are supported" << std::endl;
		std::cout << "  => timeout: <uint>: number of seconds allowed per SAT instance" << std::endl;
		std::cout << "  => threads: <uint>: number of threads allowed to use" << std::endl;
		std::cout << "  => quiet: <0/1>: suppress debug outputs by setting this to 1" << std::endl;
		std::cout << "  => minimize full adders: <0/1>: minimize the full adder count for the optimal number of adders by setting this to 1" << std::endl;
		std::cout << "  => allow post adder right shifts: <0/1>: account for the optional right shift after the addition" << std::endl;
		std::cout << "  => allow negative coefficients: <0/1>: allow the use of negative coefficients to decrease the FA count" << std::endl;
		std::cout << "  => write cnf files: <0/1>: write all SAT programs to CNF files" << std::endl;
		std::cout << "  => allow coefficient sign inversion: <0/1/-1>: 1 - allow the SAT solver to invert the sign of ANY requested coefficient to reduce the FA count; -1 - only allow it if for negative requested coefficients; 0 - never allow it" << std::endl;
		std::cout << "  => min num adders: <uint>: minimum number of adders (default: 0)" << std::endl;
		std::cout << "  => enumerate all: <0/1>: enumerate all possible solutions for optimal adder count instead of only searching for the optimum (this mode ignores the setting for <minimize full adders>; only feasible if the problem size is small enough => consider setting a timeout)" << std::endl;
		return 0;
	}
	if (argc > 1) {
		std::string s(argv[1]);
		try {
			std::stringstream c_str(s);
            std::stringstream s_str;
			std::string vector_buff;
			std::string buff;
			std::vector<std::string> v;
			while(std::getline(c_str, vector_buff, ';')) {
			    v.emplace_back(vector_buff);
			}
			for(auto entry : v){
			    std::cout << entry << std::endl;
			}
			for(int i = 0; i < v.size(); i++){
                s_str << v[i];
                std::cout << v[i] << std::endl;
                std::vector<int> vector_row;

                while(std::getline(s_str, buff, ':')) {
                    std::cout << "before emplace" << std::endl;
                    vector_row.emplace_back(std::stoi(buff));
                    std::cout << "after emplace" << std::endl;
                }
                std::cout << vector_row.size() << std::endl;
                for(int j = 0; j < vector_row.size(); j++) {
                    std::cout << vector_row[j] << std::endl;
                }
                std::cout << "----"<< std::endl;
                C.emplace_back(vector_row);
                s_str.clear();
			}
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer(s)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	//changed from constant(s): <int:int:...>: to constants(s): <int:int:int;...>
	for(auto v : C){
	    for(auto c : v){
	        std::cout << c << std::endl;
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
			auto quiet = (bool)std::stoi(s);
			if (quiet) {
				verbosity = mcm::verbosity_mode::normal_mode;
			}
			else {
				verbosity = mcm::verbosity_mode::debug_mode;
			}
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
			also_minimize_full_adders = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 7) {
		std::string s(argv[7]);
		try {
			allow_node_output_shift = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 8) {
		std::string s(argv[8]);
		try {
			allow_negative_numbers = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 9) {
		std::string s(argv[9]);
		try {
			write_cnf = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 10) {
		std::string s(argv[10]);
		try {
			allow_coefficient_sign_inversion = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to int" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 11) {
		std::string s(argv[11]);
		try {
			min_num_adders = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to int" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 12) {
		std::string s(argv[12]);
		try {
			enumerate_all = (bool)std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}

    if (C[0].size() == 1 and C.size() == 1 ) {
        std::cout << "Starting SCM for constant(s) " << C[0].front() << std::endl;
    }
    else if (C[0].size() == 1){
        std::cout << "Starting MCM for constant(s)" << std::endl;
        for (auto &v : C) {
            for (auto &c : v)
                std::cout << c << std::endl;
        }
    }
    else{
        std::cout << "Starting CMM for constant(s)" << std::endl;
        for (auto &v : C) {
            std::cout << "<";
            for (auto &c : v) {
                std::cout << "  " << c;
            }
            std::cout << " >" << std::endl;
        }
    }
    /* old console output for SCM and MCM
	std::cout << "Starting CMM for constant(s)" << (C.size()>1?"s\n":" ");
	for (auto &v : C) {
	    for (auto &c : v) {
            std::cout << (C.size() > 1 ? "  " : "") << c << (C.size() > 1 ? "\n" : " ");
        }
	}
	*/
	std::cout << "and " << timeout << " seconds timeout with solver " << solver_name << " and " << threads << " allowed threads" << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	if (solver_name == "cadical") {
#ifdef USE_CADICAL
		solver = std::make_unique<mcm_cadical>(C, timeout, verbosity, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link CaDiCaL lib to use CaDiCaL backend");
#endif
	}
	else if (solver_name == "syrup" or solver_name == "glucose" or solver_name == "glucose-syrup") {
#ifdef USE_SYRUP
		solver = std::make_unique<mcm_syrup>(C, timeout, verbosity, threads, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link Glucose-Syrup lib to use syrup backend");
#endif
	}
	else if (solver_name == "z3") {
#ifdef USE_Z3
		solver = std::make_unique<mcm_z3>(C, timeout, verbosity, threads, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link Z3 lib to use Z3 backend");
#endif
	}
	else
		throw std::runtime_error("unknown solver name '"+solver_name+"'");
	solver->set_enumerate_all(enumerate_all);
	if (also_minimize_full_adders) solver->also_minimize_full_adders();
	if (allow_node_output_shift) solver->allow_node_output_shift();
	if (allow_coefficient_sign_inversion != 0) solver->ignore_sign(allow_coefficient_sign_inversion == -1);
	if (min_num_adders >= 0) solver->set_min_add(min_num_adders);
	solver->solve();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	std::cerr << "Finished solving after " << elapsed_time << " seconds" << std::endl;
	solver->print_solution();
	auto [a,b] = solver->solution_is_optimal();
	std::cerr << "#Add optimal = " << a << std::endl;
	std::cerr << "#FAs optimal = " << b << std::endl;
	return 0;
}
