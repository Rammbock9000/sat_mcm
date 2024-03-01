#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <memory>
#include <cctype>
#include <algorithm>
#include <mcm_executable.h>

#include <mcm.h>

#ifdef USE_CADICAL
#include <mcm_cadical.h>
#endif

#ifdef USE_KISSAT
#include <mcm_kissat.h>
#endif

#ifdef USE_Z3
#include <mcm_z3.h>
#endif

#ifdef USE_SYRUP
#include <mcm_syrup.h>
#endif

int main(int argc, char **argv) {
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
	bool min_adder_depth = false;
	bool pipelining = false;
	bool eq_output_stages = false;
    std::string executable_binary;
    std::string executable_pre_cnf_params;
    std::string executable_post_cnf_params;
    std::string solver_log_filename = "temp.log";
    std::string solver_err_filename = "temp.err";
#ifdef USE_KISSAT
	solver_name = "kissat";
#endif
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
		std::cout << "Please call satmcm like this: ./satmcm \"constant(s)\" [arg1=val1 arg2=val2 ...]"
			      << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "How to specify the constant(s):" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
		std::cout << "  => constant(s) [mandatory]: \"int:int:int;...;...\" specify the constant matrix" << std::endl;
		std::cout << "      => separate columns with colons and rows with semicolons" << std::endl;
		std::cout<<R"(      => you need to put this argument into "..." on UNIX-based systems)" << std::endl;
		std::cout<<R"(      => this corresponds to "A:B:C" for SOP and "A;B;C" for MCM and "A" for SCM)" << std::endl;
		std::cout<<R"(      => e.g., specify "11:21;33:-44" to get an adder graph that computes both 11*x1+21*x2 and 33*x1-44*x2)"
			      << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Choose your optimization settings:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "  => minimize_full_adders=<0/1> [default=0]: minimize the full adder count for the optimal number of adders by setting this to 1"
                  << std::endl;
        std::cout << "  => allow_post_adder_right_shift=<0/1> [default=0]: account for the optional right shift after the addition"
                  << std::endl;
        std::cout << "  => allow_negative_coefficients=<0/1> [default=0]: allow the use of negative coefficients to decrease the FA count"
                  << std::endl;
        std::cout << "  => min_adder_depth=<0/1> [default=0]: force the solution to have minimum adder depth (useful for low-latency applications) => WORK IN PROGRESS"
                  << std::endl;
        std::cout << "  => pipelining=<0/1> [default=0]: let the solver optimize under the assumption that the adder graph will be fully pipelined after each adder stage => WORK IN PROGRESS"
                  << std::endl;
        std::cout << "  => equalize_output_stages=<0/1> [default=0]: force all outputs into the same pipeline stage (only relevant for pipelining) => WORK IN PROGRESS"
                  << std::endl;
        std::cout << "  => allow_coefficient_sign_inversion=<-1/0/1/2> [default=0]: 2 - generate coefficients EXACTLY as requested (e.g., you request a -5 with allow_negative_coefficients=1 and you get an adder graph for -5 even though it needs 1 adder more than an adder graph for +5); 1 - allow the SAT solver to invert the sign of ANY requested coefficient to reduce the FA count; -1 - only allow it for negative requested coefficients; 0 - never allow it"
                  << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Solver backend:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
		std::cout << "  => solver_name=<string> [default depends on installation]: kissat, cadical, z3, syrup, executable are currently supported" << std::endl;
		std::cout << "  => timeout=<uint> [default=300]: number of seconds allowed per SAT instance" << std::endl;
		std::cout << "  => threads=<uint> [default=1]: number of threads allowed to use" << std::endl;
		std::cout << "  => quiet=<0/1> [default=1]: enable/suppress debug outputs by setting this to 0/1" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "If you want to call some executable as your solver backend (e.g., some non-supported solver like Minisat):" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout<<R"(  => executable_binary=<string> [default=""]: path to binary used if solver_name=executable (execute an installed solver via system call instead of a C++ API); it should expect a .cnf file in DIMACS format and should produce an output "s (UN)SATISFIABLE" for (un)satisfiable instance in addition to solution literals as "v lit1 lit2 ..." according to the rules given by the SAT competition (url: satcompetition.org))" << std::endl;
        std::cout << "  => pre_cnf_params=<string> [default=\"\"]: command-line parameters that are given to the executable *before* the path to the cnf file" << std::endl;
        std::cout << "  => post_cnf_params=<string> [default=\"\"]: command-line parameters that are given to the executable *after* the path to the cnf file" << std::endl;
        std::cout << "  => solver_log_filename=<string> [default=\"temp.log\"]: the std::out file used for communication between this program and the executable solver" << std::endl;
        std::cout << "  => solver_err_filename=<string> [default=\"temp.err\"]: the std::err file used for communication between this program and the executable solver" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Some general settings:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
		std::cout << "  => write_cnf_files=<0/1>: write all SAT programs to CNF files [default=0]" << std::endl;
		std::cout << "  => min_num_adders=<uint> [default=0]: minimum number of adders" << std::endl;
		std::cout << "  => enumerate_all=<0/1> [default=0]: enumerate all possible solutions for optimal adder count instead of only searching for the optimum (this mode ignores the setting for <minimize full adders>; only feasible if the problem size is small enough => consider setting a timeout)"
			      << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Here's an example:" << std::endl;
        std::cout << "-----------------------------------------------------------------------------------" << std::endl;
        std::cout << "./satmcm \"43:11;7:23\" solver_name=cadical timeout=123 allow_negative_coefficients=1" << std::endl;
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
			while (std::getline(c_str, vector_buff, ';')) {
				v.emplace_back(vector_buff);
			}
			//for (auto &entry: v) {
			//	std::cout << entry << std::endl;
			//}
			for (int i = 0; i < v.size(); i++) {
				s_str << v[i];
				//std::cout << v[i] << std::endl;
				std::vector<int> vector_row;

				while (std::getline(s_str, buff, ':')) {
					//std::cout << "before emplace" << std::endl;
					vector_row.emplace_back(std::stoi(buff));
					//std::cout << "after emplace" << std::endl;
				}
				//std::cout << vector_row.size() << std::endl;
				//for (int j = 0; j < vector_row.size(); j++) {
				//	std::cout << vector_row[j] << std::endl;
				//}
				//std::cout << "----" << std::endl;
				C.emplace_back(vector_row);
				s_str.clear();
			}
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer(s) (arg #1)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	//changed from constant(s): <int:int:...>: to constants(s): <int:int:int;...>
	//for (auto &v: C) {
	//	for (auto c: v) {
	//		std::cout << c << std::endl;
	//	}
	//}

    for (int i=2; i<argc; i++) {
        std::string s_argv(argv[i]);
        //std::transform(s_argv.begin(), s_argv.end(), s_argv.begin(), [](unsigned char c) { return std::tolower(c); });
        std::stringstream ss;
        ss << s_argv;
        std::vector<std::string> arg_elements;
        std::string buffer;
        while (std::getline(ss, buffer, '=')) {
            arg_elements.emplace_back(buffer);
        }
        if (arg_elements.size() != 2) {
            std::cout << "UI WARNING: ignoring user argument '" << s_argv << "' -> only arguments in the form of arg=<value> are supported" << std::endl;
            continue;
        }
        std::string key = arg_elements.at(0);
        std::string val = arg_elements.at(1);
        std::string original_val = val; // copy it for case-sensitive parameters (e.g., path to executable_binary)
        std::transform(val.begin(), val.end(), val.begin(), [](unsigned char c) { return std::tolower(c); });
        if (key == "solver_name") {
            solver_name = val;
        }
        else if (key == "executable_binary") {
            executable_binary = original_val;
        }
        else if (key == "pre_cnf_params") {
            executable_pre_cnf_params = original_val;
        }
        else if (key == "post_cnf_params") {
            executable_post_cnf_params = original_val;
        }
        else if (key == "solver_log_filename") {
            solver_log_filename = original_val;
        }
        else if (key == "solver_err_filename") {
            solver_err_filename = original_val;
        }
        else if (key == "timeout") {
            try {
                timeout = std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "threads") {
            try {
                threads = std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "quiet") {
            try {
                auto quiet = (bool) std::stoi(val);
                if (quiet) {
                    verbosity = mcm::verbosity_mode::normal_mode;
                } else {
                    verbosity = mcm::verbosity_mode::debug_mode;
                }
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "minimize_full_adders") {
            try {
                also_minimize_full_adders = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "allow_post_adder_right_shift") {
            try {
                allow_node_output_shift = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "allow_negative_coefficients") {
            try {
                allow_negative_numbers = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "write_cnf_files") {
            try {
                write_cnf = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "allow_coefficient_sign_inversion") {
            try {
                allow_coefficient_sign_inversion = std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "min_num_adders") {
            try {
                min_num_adders = std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "enumerate_all") {
            try {
                enumerate_all = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "min_adder_depth") {
            try {
                min_adder_depth = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "pipelining") {
            try {
                pipelining = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else if (key == "equalize_output_stages") {
            try {
                eq_output_stages = (bool) std::stoi(val);
            }
            catch (...) {
                std::stringstream err_msg;
                err_msg << "invalid argument provided for '" << key << "'" << std::endl;
                throw std::runtime_error(err_msg.str());
            }
        }
        else {
            std::stringstream err_msg;
            err_msg << "argument of type '" << key << "=<value>' not supported" << std::endl;
            throw std::runtime_error(err_msg.str());
        }
    }

#if 0
	if (argc > 2) {
		std::string s(argv[2]);
		std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
		solver_name = s;
	}
	if (argc > 3) {
		std::string s(argv[3]);
		try {
			timeout = std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to integer (arg #3)" << std::endl;
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
			err_msg << "failed to convert " << s << " to integer (arg #4)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 5) {
		std::string s(argv[5]);
		try {
			auto quiet = (bool) std::stoi(s);
			if (quiet) {
				verbosity = mcm::verbosity_mode::normal_mode;
			} else {
				verbosity = mcm::verbosity_mode::debug_mode;
			}
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #5)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 6) {
		std::string s(argv[6]);
		try {
			also_minimize_full_adders = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #6)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 7) {
		std::string s(argv[7]);
		try {
			allow_node_output_shift = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #7)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 8) {
		std::string s(argv[8]);
		try {
			allow_negative_numbers = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #8)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 9) {
		std::string s(argv[9]);
		try {
			write_cnf = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #9)" << std::endl;
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
			err_msg << "failed to convert " << s << " to int (arg #10)" << std::endl;
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
			err_msg << "failed to convert " << s << " to int (arg #11)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 12) {
		std::string s(argv[12]);
		try {
			enumerate_all = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #12)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 13) {
		std::string s(argv[13]);
		try {
			min_adder_depth = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #13)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 14) {
		std::string s(argv[14]);
		try {
			pipelining = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #14)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
	if (argc > 15) {
		std::string s(argv[15]);
		try {
			eq_output_stages = (bool) std::stoi(s);
		}
		catch (...) {
			std::stringstream err_msg;
			err_msg << "failed to convert " << s << " to 1/0 (arg #15)" << std::endl;
			throw std::runtime_error(err_msg.str());
		}
	}
#endif

	if (C[0].size() == 1 and C.size() == 1) {
		std::cout << "Starting SCM for constant " << C[0].front() << std::endl;
	} else if (C[0].size() == 1) {
		std::cout << "Starting MCM for constants" << std::endl;
		for (auto &v: C) {
			for (auto &c: v)
				std::cout << c << std::endl;
		}
	} else if (C.size() == 1) {
		std::cout << "Starting SOP for matrix" << std::endl;
		for (auto &v: C) {
			std::cout << "<";
			for (auto &c: v) {
				std::cout << " " << c;
			}
			std::cout << " >" << std::endl;
		}
    } else {
        std::cout << "Starting CMM for matrix" << std::endl;
        for (auto &v: C) {
            std::cout << "<";
            for (auto &c: v) {
                std::cout << " " << c;
            }
            std::cout << " >" << std::endl;
        }
	}
	std::cout << "and " << timeout << " seconds timeout with solver " << solver_name << " and " << threads
						<< " allowed threads" << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	if (solver_name == "cadical") {
#ifdef USE_CADICAL
		solver = std::make_unique<mcm_cadical>(C, timeout, verbosity, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link CaDiCaL lib to use CaDiCaL backend");
#endif
	} else if (solver_name == "kissat") {
#ifdef USE_KISSAT
		solver = std::make_unique<mcm_kissat>(C, timeout, verbosity, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link kissat lib to use kissat backend");
#endif
	} else if (solver_name == "syrup" or solver_name == "glucose" or solver_name == "glucose-syrup") {
#ifdef USE_SYRUP
		solver = std::make_unique<mcm_syrup>(C, timeout, verbosity, threads, allow_negative_numbers, write_cnf);
#else
		throw std::runtime_error("Link Glucose-Syrup lib to use syrup backend");
#endif
	} else if (solver_name == "z3") {
#ifdef USE_Z3
        solver = std::make_unique<mcm_z3>(C, timeout, verbosity, threads, allow_negative_numbers, write_cnf);
#else
        throw std::runtime_error("Link Z3 lib to use Z3 backend");
#endif
    } else if (solver_name == "executable") {
        solver = std::make_unique<mcm_executable>(C, executable_binary, executable_pre_cnf_params, executable_post_cnf_params, solver_log_filename, solver_err_filename, verbosity, allow_negative_numbers, write_cnf);
	} else
		throw std::runtime_error("unknown solver name '" + solver_name + "'");
	solver->set_enumerate_all(enumerate_all);
	if (also_minimize_full_adders) solver->also_minimize_full_adders();
	if (allow_node_output_shift) solver->allow_node_output_shift();
	if (allow_coefficient_sign_inversion != 0) {
        // = 0: always implement the positive versions
        // = 1: apply to all coefficients
        // = -1: only apply to negative coefficients
        // = 2: always implement coefficients exactly as requested
        if (allow_coefficient_sign_inversion == 1 or allow_coefficient_sign_inversion == -1) {
            auto only_apply_to_negative_coeffs = allow_coefficient_sign_inversion == -1;
            solver->ignore_sign(only_apply_to_negative_coeffs);
        }
        else if (allow_coefficient_sign_inversion == 2) {
            solver->implement_signs_as_requested();
        }
        else {
            throw std::runtime_error("invalid value for allow_coefficient_sign_inversion");
        }
    }
	if (min_num_adders >= 0) solver->set_min_add(min_num_adders);
	if (min_adder_depth) solver->minimize_adder_depth();
	if (pipelining) solver->enable_pipelining();
	if (eq_output_stages) solver->equalize_output_stages();
	solver->solve();
	auto elapsed_time = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count()) / 1000.0;
	std::cerr << "Finished solving after " << elapsed_time << " seconds" << std::endl;
	solver->print_solution();
	auto [a, b] = solver->solution_is_optimal();
	std::cerr << "#Add optimal = " << a << std::endl;
	std::cerr << "#FAs optimal = " << b << std::endl;
	return 0;
}
