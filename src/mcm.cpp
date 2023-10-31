//
// Created by nfiege on 9/26/22.
//

#include "mcm.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <numeric>

#define INPUT_SELECT_MUX_OPT 0 // I have NO IDEA WHY but apparently setting this to 0 is faster...
#define FPGA_ADD 0 // try out full adders as used in FPGAs ... maybe SAT solvers like those better than normal ones?!

mcm::mcm(const std::vector<std::vector<int>> &C, int timeout, verbosity_mode verbosity, int threads, bool allow_negative_numbers, bool write_cnf)
	:	C(C), timeout(timeout), verbosity(verbosity), threads(threads), write_cnf(write_cnf) {
	// make it even and count shift
    this->calc_twos_complement = allow_negative_numbers;
	for (auto &v : this->C) {
        // ignore 0 vector
        if (std::all_of(v.begin(), v.end(), [](int i) { return i == 0; }))
            continue; // c==0 before, now all values in v == 0
        auto original_vector = v;
        int shifted_bits = 0;
        // right shift until odd
        std::cout << "entry of vector has no odd value: " << std::all_of(v.begin(), v.end(), [](int i) { return ((i & 1) == 0); }) << std::endl;
        while (std::all_of(v.begin(), v.end(), [](int i) { return ((i & 1) == 0); })) {
            for (auto &c : v) {
                std::cout << "entry of vector: " << c << std::endl;
                if (c == 0) continue;
                c = c / 2; // do not use shift operation because it is not uniquely defined for negative number
                std::cout << "entry of vector after shift: " << c << std::endl;
            }
            shifted_bits++;
        }
        //look for first non zero value in vector and check if its > or < 0
        //note: do not remove if you think it's not needed think about it twice (norm for duplicate vector)
        for (auto &c : v) {
            std::cout << "found entry: " << c << std::endl;
            if (c > 0) {
                std::cout << "first non zero entry of vector is: " << c << std::endl;
                this->inverted_coeff_requested[v] = false;
                break;
            } else if (c < 0) {
                std::cout << "first non zero entry of vector is: " << c << std::endl;
                this->inverted_coeff_requested[v] = true;
                break;
            } else {
                continue;
            }
        }
        //flip the sign of all values in the vector if its requested
        if (this->inverted_coeff_requested[v]) {
        for (auto &c : v) {
            std::cout << "flipping sign of entry  " << c ;
            if (c == 0) continue;
            c = -c;
            std::cout << " ----> " << c << std::endl;
            }
        }
        this->requested_vectors[original_vector] = {v, shifted_bits};
	}
	// set word sizes & track unique constants
	this->word_size = 1;
	std::set<std::vector<int>> non_one_unique_vectors;
	std::vector<int> absV;
	for (auto &v : this->C) {
	    for (int i = 0; i < v.size(); i++){
            absV.emplace_back(abs(v[i]));
            std::cout << "abs constant: " << abs(v[i]) << std::endl;
	    }
        std::cout << "abs accumulate: " << std::accumulate(absV.begin(), absV.end(),0) << std::endl;
        //ignore all vector that contain only 0's and the unit vector where the sum of the absolute vector is 1
        std::cout << "return value 'not only contain 0's': " << !(std::all_of(v.begin(), v.end(), [](int i) { return i==0; })) << std::endl;
        std::cout << "return value 'sum of absolute vector is not 1': " << (std::accumulate(absV.begin(), absV.end(), 0) != 1)  << std::endl;

        if (!(std::all_of(v.begin(), v.end(), [](int i) { return i==0; })) and
            std::accumulate(absV.begin(), absV.end(), 0) != 1) non_one_unique_vectors.insert(v);
		//calculate ceiling over all values
		for (auto &c : v) {
            auto w = this->ceil_log2(std::abs(c)) + 1;
            if (w > this->word_size) this->word_size = w;
        }
		absV.clear();
	}

	this->max_shift = this->word_size-1;
	if (this->calc_twos_complement) {
		// account for sign bit
		this->word_size++;
	}

	this->shift_word_size = this->ceil_log2(this->max_shift+1);
	this->num_adders = (int)non_one_unique_vectors.size()-1;
	if (this->verbosity != verbosity_mode::quiet_mode) {
		std::cout << "Min num adders = " << this->num_adders + 1 << std::endl;
	}
	// set constants matrix
	this->C.clear();
	for (auto &v : non_one_unique_vectors) {
		this->C.emplace_back(v);
	}
    std::cout << "---------------" << std::endl;
    std::cout << "optimized input" << std::endl;
    for (auto &v : this->C) {
        std::cout << "|";
        for(auto &c : v){
            std::cout << "  " << c;
        }
        std::cout << " |" << std::endl;
    }
    std::cout << "---------------" << std::endl;
    std::cout << "word size: " << this->word_size-1  << std::endl;
    std::cout << "max shift: " << this->max_shift  << std::endl;
    std::cout << "number of adders: " << this->num_adders  << std::endl;
    std::cout << "shift word size: " << this->shift_word_size  << std::endl;
	//exit(0);
}

void mcm::optimization_loop(formulation_mode mode) {
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "  Starting optimization loop (mode = " << mode << ")" << std::endl;
	auto start_time = std::chrono::steady_clock::now();
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "  Resetting backend now" << std::endl;
	this->reset_backend(mode);
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "  Constructing problem for " << this->num_adders << " adders" << (this->max_full_adders!=FULL_ADDERS_UNLIMITED?" and "+std::to_string(this->max_full_adders)+" full adders":"") << std::endl;
	this->construct_problem(mode);
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "  Start solving with " << this->variable_counter << " variables and " << this->constraint_counter << " constraints" << std::endl;
    //this->create_arbitrary_clause({{1, false}});
    //this->create_arbitrary_clause({{1, true}});
	auto [a, b] = this->check();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	//throw std::runtime_error("DEBUG");
	this->fa_minimization_timeout -= elapsed_time;
	this->found_solution = a;
	this->ran_into_timeout = b;
	if (this->found_solution) {
		if (this->verbosity != verbosity_mode::quiet_mode) std::cout << "  Found solution for #adders = " << this->num_adders << (this->max_full_adders!=FULL_ADDERS_UNLIMITED?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " after " << elapsed_time << " seconds 8-)" << std::endl;
		this->get_solution_from_backend();
		if (this->solution_is_valid()) {
			if (this->verbosity != verbosity_mode::quiet_mode) std::cout << "Solution is verified :-)" << std::endl;
		}
		else {
			throw std::runtime_error("Solution is invalid (found bug) :-(");
		}
	}
	else if (this->ran_into_timeout) {
		if (this->verbosity != verbosity_mode::quiet_mode) std::cout << "  Ran into timeout for #adders = " << this->num_adders << (this->max_full_adders!=FULL_ADDERS_UNLIMITED?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " after " << elapsed_time << " seconds :-(" << std::endl;
	}
	else {
		if (this->verbosity != verbosity_mode::quiet_mode) std::cout << "  Problem for #adders = " << this->num_adders << (this->max_full_adders!=FULL_ADDERS_UNLIMITED?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " is proven to be infeasible after " << elapsed_time << " seconds... " << (this->max_full_adders!=FULL_ADDERS_UNLIMITED?"":"keep trying :-)") << std::endl;
	}
}

void mcm::solve() {
	if (this->enumerate_all) {
		this->solve_enumeration();
	}
	else {
		this->solve_standard();
	}
}

void mcm::reset_backend(formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	this->constraint_counter = 0;
	this->variable_counter = 0;
	this->cnf_clauses.str("");
}

void mcm::construct_problem(formulation_mode mode) {
	if (mode == formulation_mode::reset_all) {
		// only construct new variables in non-incremental mode
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "    creating variables now" << std::endl;
		this->create_variables();
	}
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "    creating constraints now" << std::endl;
	this->create_constraints(mode);
	if (this->write_cnf) {
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "    creating cnf file now" << std::endl;
		this->create_cnf_file();
	}
}

void mcm::create_variables() {
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "      creating input node variables" << std::endl;
	this->create_input_node_variables();
	for (int i=idx_input_buffer() + 1; i<=(this->num_adders + idx_input_buffer()); i++) {
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "      creating variables for node " << i << std::endl;
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_select_mux_variables" << std::endl;
		this->create_input_select_mux_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_select_selection_variables" << std::endl;
		this->create_input_select_selection_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_shift_value_variables" << std::endl;
		this->create_input_shift_value_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_shift_internal_variables" << std::endl;
		this->create_shift_internal_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_negate_select_variable" << std::endl;
		this->create_input_negate_select_variable(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_negate_select_output_variables" << std::endl;
		this->create_negate_select_output_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_negate_value_variable" << std::endl;
		this->create_input_negate_value_variable(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_xor_output_variables" << std::endl;
		this->create_xor_output_variables(i);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_adder_internal_variables" << std::endl;
		this->create_adder_internal_variables(i);
		if (this->enable_node_output_shift) {
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_post_adder_input_shift_value_variables" << std::endl;
			this->create_post_adder_input_shift_value_variables(i);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_post_adder_shift_variables" << std::endl;
			this->create_post_adder_shift_variables(i);
		}
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_output_value_variables" << std::endl;
		this->create_output_value_variables(i);
		if (this->c_column_size() != 1 or (this->calc_twos_complement and this->sign_inversion_allowed[this->C[0][0]])) {
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_mcm_output_variables" << std::endl;
			this->create_mcm_output_variables(i);
		}
		// full adder variables are constructed "on the fly" and put into their containers
	}
}

void mcm::create_constraints(formulation_mode mode) {
	if (this->verbosity == verbosity_mode::debug_mode) std::cout << "      create_input_output_constraints" << std::endl;
	this->create_input_output_constraints(mode);
	for (int i=idx_input_buffer() + 1; i<=(this->num_adders + idx_input_buffer()); i++) {
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "      creating constraints for node " << i << std::endl;
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_select_constraints" << std::endl;
		this->create_input_select_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_input_select_limitation_constraints" << std::endl;
		this->create_input_select_limitation_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_shift_limitation_constraints" << std::endl;
		this->create_shift_limitation_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_shift_constraints" << std::endl;
		this->create_shift_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_negate_select_constraints" << std::endl;
		this->create_negate_select_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_xor_constraints" << std::endl;
		this->create_xor_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_adder_constraints" << std::endl;
		this->create_adder_constraints(i, mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_odd_fundamentals_constraints" << std::endl;
		this->create_odd_fundamentals_constraints(i, mode);
		if (this->enable_node_output_shift) {
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_post_adder_shift_limitation_constraints" << std::endl;
			this->create_post_adder_shift_limitation_constraints(i, mode);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_post_adder_shift_constraints" << std::endl;
			this->create_post_adder_shift_constraints(i, mode);
		}
		if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_coeff_word_size_constraints" << std::endl;
			this->create_full_adder_coeff_word_size_constraints(i, mode);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_msb_constraints" << std::endl;
			this->create_full_adder_msb_constraints(i, mode);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_coeff_word_size_sum_constraints" << std::endl;
			this->create_full_adder_coeff_word_size_sum_constraints(i, mode);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_shift_gain_constraints" << std::endl;
			this->create_full_adder_shift_gain_constraints(i, mode);
			if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_shift_sum_constraints" << std::endl;
			this->create_full_adder_shift_sum_constraints(i, mode);
		}
	}
	if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_msb_sum_constraints" << std::endl;
		this->create_full_adder_msb_sum_constraints(mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_add_subtract_inputs_constraints" << std::endl;
		this->create_full_adder_add_subtract_inputs_constraints(mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_cpa_constraints" << std::endl;
		this->create_full_adder_cpa_constraints(mode);
		if (this->verbosity == verbosity_mode::debug_mode) std::cout << "        create_full_adder_result_constraints" << std::endl;
		this->create_full_adder_result_constraints();
	}
}

void mcm::create_input_node_variables() {
    for(int v=0; v<c_row_size(); v++) {
        for (int idx = 0; idx <= idx_input_buffer(); idx++) {
            for (int i = 0; i < this->word_size; i++) {
                this->output_value_variables[{idx, i, v}] = ++this->variable_counter;
                this->create_new_variable(this->variable_counter);
            }
        }
    }
}

void mcm::create_input_select_mux_variables(int idx) {
    if (idx < 2) return;
	if (idx <= idx_input_buffer()) return;
	//example 3x3 matrix: idx <= 2
	auto select_word_size = this->ceil_log2(idx);
#if INPUT_SELECT_MUX_OPT
	auto num_muxs = idx-1;
#else
	auto num_muxs = (1 << select_word_size) - 1;
#endif
    for(int v=0; v<c_row_size(); v++) {
        for (auto &dir : input_directions) {
            for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
                for (int w = 0; w < this->word_size; w++) {
                    this->input_select_mux_variables[{idx, dir, mux_idx, w, v}] = ++this->variable_counter;
                    this->create_new_variable(this->variable_counter);
#if INPUT_SELECT_MUX_OPT
                    if (mux_idx == num_muxs-1) {
                        this->input_select_mux_output_variables[{idx, dir, w}] = this->variable_counter;
                    }
#else
                    if (mux_idx == 0) {
                        this->input_select_mux_output_variables[{idx, dir, w, v}] = this->variable_counter;
                    }
#endif
                }
            }
        }
    }
}

void mcm::create_input_select_selection_variables(int idx) {
    if (idx == 1) return;
	if (idx <= idx_input_buffer()) return;
    //example 3x3 matrix: idx <= (1 + 2)
	auto select_word_size = this->ceil_log2(idx);
	for (auto &dir : input_directions) {
		for (int w = 0; w < select_word_size; w++) {
		    //std::cout << "idx: " << idx <<std::endl;
		    //std::cout << "dir: " << dir <<std::endl;
		    //std::cout << "w: " << w <<std::endl;
			this->input_select_selection_variables[{idx, dir, w}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

void mcm::create_input_shift_value_variables(int idx) {
	for (int w = 0; w < this->shift_word_size; w++) {
		this->input_shift_value_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void mcm::create_shift_internal_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (int mux_stage = 0; mux_stage < this->shift_word_size; mux_stage++) {
            for (int w = 0; w < this->word_size; w++) {
                this->shift_internal_mux_output_variables[{idx, mux_stage, w, v}] = ++this->variable_counter;
                if (mux_stage == this->shift_word_size - 1) {
                    this->shift_output_variables[{idx, w, v}] = this->variable_counter;
                }
                this->create_new_variable(this->variable_counter);
            }
        }
    }
}

void mcm::create_post_adder_input_shift_value_variables(int idx) {
	for (int w = 0; w < this->shift_word_size; w++) {
		this->input_post_adder_shift_value_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void mcm::create_post_adder_shift_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (int mux_stage = 0; mux_stage < this->shift_word_size; mux_stage++) {
            for (int w = 0; w < this->word_size; w++) {
                this->post_adder_shift_internal_mux_output_variables[{idx, mux_stage, w, v}] = ++this->variable_counter;
                if (mux_stage == this->shift_word_size - 1) {
                    this->post_adder_shift_output_variables[{idx, w, v}] = this->variable_counter;
                }
                this->create_new_variable(this->variable_counter);
            }
        }
    }
}

void mcm::create_input_negate_select_variable(int idx) {
	this->input_negate_select_variables[idx] = ++this->variable_counter;
	this->create_new_variable(this->variable_counter);
}

void mcm::create_negate_select_output_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (auto &dir : input_directions) {
            for (int w = 0; w < this->word_size; w++) {
                this->negate_select_output_variables[{idx, dir, w, v}] = ++this->variable_counter;
                this->create_new_variable(this->variable_counter);
            }
        }
    }
}

void mcm::create_input_negate_value_variable(int idx) {
	this->input_negate_value_variables[idx] = ++this->variable_counter;
	this->create_new_variable(this->variable_counter);
}

void mcm::create_xor_output_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            this->xor_output_variables[{idx, w, v}] = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
        }
    }
}

void mcm::create_adder_internal_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            this->adder_carry_variables[{idx, w, v}] = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
            this->adder_output_value_variables[{idx, w, v}] = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
#if FPGA_ADD
            this->adder_XOR_internal_variables[{idx, w}] = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
#endif
        }
    }
}

void mcm::create_output_value_variables(int idx) {
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            if (this->enable_node_output_shift) {
                this->output_value_variables[{idx, w, v}] = this->post_adder_shift_output_variables.at({idx, w, v});
            } else {
                this->output_value_variables[{idx, w, v}] = this->adder_output_value_variables.at({idx, w, v});
            }
        }
    }
}

int mcm::ceil_log2(int n) {
	try {
		return this->ceil_log2_cache.at(n);
	}
	catch (std::out_of_range&) {
		int val;
		if (n > 0) val = std::ceil(std::log2(n));
		else val = -1;
		return this->ceil_log2_cache[n] = val;
	}
}

int mcm::floor_log2(int n) {
	try {
		return this->floor_log2_cache.at(n);
	}
	catch (std::out_of_range&) {
		int val;
		if (n > 0) val = std::floor(std::log2(n));
		else val = -1;
		return this->floor_log2_cache[n] = val;
	}
}

void mcm::create_new_variable(int idx) {
	(void) idx; // just do nothing -> should be overloaded by backend if a variable must be explicitly created
}

void mcm::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	this->constraint_counter++;
	if (!this->write_cnf) return;
	for (const auto &it : a) {
		this->cnf_clauses << (it.second?-it.first:it.first) << " ";
	}
	this->cnf_clauses << " 0" << std::endl;
}

void mcm::create_signed_shift_overflow_protection(int sel, int s_a, int a) {
	// 1)
	this->create_arbitrary_clause({{sel, true}, {s_a, true}, {a, false}});
	// 2)
	this->create_arbitrary_clause({{sel, true}, {s_a, false}, {a, true}});
}

void mcm::create_signed_add_overflow_protection(int sub, int s_a, int s_b, int s_y) {
	// 1)
	this->create_arbitrary_clause({{sub, false}, {s_a, false}, {s_b, false}, {s_y, true}});
	// 2)
	this->create_arbitrary_clause({{sub, false}, {s_a, true}, {s_b, true}, {s_y, false}});
	// 3)
	this->create_arbitrary_clause({{sub, true}, {s_a, false}, {s_b, true}, {s_y, true}});
	// 4)
	this->create_arbitrary_clause({{sub, true}, {s_a, true}, {s_b, false}, {s_y, false}});
}

void mcm::create_or(std::vector<int> &x) {
	std::vector<std::pair<int, bool>> v(x.size());
	for (auto i=0; i<x.size(); i++) {
		auto &val = x[i];
        std::cout<<"create_or val: "<< val <<std::endl;
		if (val > 0) v[i] = {val, false};
		else v[i] = {-val, true};
	}
	this->create_arbitrary_clause(v);
}

void mcm::create_1x1_implication(int a, int b) {
	this->create_arbitrary_clause({{a, true}, {b, false}});
}

void mcm::create_1x1_negated_implication(int a, int b) {
	this->create_arbitrary_clause({{a, true}, {b, true}});
}

void mcm::create_1x1_reversed_negated_implication(int a, int b) {
	this->create_arbitrary_clause({{a, false}, {b, false}});
}

void mcm::create_1xN_implication(int a, const std::vector<int> &b) {
	std::vector<std::pair<int, bool>> v(b.size()+1);
	for (auto i=0; i<b.size(); i++) {
		v[i] = {b[i], false};
	}
	v[b.size()] = {a, true};
	this->create_arbitrary_clause(v);
}

void mcm::create_MxN_implication(const std::vector<int> &a, const std::vector<int> &b) {
	std::vector<std::pair<int, bool>> v(b.size()+a.size());
	for (auto i=0; i<a.size(); i++) {
		v[i] = {a[i], true};
	}
	for (auto i=0; i<b.size(); i++) {
		v[i+a.size()] = {b[i], false};
	}
	this->create_arbitrary_clause(v);
}

void mcm::create_1x1_equivalence(int x, int y) {
	// 1)
	this->create_arbitrary_clause({{x, true}, {y, false}});
	// 2)
	this->create_arbitrary_clause({{x, false}, {y, true}});
}

void mcm::create_2x1_mux(int a, int b, int s, int y) {
	// 1)
	this->create_arbitrary_clause({{a, true}, {s, false}, {y, false}});
	// 2)
	this->create_arbitrary_clause({{b, true}, {s, true}, {y, false}});
	// 3)
	this->create_arbitrary_clause({{b, false}, {s, true}, {y, true}});
	// 4)
	this->create_arbitrary_clause({{a, false}, {s, false}, {y, true}});
	// 5)
	this->create_arbitrary_clause({{a, true}, {b, true}, {y, false}});
	// 6)
	this->create_arbitrary_clause({{a, false}, {b, false}, {y, true}});
}

void mcm::create_2x1_mux_shift_disallowed(int a, int b, int s, int y) {
	// 1)
	this->create_arbitrary_clause({{a, true}, {y, false}});
	// 2)
	this->create_arbitrary_clause({{b, true}, {s, true}, {y, false}});
	// 3)
	this->create_arbitrary_clause({{b, false}, {s, true}, {y, true}});
	// 4)
	this->create_arbitrary_clause({{a, false}, {s, false}, {y, true}});
	// 5)
	this->create_arbitrary_clause({{a, true}, {s, true}});
	// 6)
	this->create_arbitrary_clause({{a, false}, {b, false}, {y, true}});
}

void mcm::create_2x1_mux_zero_const(int a, int s, int y) {
	// 1)
	this->create_arbitrary_clause({{s, true}, {y, true}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {y, true}});
	// 3)
	this->create_arbitrary_clause({{a, true}, {s, false}, {y, false}});
}

void mcm::create_2x1_xor(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, false}, {b, false}, {y, true}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {b, true}, {y, false}});
	// 3)
	this->create_arbitrary_clause({{a, true}, {b, false}, {y, false}});
	// 4)
	this->create_arbitrary_clause({{a, true}, {b, true}, {y, true}});
}

void mcm::create_2x1_equiv(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, false}, {b, false}, {y, false}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {b, true}, {y, true}});
	// 3)
	this->create_arbitrary_clause({{a, true}, {b, false}, {y, true}});
	// 4)
	this->create_arbitrary_clause({{a, true}, {b, true}, {y, false}});
}

void mcm::create_2x1_or(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, false},{b, false},{y, true}});
	// 2)
	this->create_arbitrary_clause({{a, true},{y, false}});
	// 3)
	this->create_arbitrary_clause({{b, true},{y, false}});
}

void mcm::create_2x1_and(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, true},{b, true},{y, false}});
	// 2)
	this->create_arbitrary_clause({{a, false},{y, true}});
	// 3)
	this->create_arbitrary_clause({{b, false},{y, true}});
}

void mcm::create_2x1_and_b_inv(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, true},{b, false},{y, false}});
	// 2)
	this->create_arbitrary_clause({{a, false},{y, true}});
	// 3)
	this->create_arbitrary_clause({{b, true},{y, true}});
}

void mcm::create_add_sum(int a, int b, int c_i, int s) {
	// 1)
	this->create_arbitrary_clause({{a, false}, {b, true}, {c_i, false}, {s, false}});
	// 2)
	this->create_arbitrary_clause({{a, true}, {b, false}, {c_i, false}, {s, false}});
	// 3)
	this->create_arbitrary_clause({{a, false}, {b, false}, {c_i, false}, {s, true}});
	// 4)
	this->create_arbitrary_clause({{a, true}, {b, true}, {c_i, false}, {s, true}});
	// 5)
	this->create_arbitrary_clause({{a, false}, {b, true}, {c_i, true}, {s, true}});
	// 6)
	this->create_arbitrary_clause({{a, true}, {b, false}, {c_i, true}, {s, true}});
	// 7)
	this->create_arbitrary_clause({{a, false}, {b, false}, {c_i, true}, {s, false}});
	// 8)
	this->create_arbitrary_clause({{a, true}, {b, true}, {c_i, true}, {s, false}});
}

void mcm::create_add_carry(int a, int b, int c_i, int c_o) {
	// 1)
	this->create_arbitrary_clause({{a, true}, {b, true}, {c_o, false}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {c_i, false}, {c_o, true}});
	// 3)
	this->create_arbitrary_clause({{b, false}, {c_i, false}, {c_o, true}});
	// 4)
	this->create_arbitrary_clause({{a, false}, {b, false}, {c_o, true}});
	// 5)
	this->create_arbitrary_clause({{b, true}, {c_i, true}, {c_o, false}});
	// 6)
	this->create_arbitrary_clause({{a, true}, {c_i, true}, {c_o, false}});
}

void mcm::create_add_redundant(int a, int b, int c_i, int s, int c_o) {
	// 1)
	this->create_arbitrary_clause({{a, false}, {s, true}, {c_o, true}});
	// 2)
	this->create_arbitrary_clause({{b, false}, {s, true}, {c_o, true}});
	// 3)
	this->create_arbitrary_clause({{c_i, false}, {s, true}, {c_o, true}});
	// 4)
	this->create_arbitrary_clause({{a, true}, {s, false}, {c_o, false}});
	// 5)
	this->create_arbitrary_clause({{b, true}, {s, false}, {c_o, false}});
	// 6)
	this->create_arbitrary_clause({{c_i, true}, {s, false}, {c_o, false}});
}

void mcm::force_bit(int x, int val) {
	this->create_arbitrary_clause({{x, val != 1}});
}

void mcm::forbid_number(const std::vector<int> &x, int val) {
	auto num_bits = (int)x.size();
	std::vector<std::pair<int, bool>> v(num_bits);
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			v[i] = {x[i], true};
		}
		else {
			v[i] = {x[i], false};
		}
	}
	this->create_arbitrary_clause(v);
}

void mcm::force_number(const std::vector<int> &x, int val) {
	auto num_bits = (int)x.size();
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			this->create_arbitrary_clause({{x[i], false}});
		}
		else {
			this->create_arbitrary_clause({{x[i], true}});
		}
	}
}

std::pair<bool, bool> mcm::check() {
	throw std::runtime_error("check is impossible in base class");
}

int mcm::get_result_value(int var_idx) {
	throw std::runtime_error("get_result_value is impossible in base class");
}

void mcm::create_input_output_constraints(formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	//inputs and outputs for SCM and MCM where input is only 1
	std::vector<int> input_bits(this->word_size);
	std::vector<int> output_bits(this->word_size);

	//non loop input/output_bits for SCM and parts of MCM and SOP
	for (auto w = 0; w < this->word_size; w++) {
	    input_bits[w] = this->output_value_variables.at({0, w, 0});
	    output_bits[w] = this->output_value_variables.at({this->num_adders, w, 0});
        //std::cout<<"input bits: "<< input_bits[w] <<std::endl;
        //std::cout<<"output bits: "<< output_bits[w] <<std::endl;
	}

	// force input to 1 and output to C[0][0]
	if (this->c_row_size() == 1 and this->c_column_size() == 1 and (!this->calc_twos_complement or !this->sign_inversion_allowed[this->C[0][0]])) {
		// SCM
        this->force_number(input_bits, 1);
		this->force_number(output_bits, this->C[0][0]);
	}
	else {
	    if(this->c_row_size() == 1) {
            // force input to 1 and build MCM output constraints
            // MCM
            this->force_number(input_bits, 1);
            this->create_mcm_output_constraints(mode);
        } else if(this->c_column_size() == 1) {
            // force input to unit vectors and outputs to C[0][0]
	        // SOP
            this->create_mcm_input_constraints(mode);
            this->force_number(output_bits, this->C[0][0]);
        } else {
            // force input to unit vectors and build CMM output constraints
            // CMM
            this->create_mcm_input_constraints(mode);
            this->create_mcm_output_constraints(mode);
        }
	}
}

void mcm::create_input_select_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	// stage 1 has no input MUX because it can only be connected to the input node with idx=0
    if (idx < 2) return;
	if (idx <= idx_input_buffer()) return;

	//create constraints for all muxs
	if (this->verbosity == verbosity_mode::debug_mode) {
		std::cout << "creating input select constraints for node #" << idx << std::endl;
	}
	auto select_word_size = this->ceil_log2(idx);
	auto next_pow_two = (1 << select_word_size);
    for(int v=0; v<c_row_size(); v++) {
	    for (auto &dir : this->input_directions) {
		    int mux_idx = 0;
		    std::map<std::pair<int, int>, int> signal_variables;
		    for (int i=0; i<idx; i++) {
			    for (int w=0; w<this->word_size; w++) {
				    signal_variables[{i, w}] = this->output_value_variables.at({i, w, v});
			    }
		    }
		    std::map<std::pair<int, int>, int> next_signal_variables;
		    auto num_signals = idx;
		    auto next_num_signals = 0;
#if INPUT_SELECT_MUX_OPT
		for (int mux_stage = 0; mux_stage < select_word_size; mux_stage++) {
			auto num_muxs_per_stage = next_pow_two >> (mux_stage+1);
			for (int mux_idx_per_stage = 0; mux_idx_per_stage < num_muxs_per_stage; mux_idx_per_stage++) {
				if (num_signals >= 2*(mux_idx_per_stage+1)) {
					// connect two signals from last stage to mux
					auto select_signal = this->input_select_selection_variables.at({idx, dir, mux_stage});
					for (int w = 0; w < this->word_size; w++) {
						auto zero_input = signal_variables.at({2*mux_idx_per_stage, w});
						auto one_input = signal_variables.at({2*mux_idx_per_stage+1, w});
						auto mux_output = this->input_select_mux_variables.at({idx, dir, mux_idx, w});
						next_signal_variables[{mux_idx_per_stage, w}] = mux_output;
						this->create_2x1_mux(zero_input, one_input, select_signal, mux_output);
					}
					next_num_signals++;
					mux_idx++;
				}
				else if (num_signals == 2*(mux_idx_per_stage+1)-1) {
					// only 1 signal left -> use it as an input to the next stage
					for (int w = 0; w < this->word_size; w++) {
						next_signal_variables[{mux_idx_per_stage, w}] = signal_variables.at({2*mux_idx_per_stage, w});
					}
					next_num_signals++;
				}
				else {
					// handled all signals of this stage
					break;
				}
			}
			// update counter & container
			num_signals = next_num_signals;
			signal_variables = next_signal_variables;
		}
#else
		    for (int mux_stage = 0; mux_stage < select_word_size; mux_stage++) {
			    auto num_muxs_per_stage = (1 << mux_stage);
			    auto mux_select_var_idx = this->input_select_selection_variables.at({idx, dir, select_word_size-mux_stage-1});
			    for (int mux_idx_in_stage = 0; mux_idx_in_stage < num_muxs_per_stage; mux_idx_in_stage++) {
				    if (mux_stage == select_word_size-1) {
					    // connect with another node output
					    auto zero_input_node_idx = 2 * mux_idx_in_stage;
					    auto one_input_node_idx = zero_input_node_idx + 1;
					    if (zero_input_node_idx >= idx) zero_input_node_idx = idx-1;
					    if (one_input_node_idx >= idx) one_input_node_idx = idx-1;
					    for (int w = 0; w < this->word_size; w++) {
						    auto mux_output_var_idx = this->input_select_mux_variables.at({idx, dir, mux_idx, w, v});
						    auto zero_input_var_idx = this->output_value_variables.at({zero_input_node_idx, w, v});
						    auto one_input_var_idx = this->output_value_variables.at({one_input_node_idx, w, v});
						    if (zero_input_node_idx == one_input_node_idx) {
							    // both inputs are equal -> mux output == mux input (select line does not matter...)
							    this->create_1x1_equivalence(zero_input_var_idx, mux_output_var_idx);
						    }
						    else {
							    this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, mux_select_var_idx, mux_output_var_idx);
						    }
					    }
				    }
				    else {
					    // connect with mux from higher stage
					    auto num_muxs_in_next_stage = (1 << (mux_stage + 1));
					    auto zero_mux_idx_in_next_stage = 2 * mux_idx_in_stage;
					    auto zero_input_mux_idx = num_muxs_in_next_stage - 1 + zero_mux_idx_in_next_stage;
					    auto one_input_mux_idx = zero_input_mux_idx + 1;
					    for (int w = 0; w < this->word_size; w++) {
						    auto mux_output_var_idx = this->input_select_mux_variables.at({idx, dir, mux_idx, w, v});
						    auto zero_input_var_idx = this->input_select_mux_variables.at({idx, dir, zero_input_mux_idx, w, v});
						    auto one_input_var_idx = this->input_select_mux_variables.at({idx, dir, one_input_mux_idx, w, v});
						    this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, mux_select_var_idx, mux_output_var_idx);
					    }
				    }
				    // increment current mux idx
				    mux_idx++;
			    }
		    }
#endif
	    }
    }
}

void mcm::create_shift_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
    for(int v=0; v<c_row_size(); v++) {
        for (auto stage = 0; stage < this->shift_word_size; stage++) {
            auto shift_width = (1 << stage);
            auto select_input_var_idx = this->input_shift_value_variables.at({idx, stage});
            auto first_disallowed_shift_bit = this->word_size - shift_width;
            for (auto w = 0; w < this->word_size; w++) {
                auto w_prev = w - shift_width;
                auto connect_zero_const = w_prev < 0;
                int zero_input_var_idx;
                int zero_input_sign_bit_idx;
                int one_input_var_idx;
                auto mux_output_var_idx = this->shift_internal_mux_output_variables.at({idx, stage, w, v});
                if (stage == 0) {
                    // connect shifter inputs
                    if (idx <= 1){
                        // shifter input is the output of the input node with idx = 0
                        zero_input_var_idx = this->output_value_variables.at({0, w, v});
                        zero_input_sign_bit_idx = this->output_value_variables.at({0, this->word_size - 1, v});
                        if (!connect_zero_const) {
                            one_input_var_idx = this->output_value_variables.at({0, w_prev, v});
                        }
                    } else {
                        // shifter input is the left input value
                        zero_input_var_idx = this->input_select_mux_output_variables.at({idx, mcm::left, w, v});
                        zero_input_sign_bit_idx = this->input_select_mux_output_variables.at(
                                {idx, mcm::left, this->word_size - 1, v});
                        if (!connect_zero_const) {
                            one_input_var_idx = this->input_select_mux_output_variables.at({idx, mcm::left, w_prev, v});
                        }
                    }
                } else {
                    // connect output of previous stage
                    zero_input_var_idx = this->shift_internal_mux_output_variables.at({idx, stage - 1, w, v});
                    zero_input_sign_bit_idx = this->shift_internal_mux_output_variables.at(
                            {idx, stage - 1, this->word_size - 1, v});
                    if (!connect_zero_const) {
                        one_input_var_idx = this->shift_internal_mux_output_variables.at({idx, stage - 1, w_prev, v});
                    }
                }
                if (w >= first_disallowed_shift_bit) {
                    if (this->calc_twos_complement) {
                        if (connect_zero_const) {
                            this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx,
                                                            mux_output_var_idx);
                        } else {
                            this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx,
                                                 mux_output_var_idx);
                        }
                        if (w == this->word_size - 1) {
                            // these clauses are different for the sign bit
                            this->create_signed_shift_overflow_protection(select_input_var_idx, zero_input_sign_bit_idx,
                                                                          one_input_var_idx);
                        } else {
                            this->create_signed_shift_overflow_protection(select_input_var_idx, zero_input_sign_bit_idx,
                                                                          zero_input_var_idx);
                        }
                    } else {
                        if (connect_zero_const) {
                            this->create_1x1_equivalence(zero_input_var_idx, mux_output_var_idx);
                            this->create_1x1_negated_implication(zero_input_var_idx, select_input_var_idx);
                            this->create_1x1_negated_implication(mux_output_var_idx, select_input_var_idx);
                        } else {
                            this->create_2x1_mux_shift_disallowed(zero_input_var_idx, one_input_var_idx,
                                                                  select_input_var_idx, mux_output_var_idx);
                        }
                    }
                } else {
                    if (connect_zero_const) {
                        this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx, mux_output_var_idx);
                    } else {
                        this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx,
                                             mux_output_var_idx);
                    }
                }
            }
            // sign bits before and after shifting must be identical if calculating in 2's complement
            if (this->calc_twos_complement) {
                int shift_input_sign_bit_idx;
                int shift_output_sign_bit_idx = this->shift_output_variables.at({idx, this->word_size - 1, v});
                if (idx == 1) {
                    shift_input_sign_bit_idx = this->output_value_variables.at({0, this->word_size - 1, v});
                } else {
                    shift_input_sign_bit_idx = this->input_select_mux_output_variables.at(
                            {idx, mcm::left, this->word_size - 1, v});
                }
                this->create_1x1_equivalence(shift_input_sign_bit_idx, shift_output_sign_bit_idx);
            }
        }
    }
}

void mcm::create_post_adder_shift_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
    for(int v=0; v<c_row_size(); v++) {
        for (auto stage = 0; stage < this->shift_word_size; stage++) {
            auto shift_width = (1 << stage);
            auto select_input_var_idx = this->input_post_adder_shift_value_variables.at({idx, stage});
            auto last_disallowed_shift_bit = shift_width - 1;
            for (auto w = 0; w < this->word_size; w++) {
                auto w_prev = w + shift_width;
                auto connect_zero_const = w_prev >= this->word_size;
                int zero_input_var_idx;
                int zero_input_sign_bit_idx;
                int one_input_var_idx;
                auto mux_output_var_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage, w, v});
                if (stage == 0) {
                    // connect shifter inputs
                    // shifter input is the adder output
                    zero_input_var_idx = this->adder_output_value_variables.at({idx, w, v});
                    zero_input_sign_bit_idx = this->adder_output_value_variables.at({idx, this->word_size - 1, v});
                    if (!connect_zero_const) {
                        one_input_var_idx = this->adder_output_value_variables.at({idx, w_prev, v});
                    }
                } else {
                    // connect output of previous stage
                    zero_input_var_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage - 1, w, v});
                    zero_input_sign_bit_idx = this->post_adder_shift_internal_mux_output_variables.at(
                            {idx, stage - 1, this->word_size - 1, v});
                    if (!connect_zero_const) {
                        one_input_var_idx = this->post_adder_shift_internal_mux_output_variables.at(
                                {idx, stage - 1, w_prev, v});
                    }
                }
                if (w <= last_disallowed_shift_bit) {
                    // shifting out 1s is not allowed in these places
                    if (connect_zero_const) {
                        if (this->calc_twos_complement) {
                            // connect the sign bit instead of a constant zero
                            this->create_2x1_mux_shift_disallowed(zero_input_var_idx, zero_input_sign_bit_idx,
                                                                  select_input_var_idx, mux_output_var_idx);
                        } else {
                            this->create_1x1_equivalence(zero_input_var_idx, mux_output_var_idx);
                            this->create_1x1_negated_implication(zero_input_var_idx, select_input_var_idx);
                            this->create_1x1_negated_implication(mux_output_var_idx, select_input_var_idx);
                        }
                    } else {
                        this->create_2x1_mux_shift_disallowed(zero_input_var_idx, one_input_var_idx,
                                                              select_input_var_idx, mux_output_var_idx);
                    }
                } else {
                    // we can shift bits around however we like
                    if (connect_zero_const) {
                        if (this->calc_twos_complement) {
                            // connect the sign bit instead of a constant zero
                            this->create_2x1_mux(zero_input_var_idx, zero_input_sign_bit_idx, select_input_var_idx,
                                                 mux_output_var_idx);
                        } else {
                            this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx,
                                                            mux_output_var_idx);
                        }
                    } else {
                        this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx,
                                             mux_output_var_idx);
                    }
                }
            }
            // sign bits before and after shifting must be identical if calculating in 2's complement
            if (this->calc_twos_complement) {
                int shift_output_sign_bit_idx = this->post_adder_shift_output_variables.at({idx, this->word_size - 1, v});
                int shift_input_sign_bit_idx = this->adder_output_value_variables.at({idx, this->word_size - 1, v});
                this->create_1x1_equivalence(shift_input_sign_bit_idx, shift_output_sign_bit_idx);
            }
        }
    }
}

void mcm::create_negate_select_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	auto select_var_idx = this->input_negate_select_variables.at(idx);
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            auto left_input_var_idx = this->shift_output_variables.at({idx, w, v});
            int right_input_var_idx;
            if (idx == 1){
                // right input is the output of the input node with idx = 0
                right_input_var_idx = this->output_value_variables.at({0, w, v});
            } else {
                // right input is the output of the right input select mux
                right_input_var_idx = this->input_select_mux_output_variables.at({idx, mcm::right, w, v});
            }
            for (auto &dir : this->input_directions) {
                auto mux_output_var_idx = this->negate_select_output_variables.at({idx, dir, w, v});
                if (dir == mcm::left) {
                    this->create_2x1_mux(right_input_var_idx, left_input_var_idx, select_var_idx, mux_output_var_idx);
                } else {
                    this->create_2x1_mux(left_input_var_idx, right_input_var_idx, select_var_idx, mux_output_var_idx);
                }
            }
        }
    }
}

void mcm::create_xor_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	auto negate_var_idx = this->input_negate_value_variables.at(idx);
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            auto input_var_idx = this->negate_select_output_variables.at({idx, mcm::right, w, v});
            auto output_var_idx = this->xor_output_variables.at({idx, w, v});
            this->create_2x1_xor(negate_var_idx, input_var_idx, output_var_idx);
        }
    }
}

void mcm::create_adder_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
    for(int v=0; v<c_row_size(); v++) {
        for (int w = 0; w < this->word_size; w++) {
            int c_i;
            if (w == 0) {
                // carry input = input negate value
                c_i = this->input_negate_value_variables.at(idx);
            } else {
                // carry input = carry output of last stage
                c_i = this->adder_carry_variables.at({idx, w - 1, v});
            }
            // in/out variables
            int a = this->negate_select_output_variables.at({idx, mcm::left, w, v});
            int b = this->xor_output_variables.at({idx, w, v});
            int s = this->adder_output_value_variables.at({idx, w, v});
            int c_o = this->adder_carry_variables.at({idx, w, v});
#if FPGA_ADD
            int xor_int = this->adder_XOR_internal_variables.at({idx, w});
            // build first XOR
            this->create_2x1_xor(a, b, xor_int);
            // build second XOR
            this->create_2x1_xor(xor_int, c_i, s);
            // build MUX
            this->create_2x1_mux(a, c_i, xor_int, c_o);
#else
            // build sum
            this->create_add_sum(a, b, c_i, s);
            // build carry
            this->create_add_carry(a, b, c_i, c_o);
            // build redundant clauses to increase strength of unit propagation
            // note (nfiege): this doesn't bring any speedup
            //this->create_add_redundant(a, b, c_i, s, c_o);
#endif
        }
        // disallow overflows
        if (this->calc_twos_complement) {
            this->create_signed_add_overflow_protection(this->input_negate_value_variables.at(idx),
                                                        this->negate_select_output_variables.at(
                                                                {idx, mcm::left, this->word_size - 1, v}),
                                                        this->negate_select_output_variables.at(
                                                                {idx, mcm::right, this->word_size - 1, v}),
                                                        this->output_value_variables.at({idx, this->word_size - 1, v}));
        } else {
            this->create_1x1_equivalence(this->adder_carry_variables.at({idx, this->word_size - 1, v}),
                                         this->input_negate_value_variables.at(idx));
        }
    }
}

void mcm::create_input_select_limitation_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
    if (idx <= 1) return;
    if (idx <= idx_input_buffer()) return;
    auto select_input_word_size = this->ceil_log2(idx);
	int max_representable_input_select = (1 << select_input_word_size) - 1;
    for (auto &dir : this->input_directions) {
        std::vector<int> x(select_input_word_size);
        for (int w = 0; w < select_input_word_size; w++) {
            x[w] = this->input_select_selection_variables.at({idx, dir, w});
        }
		for (int forbidden_number = max_representable_input_select; forbidden_number >= idx; forbidden_number--) {
			this->forbid_number(x, forbidden_number);
        }
	}
}

void mcm::create_shift_limitation_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	int max_representable_shift = (1 << this->shift_word_size) - 1;
	std::vector<int> x(this->shift_word_size);
	for (int w = 0; w < this->shift_word_size; w++) {
		x[w] = this->input_shift_value_variables.at({idx, w});
    }
	for (int forbidden_number = max_representable_shift; forbidden_number > this->max_shift; forbidden_number--) {
		this->forbid_number(x, forbidden_number);
	}
}

void mcm::create_post_adder_shift_limitation_constraints(int idx, formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	int max_representable_shift = (1 << this->shift_word_size) - 1;
	std::vector<int> x(this->shift_word_size);
    for (int w = 0; w < this->shift_word_size; w++) {
		x[w] = this->input_post_adder_shift_value_variables.at({idx, w});
    }
	for (int forbidden_number = max_representable_shift; forbidden_number > this->max_shift; forbidden_number--) {
		this->forbid_number(x, forbidden_number);
	}
}

void mcm::get_solution_from_backend() {
	// clear containers
	this->input_select.clear();
	this->input_select_mux_output.clear();
	this->shift_value.clear();
	this->negate_select.clear();
	this->subtract.clear();
	this->post_adder_shift_value.clear();
	this->add_result_values.clear();
	this->output_values.clear();
	this->coeff_word_size_values.clear();
	this->can_cut_msb_values.clear();
	this->coeff_word_size_sum_values.clear();
	this->shift_sum_values.clear();
	this->num_FAs_value = 0;
	// get solution
    for(int idx = 0; idx <= (this->num_adders + idx_input_buffer()); idx++) {
        for (int v=0; v<c_row_size(); v++) {
            // output_values
            this->output_values[{idx,v}] = 0;
            for (int w = 0; w < this->word_size; w++) {
                this->output_values[{idx,v}] += (this->get_result_value(this->output_value_variables.at({idx, w, v})) << w);
            }
            if (this->calc_twos_complement)
                this->output_values[{idx,v}] = sign_extend(this->output_values[{idx,v}], this->word_size);
            if (idx > (0 + idx_input_buffer())) {
                if (idx > 1 ) {
                    // input_select
                    for (auto &dir : this->input_directions) {
                        auto input_select_width = this->ceil_log2(idx);
                        this->input_select[{idx, dir}] = 0;
                        for (auto w = 0; w < input_select_width; w++) {
                            this->input_select[{idx, dir}] += (
                                    this->get_result_value(this->input_select_selection_variables[{idx, dir, w}]) << w);
                        }
                    }
                } else {
                    this->input_select[{idx, input_direction::left}] = 0;
                    this->input_select[{idx, input_direction::right}] = 0;
                }
                // shift_value
                this->shift_value[idx] = 0;
                for (auto w = 0; w < this->shift_word_size; w++) {
                    this->shift_value[idx] += (this->get_result_value(this->input_shift_value_variables[{idx, w}])
                            << w);
                }
                // negate_select
                this->negate_select[idx] = this->get_result_value(this->input_negate_select_variables[idx]);
                // subtract
                this->subtract[idx] = this->get_result_value(this->input_negate_value_variables[idx]);
                // output shift
                if (this->enable_node_output_shift) {
                    this->post_adder_shift_value[idx] = 0;
                    for (auto w = 0; w < this->shift_word_size; w++) {
                        this->post_adder_shift_value[idx] += (
                                this->get_result_value(this->input_post_adder_shift_value_variables[{idx, w}]) << w);
                    }
                }
                // add result
                this->add_result_values[{idx,v}] = 0;
                for (auto w = 0; w < this->word_size; w++) {
                    this->add_result_values[{idx,v}] += (
                            this->get_result_value(this->adder_output_value_variables[{idx, w, v}]) << w);
                }
                if (this->calc_twos_complement)
                    this->add_result_values[{idx,v}] = sign_extend(this->add_result_values[{idx,v}], this->word_size);
                if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
                    // coeff word size internal
                    for (auto w = this->word_size - 1; w >= 0; w--) {
                        auto max_val = this->word_size - w;
                        auto max_val_w = this->ceil_log2(max_val + 1);
                        auto internal_val = 0;
                        for (auto x = 0; x < max_val_w; x++) {
                            auto bit_val = this->get_result_value(
                                    this->full_adder_coeff_word_size_internal_variables.at({idx, w, x}));
                            internal_val += (bit_val << x);
                        }
                    }
                    // coeff word size
                    this->coeff_word_size_values[idx] = 0;
                    auto num_bits_word_size = this->ceil_log2(this->word_size + 1);
                    for (auto w = 0; w < num_bits_word_size; w++) {
                        this->coeff_word_size_values[idx] += (
                                this->get_result_value(this->full_adder_coeff_word_size_variables.at({idx, w})) << w);
                    }
                    // cut msb
                    this->can_cut_msb_values[idx] = this->get_result_value(this->full_adder_msb_variables.at(idx));
                    // coeff word size sum
                    this->coeff_word_size_sum_values[idx] = 0;
                    auto coeff_sum_output_word_size = this->ceil_log2((idx * this->word_size) + 1);
                    for (auto w = 0; w < coeff_sum_output_word_size; w++) {
                        this->coeff_word_size_sum_values[idx] += (
                                this->get_result_value(this->full_adder_word_size_sum_variables.at({idx, w})) << w);
                    }
                    // shift gain
                    this->shift_gain_values[idx] = 0;
                    for (auto w = 0; w < this->shift_word_size; w++) {
                        this->shift_gain_values[idx] += (
                                this->get_result_value(this->full_adder_shift_gain_variables.at({idx, w})) << w);
                    }
                    // shift sum
                    this->shift_sum_values[idx] = 0;
                    auto shift_sum_output_word_size = this->ceil_log2((idx * this->max_shift) + 1);
                    for (auto w = 0; w < shift_sum_output_word_size; w++) {
                        this->shift_sum_values[idx] += (
                                this->get_result_value(this->full_adder_shift_sum_variables.at({idx, w})) << w);
                    }
                }
            }
        }
        if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
            // number of FAs
            this->num_FAs_value = 0;
            auto input_word_size_add = this->ceil_log2(this->word_size * this->num_adders + 1);
            auto input_word_size_sub = this->ceil_log2((this->num_adders + 1) * this->max_shift + 1);
            auto output_word_size = std::max(input_word_size_add, input_word_size_sub) + 1;
            for (auto w = 0; w < output_word_size; w++) {
                this->num_FAs_value += (this->get_result_value(this->full_adder_result_variables.at(w)) << w);
            }
        }
    }
}

void mcm::print_solution() {
	if (this->found_solution) {
		std::cout << "Solution for Vector" << std::endl;
		for (auto &v : this->C) {
            std::cout << "  V = <";
		    for (auto c : v){
                std::cout << " " << c;
		    }
			std::cout << " >" <<std::endl;
		}
		//ToDo make it look good for vector inputs
        for (int v=0; v<c_row_size(); v++) {
            std::cout << "#adders = " << this->num_adders << ", word size = " << this->word_size << std::endl;

            //print input nodes
            for(int i = 0; i <= idx_input_buffer(); i++){
                std::cout << "  node #" << i << " = "
                          << (this->calc_twos_complement ? sign_extend(this->output_values[{i, v}], this->word_size)
                                                         : this->output_values[{i, v}]) << std::endl;
            }

            //print following nodes
            for (auto idx = (1 + idx_input_buffer()); idx <= (this->num_adders + idx_input_buffer()); idx++) {
                std::cout << "  node #" << idx << " = "
                          << (this->calc_twos_complement ? sign_extend((int64_t) this->output_values[{idx, v}],
                                                                       this->word_size) : this->output_values[{idx, v}])
                          << std::endl;
                std::cout << "    left input: node " << this->input_select[{idx, mcm::left}] << std::endl;
                std::cout << "    right input: node " << this->input_select[{idx, mcm::right}] << std::endl;
                std::cout << "    shift value: " << this->shift_value[idx] << std::endl;
                std::cout << "    negate select: " << this->negate_select[idx]
                          << (this->negate_select[idx] == 1 ? " (non-shifted)" : " (shifted)") << std::endl;
                std::cout << "    subtract: " << this->subtract[idx] << std::endl;
                if (this->enable_node_output_shift) {
                    std::cout << "    post adder right shift value: " << this->post_adder_shift_value[idx] << std::endl;
                }
            }
        }
		std::cerr << "Adder graph: " << this->get_adder_graph_description() << std::endl;
	}
	else {
        std::cout << "Failed to find solution for Vector" << std::endl;
	    for (auto &v : this->C) {
	        std::cout << "  V = <";
	        for (auto c : v){
	            std::cout << " " << c;
	        }
	        std::cout << " >" <<std::endl;
		}
	}
}

bool mcm::solution_is_valid() {
	bool valid = true;
    for(int idx = idx_input_buffer() + 1; idx <= (this->num_adders + idx_input_buffer()); idx++) {
        for (int v=0; v<c_row_size(); v++) {

            std::cout << "test where I left the function 1" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

            // verify node inputs
            int64_t input_node_idx_l = 0;
            int64_t input_node_idx_r = 0;

            int64_t actual_input_value_l = 1;
            int64_t actual_input_value_r = 1;

            if (idx > 1 ) {
                for (auto &dir : this->input_directions) {
                    for (auto w = 0; w < this->word_size; w++) {
                        this->input_select_mux_output[{idx, dir, v}] += (
                                this->get_result_value(this->input_select_mux_output_variables[{idx, dir, w, v}]) << w);
                    }
                }
                if (this->calc_twos_complement)
                    this->input_select_mux_output[{idx, mcm::left, v}] = sign_extend(
                            this->input_select_mux_output[{idx, mcm::left, v}], this->word_size);
                if (this->calc_twos_complement)
                    this->input_select_mux_output[{idx, mcm::right, v}] = sign_extend(
                            this->input_select_mux_output[{idx, mcm::right, v}], this->word_size);
                input_node_idx_l = this->input_select[{idx, mcm::left}];
                input_node_idx_r = this->input_select[{idx, mcm::right}];
                actual_input_value_l = this->input_select_mux_output[{idx, mcm::left, v}];
                actual_input_value_r = this->input_select_mux_output[{idx, mcm::right, v}];
            } else {
                this->input_select_mux_output[{idx, mcm::left, v}] = this->input_select_mux_output[{idx, mcm::right, v}] = 1;
            }
            std::cout << "test where I left the function 2" << std::endl;
            std::cout << "valid value: " << valid << std::endl;
            std::cout << "left_input_value: " << this->output_values[{input_node_idx_l,v}] << std::endl;
            std::cout << "right_input_value: " << this->output_values[{input_node_idx_r,v}] << std::endl;

            int64_t left_input_value = this->output_values[{input_node_idx_l,v}];
            int64_t right_input_value = this->output_values[{input_node_idx_r,v}];
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " left input" << std::endl;
                std::cout << "  input select = " << input_node_idx_l << std::endl;
                std::cout << "  value = " << actual_input_value_l << std::endl;
            }
            if (left_input_value != actual_input_value_l) {
                std::cout << "node #" << idx << " has invalid left input" << std::endl;
                std::cout << "  input select = " << input_node_idx_l << std::endl;
                std::cout << "  expected value " << left_input_value << " but got " << actual_input_value_l
                          << std::endl;
                auto num_muxs = (1 << this->ceil_log2(idx)) - 1;
                for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
                    int64_t mux_output = 0;
                    for (auto w = 0; w < this->word_size; w++) {
                        mux_output += (this->get_result_value(
                                this->input_select_mux_variables[{idx, mcm::left, mux_idx, w, v}]) << w);
                    }
                    std::cout << "    mux #" << mux_idx << " output: " << mux_output << std::endl;
                }
                valid = false;
            }
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " right input" << std::endl;
                std::cout << "  input select = " << input_node_idx_r << std::endl;
                std::cout << "  value = " << actual_input_value_r << std::endl;
            }
            if (right_input_value != actual_input_value_r) {
                std::cout << "node #" << idx << " has invalid right input" << std::endl;
                std::cout << "  input select = " << input_node_idx_r << std::endl;
                std::cout << "  expected value " << right_input_value << " but got " << actual_input_value_r
                          << std::endl;
                int64_t num_muxs = (1 << this->ceil_log2(idx)) - 1;
                for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
                    int64_t mux_output = 0;
                    for (auto w = 0; w < this->word_size; w++) {
                        mux_output += (this->get_result_value(
                                this->input_select_mux_variables[{idx, mcm::right, mux_idx, w, v}]) << w);
                    }
                    std::cout << "    mux #" << mux_idx << " output: " << mux_output << std::endl;
                }
                valid = false;
            }

            std::cout << "test where I left the function 3" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

            // verify shifter output
            int64_t expected_shift_output = (((int64_t) left_input_value)
                    << this->shift_value[idx]);// % (int64_t)(1 << this->word_size);
            if (this->calc_twos_complement) expected_shift_output = sign_extend(expected_shift_output, this->word_size);
            int64_t actual_shift_output = 0;
            for (int w = 0; w < this->word_size; w++) {
                actual_shift_output += (this->get_result_value(this->shift_output_variables[{idx, w, v}]) << w);
                std::cout << "shift_output_variables[{idx, w, v}: " << shift_output_variables[{idx, w, v}] << std::endl;
                std::cout << "this->get_result_value(shift_output_variables[{idx, w, v}]): " << (this->get_result_value(this->shift_output_variables[{idx, w, v}]) << w)<< std::endl;
            }
            if (this->calc_twos_complement) actual_shift_output = sign_extend(actual_shift_output, this->word_size);
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " shift output" << std::endl;
                std::cout << "  input value = " << left_input_value << std::endl;
                std::cout << "  shift value = " << this->shift_value[idx] << std::endl;
                std::cout << "  output value = " << actual_shift_output << std::endl;
            }
            if (expected_shift_output != actual_shift_output) {
                std::cout << "node #" << idx << " has invalid shift output" << std::endl;
                std::cout << "  input value = " << left_input_value << std::endl;
                std::cout << "  shift value = " << this->shift_value[idx] << std::endl;
                std::cout << "  expected output value = " << expected_shift_output << std::endl;
                std::cout << "  actual output value = " << actual_shift_output << std::endl;
                valid = false;
            }
            // verify negate mux outputs
            int64_t negate_mux_output_l = actual_shift_output;
            int64_t negate_mux_output_r = right_input_value;
            if (this->get_result_value(this->input_negate_select_variables[idx]) == 0) {
                negate_mux_output_l = right_input_value;
                negate_mux_output_r = actual_shift_output;
            }
            std::map<mcm::input_direction, int> actual_negate_mux_output;
            for (auto &dir : this->input_directions) {
                for (auto w = 0; w < this->word_size; w++) {
                    actual_negate_mux_output[dir] += (
                            this->get_result_value(this->negate_select_output_variables[{idx, dir, w, v}]) << w);
                }
            }
            if (this->calc_twos_complement)
                actual_negate_mux_output[mcm::left] = sign_extend(actual_negate_mux_output[mcm::left], this->word_size);
            if (this->calc_twos_complement)
                actual_negate_mux_output[mcm::right] = sign_extend(actual_negate_mux_output[mcm::right],
                                                                   this->word_size);
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " left negate select mux output" << std::endl;
                std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx])
                          << std::endl;
                std::cout << "  output value = " << actual_negate_mux_output[mcm::left] << std::endl;
            }
            if (negate_mux_output_l != actual_negate_mux_output[mcm::left]) {
                std::cout << "node #" << idx << " has invalid left negate select mux output" << std::endl;
                std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx])
                          << std::endl;
                std::cout << "  actual value = " << actual_negate_mux_output[mcm::left] << std::endl;
                std::cout << "  expected value = " << negate_mux_output_l << std::endl;
                valid = false;
            }

            std::cout << "test where I left the function 4" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " right negate select mux output" << std::endl;
                std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx])
                          << std::endl;
                std::cout << "  output value = " << actual_negate_mux_output[mcm::right] << std::endl;
            }
            if (negate_mux_output_r != actual_negate_mux_output[mcm::right]) {
                std::cout << "node #" << idx << " has invalid right negate select mux output" << std::endl;
                std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx])
                          << std::endl;
                std::cout << "  actual value = " << actual_negate_mux_output[mcm::right] << std::endl;
                std::cout << "  expected value = " << negate_mux_output_r << std::endl;
                valid = false;
            }
            // verify xor output
            int64_t sub = this->get_result_value(this->input_negate_value_variables[idx]);
            int64_t expected_xor_output =
                    sub == 1 ? (~negate_mux_output_r) & ((((int64_t) 1) << this->word_size) - 1) : negate_mux_output_r;
            if (this->calc_twos_complement) expected_xor_output = sign_extend(expected_xor_output, this->word_size);
            int64_t actual_xor_output = 0;
            for (int w = 0; w < this->word_size; w++) {
                actual_xor_output += (this->get_result_value(this->xor_output_variables[{idx, w, v}]) << w);
            }
            if (this->calc_twos_complement) actual_xor_output = sign_extend(actual_xor_output, this->word_size);
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " xor output" << std::endl;
                std::cout << "  sub = " << sub << std::endl;
                std::cout << "  input value = " << negate_mux_output_r << std::endl;
                std::cout << "  output value = " << actual_xor_output << std::endl;
            }
            if (expected_xor_output != actual_xor_output) {
                std::cout << "node #" << idx << " has invalid xor output" << std::endl;
                std::cout << "  sub = " << sub << std::endl;
                std::cout << "  input value = " << negate_mux_output_r << std::endl;
                std::cout << "  actual output value = " << actual_xor_output << std::endl;
                std::cout << "  expected output value = " << expected_xor_output << std::endl;
                valid = false;
            }
            // verify adder output
            int64_t expected_adder_output = (sub == 1) ? (negate_mux_output_l - negate_mux_output_r) : (
                    negate_mux_output_l + negate_mux_output_r);
            if (this->calc_twos_complement) expected_adder_output = sign_extend(expected_adder_output, this->word_size);
            int64_t actual_adder_output = 0;
            for (int w = 0; w < this->word_size; w++) {
                actual_adder_output += (this->get_result_value(this->adder_output_value_variables[{idx, w, v}]) << w);
            }
            if (this->calc_twos_complement) actual_adder_output = sign_extend(actual_adder_output, this->word_size);
            if (this->verbosity == verbosity_mode::debug_mode) {
                std::cout << "node #" << idx << " adder output value" << std::endl;
                std::cout << "  sub = " << sub << std::endl;
                std::cout << "  left input value = " << negate_mux_output_l << std::endl;
                std::cout << "  right input value = " << actual_xor_output << std::endl;
                std::cout << "  actual output value = " << actual_adder_output << std::endl;
            }
            if (expected_adder_output != actual_adder_output) {
                std::cout << "node #" << idx << " has invalid adder output value" << std::endl;
                std::cout << "  sub = " << sub << std::endl;
                std::cout << "  left input value = " << negate_mux_output_l << std::endl;
                std::cout << "  right input value = " << actual_xor_output << std::endl;
                std::cout << "  expected output value = " << expected_adder_output << std::endl;
                std::cout << "  actual output value = " << actual_adder_output << std::endl;
                valid = false;
            }

            std::cout << "test where I left the function 5" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

            if (this->enable_node_output_shift) {
                // verify post adder shift output
                int64_t expected_post_adder_shift_output = actual_adder_output >> this->post_adder_shift_value.at(idx);
                if (this->calc_twos_complement)
                    expected_post_adder_shift_output = sign_extend(expected_post_adder_shift_output, this->word_size);
                int64_t actual_post_adder_shift_output = 0;
                for (int w = 0; w < this->word_size; w++) {
                    actual_post_adder_shift_output += (
                            this->get_result_value(this->post_adder_shift_output_variables[{idx, w, v}]) << w);
                }
                if (this->calc_twos_complement)
                    actual_post_adder_shift_output = sign_extend(actual_post_adder_shift_output, this->word_size);
                if (this->verbosity == verbosity_mode::debug_mode) {
                    std::cout << "node #" << idx << " post adder shift output value" << std::endl;
                    std::cout << "  shift value = " << this->post_adder_shift_value.at(idx) << std::endl;
                    std::cout << "  input value = " << actual_adder_output << std::endl;
                    std::cout << "  actual output value = " << actual_post_adder_shift_output << std::endl;
                }
                if (expected_post_adder_shift_output != actual_post_adder_shift_output) {
                    std::cout << "node #" << idx << " has invalid post adder shift output value" << std::endl;
                    std::cout << "  shift value = " << this->post_adder_shift_value.at(idx) << std::endl;
                    std::cout << "  input value = " << actual_adder_output << std::endl;
                    std::cout << "  expected output value = " << expected_post_adder_shift_output << std::endl;
                    std::cout << "  actual output value = " << actual_post_adder_shift_output << std::endl;
                    valid = false;
                }
            }

            std::cout << "test where I left the function 6" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

            if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
                // coeff word size
                auto add_result = this->add_result_values.at({idx,v});
                auto expected_add_result_word_size = this->ceil_log2(std::abs(add_result) + 1);
                auto actual_add_result_word_size = this->coeff_word_size_values.at(idx);
                if (this->verbosity == verbosity_mode::debug_mode) {
                    std::cout << "node #" << idx << " add result word size" << std::endl;
                    std::cout << "  add result = " << add_result << std::endl;
                    std::cout << "  actual word size = " << actual_add_result_word_size << std::endl;
                }
                if (expected_add_result_word_size != actual_add_result_word_size) {
                    std::cout << "node #" << idx << " has invalid add result word size" << std::endl;
                    std::cout << "  add result = " << add_result << std::endl;
                    std::cout << "  expected word size = " << expected_add_result_word_size << std::endl;
                    std::cout << "  actual word size = " << actual_add_result_word_size << std::endl;
                    valid = false;
                }
                // coeff word size sum
                auto expected_sum = 0;
                for (auto prev_idx = 1; prev_idx <= idx; prev_idx++) {
                    expected_sum += this->coeff_word_size_values.at(prev_idx);
                }
                auto actual_sum = this->coeff_word_size_sum_values.at(idx);
                if (this->verbosity == verbosity_mode::debug_mode) {
                    std::cout << "node #" << idx << " add result word size sum" << std::endl;
                    for (auto prev_idx = 1; prev_idx <= idx; prev_idx++) {
                        std::cout << "  value " << prev_idx << " = " << this->coeff_word_size_values.at(prev_idx)
                                  << std::endl;
                    }
                    std::cout << "  actual sum = " << actual_sum << std::endl;
                }
                if (expected_sum != actual_sum) {
                    std::cout << "node #" << idx << " has invalid add result word size sum" << std::endl;
                    for (auto prev_idx = 1; prev_idx <= idx; prev_idx++) {
                        std::cout << "  value " << prev_idx << " = " << this->coeff_word_size_values.at(prev_idx)
                                  << std::endl;
                    }
                    std::cout << "  expected sum = " << expected_sum << std::endl;
                    std::cout << "  actual sum = " << actual_sum << std::endl;
                    valid = false;
                }
            }
            std::cout << "test where I left the function 7" << std::endl;
            std::cout << "valid value: " << valid << std::endl;

        }
    }
	return valid;
	//return true;
}

void mcm::create_cnf_file() {
	std::ofstream f;
	std::stringstream constants;
	for (int i = 0; i < this->c_column_size(); i++) {
	    if (i != 0) constants << "_";
	    constants << this->C[i][0];
	    //TODO note Christoph: take a look later just naming of the .cnf file
	}
	std::string filename;
	if (this->max_full_adders != FULL_ADDERS_UNLIMITED) {
		filename = constants.str() + "-" + std::to_string(this->num_adders) + "-" + std::to_string(this->max_full_adders) + ".cnf";
	}
	else {
		filename = constants.str() + "-" + std::to_string(this->num_adders) + ".cnf";
	}
	f.open(filename.c_str());
	f << "p cnf " << this->variable_counter << " " << this->constraint_counter << std::endl;
	f << this->cnf_clauses.str();
	f.close();
}
void mcm::create_mcm_input_constraints(mcm::formulation_mode mode) {
    if (mode != formulation_mode::reset_all) return;

    std::vector<int> input_bits(this->word_size);

    for(int i=0; i<=c_row_size()-1; i++){
        for (int j=0; j<=c_row_size()-1; j++) {
            for (auto w = 0; w < this->word_size; w++) {
                input_bits[w] = this->output_value_variables.at({i, w, j});
                std::cout<<"input bits: "<< input_bits[w] <<std::endl;
            }
            if(i==j) {
                force_number(input_bits, 1);
            }else{
                force_number(input_bits, 0);
            }
        }
    }
}
void mcm::create_mcm_output_constraints(formulation_mode mode) {
	if (mode != formulation_mode::reset_all) return;
	for(int m=1; m<=c_column_size(); m++){
	    std::vector<int> or_me;
	    for (int idx = idx_input_buffer() + 1; idx <= (this->num_adders + idx_input_buffer()); idx++) {
            std::cout<<"mcm_output_variables[{idx, m}: "<< this->mcm_output_variables[{idx, m}] <<std::endl;
	        or_me.emplace_back(this->mcm_output_variables[{idx, m}]);
	        for(int v=0; v<c_row_size(); v++) {
	            for (int w = 0; w < this->word_size; w++) {
	                if (((C[m-1][v] >> w) & 1) == 1) {
	                    this->create_1x1_implication(this->mcm_output_variables[{idx, m}],
                                                  this->output_value_variables[{idx, w, v}]);
	                } else {
	                    this->create_1x1_negated_implication(this->mcm_output_variables[{idx, m}],
                                                          this->output_value_variables[{idx, w, v}]);
	                }
	            }
	        }
	    }
	    if (this->calc_twos_complement and this->sign_inversion_allowed[m]) {
	        // also allow the solver to choose -c instead of c if it's easier to implement
	        for (int idx = idx_input_buffer() + 1; idx <= (this->num_adders + idx_input_buffer()); idx++) {
                std::cout<<"mcm_output_variables[{idx, m}: "<< this->mcm_output_variables[{idx, -m}] <<std::endl;
	            or_me.emplace_back(this->mcm_output_variables[{idx, -m}]);
	            for(int v=0; v<c_row_size(); v++) {
	                for (int w = 0; w < this->word_size; w++) {
	                    if ((((-C[m-1][v]) >> w) & 1) == 1) {
	                        this->create_1x1_implication(this->mcm_output_variables[{idx, -m}],
                                                      this->output_value_variables[{idx, w, v}]);
	                    } else {
	                        this->create_1x1_negated_implication(this->mcm_output_variables[{idx, -m}],
                                                              this->output_value_variables[{idx, w, v}]);
	                    }
	                }
	            }
	        }
	    }
	    this->create_or(or_me);
    }
}

void mcm::create_mcm_output_variables(int idx) {
    for(int m=1; m<=c_column_size(); m++) {
        this->mcm_output_variables[{idx, m}] = ++this->variable_counter;
        this->create_new_variable(this->variable_counter);
        if (this->calc_twos_complement and this->sign_inversion_allowed[m]) {
            this->mcm_output_variables[{idx, -m}] = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
        }
    }
}

void mcm::create_odd_fundamentals_constraints(int idx, formulation_mode mode) {
    if(idx_input_buffer() > 0) return;
    //this constraint is not needed for CMM and SOP
	if (mode != formulation_mode::reset_all) return;
    for(int v=0; v<c_row_size(); v++) {
        this->force_bit(this->output_value_variables.at({idx, 0, v}), 1);
    }
}

int64_t mcm::sign_extend(int64_t x, int w) {
	auto sign_bit = (x >> (w-1)) & 1;
	if (sign_bit == 0) return x; // x >= 0 -> no conversion needed
	auto mask = (1 << w) - 1;
	mask = ~mask;
	x = x | mask;
	return x;
}

std::string mcm::get_adder_graph_description() {
	std::stringstream s;
	if (!this->found_solution) return s.str();
	s << "{";
	std::map<int, int> stage;
	//build initial stage for inputs
	for(int i=0; i <= idx_input_buffer(); i++){
        stage[i] = i;
	}

	for (int idx = (1 + idx_input_buffer()); idx <= (this->num_adders + idx_input_buffer()); idx++) {
        // get left and right inputs and their shift
        int left_idx;
        int right_idx;
        int left_input;
        int right_input;
        int left_shift;
        int right_shift;
        int left_stage;
        int right_stage;
        int current_stage;
        // for (int v=0; v<c_row_size(); v++) {
        //     // insert comma if needed
        //     if (idx > (1 + idx_input_buffer())) s << ",";
        //     left_idx = this->input_select.at({idx, mcm::left});
        //     right_idx = this->input_select.at({idx, mcm::right});
        //     left_input = this->output_values.at({left_idx,v});
        //     right_input = this->output_values.at({right_idx,v});
        //     left_shift = this->shift_value.at(idx);
        //     right_shift = 0;
        //     // add/sub?
        //     if (this->negate_select.at(idx) == 0) {
        //         // swap again for subtract
        //         int idx_cpy = left_idx;
        //         int input_cpy = left_input;
        //         int shift_cpy = left_shift;
        //         left_idx = right_idx;
        //         left_input = right_input;
        //         left_shift = right_shift;
        //         right_idx = idx_cpy;
        //         right_input = input_cpy;
        //         right_shift = shift_cpy;
        //     }
        //     if (this->subtract.at(idx) == 1) {
        //         right_input *= -1;
        //     }
        //     // calc stage
        //     left_stage = stage.at(left_idx);
        //     right_stage = stage.at(right_idx);
        //     current_stage = std::max(left_stage, right_stage) + 1;
        //     stage[idx] = current_stage;
        // }

        //TODO there has to be a way to do this more efficient this is pretty ugly

        // basic node info
        s << "{'A',[" ;
        for (int v=0; v<c_row_size(); v++) {
            left_idx = this->input_select.at({idx, mcm::left});
            right_idx = this->input_select.at({idx, mcm::right});
            left_input = this->output_values.at({left_idx,v});
            right_input = this->output_values.at({right_idx,v});
            left_shift = this->shift_value.at(idx);
            right_shift = 0;
            // add/sub?
            if (this->negate_select.at(idx) == 0) {
                // swap again for subtract
                int idx_cpy = left_idx;
                int input_cpy = left_input;
                int shift_cpy = left_shift;
                left_idx = right_idx;
                left_input = right_input;
                left_shift = right_shift;
                right_idx = idx_cpy;
                right_input = input_cpy;
                right_shift = shift_cpy;
            }
            if (this->subtract.at(idx) == 1) {
                right_input *= -1;
            }
            // calc stage
            left_stage = stage.at(left_idx);
            right_stage = stage.at(right_idx);
            current_stage = std::max(left_stage, right_stage) + 1;
            stage[idx] = current_stage;

            s << this->output_values.at({idx,v});
            if (v!= c_row_size()-1) s << ",";
        }
        s << "]," << current_stage;
        if (this->enable_node_output_shift) {
            s << "," << this->post_adder_shift_value.at(idx);
        }

        // left input
        s << ",[";
        for (int v=0; v<c_row_size(); v++) {
            left_idx = this->input_select.at({idx, mcm::left});
            right_idx = this->input_select.at({idx, mcm::right});
            left_input = this->output_values.at({left_idx,v});
            right_input = this->output_values.at({right_idx,v});
            left_shift = this->shift_value.at(idx);
            right_shift = 0;
            // add/sub?
            if (this->negate_select.at(idx) == 0) {
                // swap again for subtract
                int idx_cpy = left_idx;
                int input_cpy = left_input;
                int shift_cpy = left_shift;
                left_idx = right_idx;
                left_input = right_input;
                left_shift = right_shift;
                right_idx = idx_cpy;
                right_input = input_cpy;
                right_shift = shift_cpy;
            }
            if (this->subtract.at(idx) == 1) {
                right_input *= -1;
            }
            // calc stage
            left_stage = stage.at(left_idx);
            right_stage = stage.at(right_idx);
            s << left_input;
            if (v != c_row_size()-1) s << ",";
        }
        s << "]," << left_stage << "," << left_shift;

        // right input
        s << ",[";
        for (int v=0; v<c_row_size(); v++) {
            left_idx = this->input_select.at({idx, mcm::left});
            right_idx = this->input_select.at({idx, mcm::right});
            left_input = this->output_values.at({left_idx,v});
            right_input = this->output_values.at({right_idx,v});
            left_shift = this->shift_value.at(idx);
            right_shift = 0;
            // add/sub?
            if (this->negate_select.at(idx) == 0) {
                // swap again for subtract
                int idx_cpy = left_idx;
                int input_cpy = left_input;
                int shift_cpy = left_shift;
                left_idx = right_idx;
                left_input = right_input;
                left_shift = right_shift;
                right_idx = idx_cpy;
                right_input = input_cpy;
                right_shift = shift_cpy;
            }
            if (this->subtract.at(idx) == 1) {
                right_input *= -1;
            }
            // calc stage
            left_stage = stage.at(left_idx);
            right_stage = stage.at(right_idx);
            s << right_input;
            if (v != c_row_size()-1) s << ",";
        }
        s << "]," << right_stage << "," << right_shift;
        // close bracket
        s << "}";

        // s <<"]";
        // // basic node info
        // s << "{'A',[" << this->output_values.at({idx,v}) << "]," << current_stage;
        // if (this->enable_node_output_shift) {
        //     s << "," << this->post_adder_shift_value.at(idx);
        // }
        // // left input
        // s << ",[" << left_input << "]," << left_stage << "," << left_shift;
        // // right input
        // s << ",[" << right_input << "]," << right_stage << "," << right_shift;
        // // close bracket
        // s << "}";
    }
	s << "}";
	return s.str();
}

void mcm::set_min_add(int new_min_add) {
	this->num_adders = std::max(this->num_adders, new_min_add-1);
	this->num_adders = std::max(this->num_adders, 0);
}

void mcm::also_minimize_full_adders() {
	this->minimize_full_adders = true;
}

void mcm::allow_node_output_shift() {
	this->enable_node_output_shift = true;
}

std::pair<int, int> mcm::solution_is_optimal() {
	return {this->num_add_opt, this->num_FA_opt};
}

void mcm::ignore_sign(bool only_apply_to_negative_coefficients) {
    for(int m=1; m<=c_column_size(); m++) {
        if (only_apply_to_negative_coefficients) {
            // only allow sign inversion for all negative coefficients
            if (this->inverted_coeff_requested[C[m - 1]]) {
                if (this->vector_all_positive(C[m - 1])) {
                    this->sign_inversion_allowed[m] = true;
                }
            }
        } else {
            // only the solver to choose sign for all coefficients
            this->sign_inversion_allowed[m] = true;
        }
	}
}

bool mcm::vector_all_positive(std::vector<int> v){
    bool all_positive = true;
    for(auto &c : v){
        //skip over negative values
        if(c>=0)continue;
        else{
            //break out for the first non negative value and return false
            all_positive = false;
            break;
        }
    }
    //only negative elements => return true
    return all_positive;
}

void mcm::create_full_adder_coeff_word_size_constraints(int idx, formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	std::vector<int> abs_coeff_bits(this->word_size);
    for(int v=0; v<c_row_size(); v++) {
        if (this->calc_twos_complement) {
            // compute abs(c) using a MUX and an inversion and a +1 adder
            /*int carry_bit = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
            this->force_bit(carry_bit, 1);*/
            int carry_bit = this->init_const_one_bit();
            std::vector<int> inv_c_bits(this->word_size);
            auto &sign_bit = this->adder_output_value_variables.at({idx, this->word_size - 1, v});
            for (int w = 0; w < this->word_size; w++) {
                // first, implement c*(-1)
                auto sum_bit = inv_c_bits[w] = ++this->variable_counter; // sum bit
                this->create_new_variable(this->variable_counter);
                auto &add_bit = this->adder_output_value_variables.at({idx, w, v});
                if (w < this->word_size - 1) {
                    int carry_out_bit = ++this->variable_counter;
                    this->create_new_variable(this->variable_counter);
                    // create clauses for sum and carry bits
                    //this->create_half_adder_inv_b(carry_bit, add_bit, sum_bit, carry_out_bit);
                    this->create_half_adder({carry_bit, false}, {add_bit, true}, {sum_bit, false},
                                            {carry_out_bit, false});
                    // pass carry bit to next stage
                    carry_bit = carry_out_bit;
                } else {
                    // only create clauses for sum bit
                    //this->create_half_adder_inv_b(carry_bit, add_bit, sum_bit);
                    this->create_half_adder({carry_bit, false}, {add_bit, true}, {sum_bit, false});
                }
                // now, implement the MUX, controlled by the sign bit
                auto mux_bit = abs_coeff_bits[w] = ++this->variable_counter;
                this->create_new_variable(this->variable_counter);
                this->create_2x1_mux(add_bit, sum_bit, sign_bit, mux_bit);
            }
        } else {
            // abs(c) = c because c is an unsigned number
            for (int w = 0; w < this->word_size; w++) {
                abs_coeff_bits[w] = this->adder_output_value_variables.at({idx, w, v});
            }
        }

        // compute word size of abs(c)
        // carry-input = 0
        /*int carry_bit = ++this->variable_counter;
        this->create_new_variable(this->variable_counter);
        this->force_bit(carry_bit, 0);*/
        int carry_bit = this->init_const_zero_bit();
        // initial value = 0
        std::vector<int> val(1);
        /*val[0] = ++this->variable_counter;
        this->create_new_variable(this->variable_counter);
        this->force_bit(val[0], 0);*/
        val[0] = this->init_const_zero_bit();
        for (int w = this->word_size - 1; w >= 0; w--) {
            this->full_adder_coeff_word_size_internal_carry_input_variables[{idx, w}] = carry_bit;
            auto w_in = val.size();
            auto max_val = this->word_size - w;
            auto w_out = this->ceil_log2(max_val + 1);
            auto &bit_value = abs_coeff_bits[w];
            // temp_or = carry OR bit_value
            auto temp_or = ++this->variable_counter;
            this->create_new_variable(this->variable_counter);
            this->create_2x1_or(carry_bit, bit_value, temp_or);
            // temp_and[x] = carry AND val[x]
            std::vector<int> temp_and(w_in);
            for (int x = 0; x < w_in; x++) {
                temp_and[x] = ++this->variable_counter;
                this->create_new_variable(this->variable_counter);
                this->create_2x1_and(carry_bit, val[x], temp_and[x]);
            }
            // val_new = temp_or + temp_and
            auto add_carry = temp_or;
            std::vector<int> val_new(w_out);
            for (int x = 0; x < w_in; x++) {
                this->full_adder_coeff_word_size_internal_variables[{idx, w,
                                                                     x}] = val_new[x] = ++this->variable_counter;
                this->create_new_variable(this->variable_counter);
                if (x == w_in - 1) {
                    if (w_in == w_out) {
                        // no carry output required
                        //this->create_half_adder(add_carry, val[x], val_new[x]);
                        this->create_half_adder({add_carry, false}, {val[x], false}, {val_new[x], false});
                    } else {
                        // use carry output as sum output for result MSB
                        int carry_out = ++this->variable_counter;
                        this->create_new_variable(this->variable_counter);
                        //this->create_half_adder(add_carry, val[x], val_new[x], carry_out);
                        this->create_half_adder({add_carry, false}, {val[x], false}, {val_new[x], false},
                                                {carry_out, false});
                        this->full_adder_coeff_word_size_internal_variables[{idx, w, x + 1}] = val_new[x +
                                                                                                       1] = carry_out;
                    }
                } else {
                    // normal half adder
                    int carry_out = ++this->variable_counter;
                    this->create_new_variable(this->variable_counter);
                    //this->create_half_adder(add_carry, val[x], val_new[x], carry_out);
                    this->create_half_adder({add_carry, false}, {val[x], false}, {val_new[x], false},
                                            {carry_out, false});
                    add_carry = carry_out;
                }
            }
            val = val_new;
            carry_bit = temp_or;
        }
        for (int w = 0; w < val.size(); w++) {
            this->full_adder_coeff_word_size_variables[{idx, w}] = val[w];
        }
    }
}

void mcm::create_full_adder_msb_constraints(int idx, formulation_mode mode) {
    if (mode == formulation_mode::only_FA_limit) return;
    for (int v = 0; v <c_row_size(); v++) {
        if (!this->calc_twos_complement) {
            // can always cut MSB because all coefficients are positive
            // -> just set the m to 1 and count on unit propagation within the solver :)
            this->full_adder_msb_variables[idx] = this->init_const_one_bit();
            this->force_bit(this->const_one_bit, 1);
            return;
        }
        auto m = this->full_adder_msb_variables[idx] = ++this->variable_counter;
        this->create_new_variable(this->variable_counter);
        int s_c = this->adder_output_value_variables.at({idx, this->word_size - 1, v});
        if (idx <= (1 + idx_input_buffer()) && idx >0) {
            // m = not s_c
            this->create_1x1_negated_implication(s_c, m);
            this->create_1x1_reversed_negated_implication(s_c, m);
            return;
        }
        int s_x = this->input_select_mux_output_variables.at({idx, mcm::left, this->word_size - 1, v});
        int s_y = this->input_select_mux_output_variables.at({idx, mcm::right, this->word_size - 1, v});
        // create clauses to decide whether the sign m can be copied from one of the inputs
        // 1)
        this->create_arbitrary_clause({
                                              {s_y, false},
                                              {s_c, false},
                                              {m,   false},
                                      });
        // 2)
        this->create_arbitrary_clause({
                                              {s_x, false},
                                              {s_y, true},
                                              {m,   false},
                                      });
        // 3)
        this->create_arbitrary_clause({
                                              {s_x, true},
                                              {s_c, true},
                                              {m,   false},
                                      });
        // 4) -> redundant
        // 5) -> redundant
        // 6) -> redundant
        // 7)
        this->create_arbitrary_clause({
                                              {s_x, true},
                                              {s_y, true},
                                              {s_c, false},
                                              {m,   true},
                                      });
        // 8)
        this->create_arbitrary_clause({
                                              {s_x, false},
                                              {s_y, false},
                                              {s_c, true},
                                              {m,   true},
                                      });
    }
}

void mcm::create_full_adder_coeff_word_size_sum_constraints(int idx, formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	auto num_bits_word_size = this->ceil_log2(this->word_size+1);
	if (idx <= (1 + idx_input_buffer()) && idx > 0) {
		// no addition necessary -> just set container with variables
		for (int w=0; w<num_bits_word_size; w++) {
			this->full_adder_word_size_sum_variables[{idx, w}] = this->full_adder_coeff_word_size_variables.at({idx, w});
		}
		return;
	}
	// result[idx] = num_bits[idx] + result[idx-1]
	auto input_word_size = this->ceil_log2(((idx-1)*this->word_size)+1);
	auto output_word_size = this->ceil_log2((idx*this->word_size)+1);
	std::vector<std::pair<std::vector<int>, bool>> x(2);
	x[0].second = x[1].second = false; // add both bit vectors
	// first input: num_bits[idx]
	x[0].first.resize(num_bits_word_size);
	for (int w=0; w<num_bits_word_size; w++) {
		x[0].first[w] = this->full_adder_coeff_word_size_variables.at({idx, w});
	}
	// second input: result[idx-1]
	x[1].first.resize(input_word_size);
	for (int w=0; w<input_word_size; w++) {
		x[1].first[w] = this->full_adder_word_size_sum_variables.at({idx-1, w});
	}
	// output: result[idx]
	auto output_bits = this->create_bitheap(x);
	for (int w=0; w<output_word_size; w++) {
		this->full_adder_word_size_sum_variables[{idx, w}] = output_bits.at(w);
	}
}

void mcm::create_full_adder_shift_gain_constraints(int idx, formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	for (int w = 0; w < this->shift_word_size; w++) {
		auto var = this->full_adder_shift_gain_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
		/*
		this->create_arbitrary_clause({
																		{this->shift_output_variables.at({idx, w}), false},
																		{var, true},
																	});
		this->create_arbitrary_clause({
																		{this->input_negate_value_variables.at(idx), true},
																		{this->input_negate_select_variables.at(idx), true},
																		{var, true},
																	});
		this->create_arbitrary_clause({
																		{this->shift_output_variables.at({idx, w}), true},
																		{this->input_negate_value_variables.at(idx), false},
																		{var, false},
																	});
		this->create_arbitrary_clause({
																		{this->shift_output_variables.at({idx, w}), true},
																		{this->input_negate_select_variables.at(idx), false},
																		{var, false},
																	});
																	*/
		this->create_arbitrary_clause({
																		{this->input_shift_value_variables.at({idx, w}), false},
																		{var, true},
																	});
		this->create_arbitrary_clause({
																		{this->input_negate_value_variables.at(idx), true},
																		{this->input_negate_select_variables.at(idx), true},
																		{var, true},
																	});
		this->create_arbitrary_clause({
																		{this->input_shift_value_variables.at({idx, w}), true},
																		{this->input_negate_value_variables.at(idx), false},
																		{var, false},
																	});
		this->create_arbitrary_clause({
																		{this->input_shift_value_variables.at({idx, w}), true},
																		{this->input_negate_select_variables.at(idx), false},
																		{var, false},
																	});
	}
}

void mcm::create_full_adder_shift_sum_constraints(int idx, formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	auto num_bits_word_size = this->ceil_log2(this->max_shift+1);
	if (idx <= (1 + idx_input_buffer()) && idx > 0) {
		// no addition necessary -> just set container with variables
		for (int w=0; w<num_bits_word_size; w++) {
			this->full_adder_shift_sum_variables[{idx, w}] = this->full_adder_shift_gain_variables.at({idx, w});
		}
	}
	// result[idx] = shift[idx] + result[idx-1]
	auto input_word_size = this->ceil_log2(((idx-1)*this->max_shift)+1);
	auto output_word_size = this->ceil_log2((idx*this->max_shift)+1);
	std::vector<std::pair<std::vector<int>, bool>> x(2);
	x[0].second = x[1].second = false; // add both bit vectors
	// first input: shift[idx]
	x[0].first.resize(num_bits_word_size);
	for (int w=0; w<num_bits_word_size; w++) {
		x[0].first[w] = this->full_adder_shift_gain_variables.at({idx, w});
	}
	// second input: result[idx-1]
	x[1].first.resize(input_word_size);
	for (int w=0; w<input_word_size; w++) {
		x[1].first[w] = this->full_adder_shift_sum_variables.at({idx-1, w});
	}
	// output: result[idx]
	auto output_bits = this->create_bitheap(x);
	for (int w=0; w<output_word_size; w++) {
		this->full_adder_shift_sum_variables[{idx, w}] = output_bits.at(w);
	}
}

void mcm::create_full_adder_msb_sum_constraints(formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	if (this->num_adders == 1) {
		// no addition necessary -> just set container with variable
		this->full_adder_msb_sum_variables[0] = this->full_adder_msb_variables.at(1);
		return;
	}
	// sum all bits up
	std::vector<std::pair<std::vector<int>, bool>> x(this->num_adders, std::pair<std::vector<int>, bool>(std::vector<int>(1), false));
	for (int idx=1; idx<= (this->num_adders + idx_input_buffer()); idx++) {
		x[idx-1].first[0] = this->full_adder_msb_variables.at(idx);
	}
	auto output_bits = this->create_bitheap(x);
	auto output_word_size = this->ceil_log2(this->num_adders+1);
	for (int w=0; w<output_word_size; w++) {
		this->full_adder_msb_sum_variables[w] = output_bits.at(w);
	}
}

void mcm::create_full_adder_add_subtract_inputs_constraints(formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	// result = shift_result[last_stage] + msb_sum
	auto word_size_left_input = this->ceil_log2(this->num_adders * this->max_shift + 1);
	auto word_size_right_input = this->ceil_log2(this->num_adders + 1);
	auto output_word_size = this->ceil_log2((this->num_adders + 1) * this->max_shift + 1);
	std::vector<std::pair<std::vector<int>, bool>> x(2);
	x[0].second = x[1].second = false;
	x[0].first.resize(word_size_left_input);
	for (int w=0; w<word_size_left_input; w++) {
		x[0].first[w] = this->full_adder_shift_sum_variables.at({this->num_adders, w});
	}
	x[1].first.resize(word_size_right_input);
	for (int w=0; w<word_size_right_input; w++) {
		x[1].first[w] = this->full_adder_msb_sum_variables.at(w);
	}
	auto output_bits = this->create_bitheap(x);
	for (int w=0; w<output_word_size; w++) {
		this->full_adder_add_subtract_inputs_variables[w] = output_bits.at(w);
	}
}

void mcm::create_full_adder_cpa_constraints(formulation_mode mode) {
	if (mode == formulation_mode::only_FA_limit) return;
	// result = num_bits[last_stage] - (shift_result[last_stage] + msb_sum)
	auto input_word_size_add = this->ceil_log2(this->word_size * this->num_adders + 1);
	auto input_word_size_sub = this->ceil_log2((this->num_adders + 1) * this->max_shift + 1);
	auto output_word_size = std::max(input_word_size_add, input_word_size_sub)+1;
	std::vector<std::pair<std::vector<int>, bool>> x(2);
	// decide add/sub for the two inputs
	x[0].second = false;
	x[1].second = true;
	// prepare add input bits
	x[0].first.resize(output_word_size);
	for (int w=0; w<input_word_size_add; w++) {
		x[0].first[w] = this->full_adder_word_size_sum_variables.at({this->num_adders, w});
	}
	// sign extend add input with zeros
	for (int w=input_word_size_add; w<output_word_size; w++) {
		x[0].first[w] = this->init_const_zero_bit();
	}
	// prepare sub input bits
	x[1].first.resize(output_word_size);
	for (int w=0; w<input_word_size_sub; w++) {
		x[1].first[w] = this->full_adder_add_subtract_inputs_variables.at(w);
	}
	// sign extend sub input with zeros
	for (int w=input_word_size_sub; w<output_word_size; w++) {
		x[1].first[w] = this->init_const_zero_bit();
	}
	// compute output
	auto output_bits = this->create_bitheap(x);
	for (int w=0; w<output_word_size; w++) {
		this->full_adder_result_variables[w] = output_bits.at(w);
	}
}

void mcm::create_full_adder_result_constraints() {
	// force num_full_adders <= max_full_adders
	auto input_word_size_add = this->ceil_log2(this->word_size * this->num_adders + 1);
	auto input_word_size_sub = this->ceil_log2((this->num_adders + 1) * this->max_shift + 1);
	auto output_word_size = std::max(input_word_size_add, input_word_size_sub)+1;
	int ok_last = -1;
	int carry_last = -1;
	for (int w=output_word_size-1; w>=0; w--) {
		auto c = (int)((this->max_full_adders >> w) & 1);
		int x = this->full_adder_result_variables.at(w);
		int ok_new = -1;
		int carry_new = -1;
		if (w == output_word_size-1) {
			// sign bit
			if (c) {
				ok_new = this->full_adder_result_variables.at(w);
				if (w != 0) {
					carry_new = this->full_adder_result_variables.at(w);
				}
			}
			else {
				ok_new = this->init_const_one_bit();
				if (w != 0) {
					carry_new = ++this->variable_counter;
					this->create_new_variable(this->variable_counter);
					// carry_new = not x (via negated implications)
					this->create_1x1_negated_implication(x, carry_new);
					this->create_1x1_reversed_negated_implication(x, carry_new);
				}
			}
		}
		else {
			// regular bit
			ok_new = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
			if (w != 0) {
				carry_new = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
			}
			if (c) {
				this->create_2x1_or(ok_last, carry_last, ok_new);
				if (w != 0) {
					this->create_2x1_and(x, carry_last, carry_new);
				}
			}
			else {
				auto not_x = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
				this->create_1x1_negated_implication(x, not_x);
				this->create_1x1_reversed_negated_implication(x, not_x);
				this->create_2x1_mux(ok_last, not_x, carry_last, ok_new);
				if (w != 0) {
					this->create_2x1_and(carry_last, not_x, carry_new);
				}
			}
		}
		// force ok bit of this stage to 1
		this->force_bit(ok_new, 1);
		// pass ok and carry bits to next stage
		this->full_adder_comparator_ok_variables[w] = ok_last = ok_new;
		this->full_adder_comparator_carry_variables[w] = carry_last = carry_new;
	}
}

void mcm::create_full_adder(std::pair<int, bool> a, std::pair<int, bool> b, std::pair<int, bool> c_i, std::pair<int, bool> sum, std::pair<int, bool> c_o) {
	//this->create_add_sum(a, b, c_i, sum);
	// 1)
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, not b.second}, {c_i.first, c_i.second}, {sum.first, sum.second}});
	// 2)
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, b.second}, {c_i.first, c_i.second}, {sum.first, sum.second}});
	// 3)
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, b.second}, {c_i.first, c_i.second}, {sum.first, not sum.second}});
	// 4)
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, not b.second}, {c_i.first, c_i.second}, {sum.first, not sum.second}});
	// 5)
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, not b.second}, {c_i.first, not c_i.second}, {sum.first, not sum.second}});
	// 6)
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, b.second}, {c_i.first, not c_i.second}, {sum.first, not sum.second}});
	// 7)
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, b.second}, {c_i.first, not c_i.second}, {sum.first, sum.second}});
	// 8)
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, not b.second}, {c_i.first, not c_i.second}, {sum.first, sum.second}});

	if (c_o.first <= 0) return;
	//this->create_add_carry(a, b, c_i, c_o);
	// 1)
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, not b.second}, {c_o.first, c_o.second}});
	// 2)
	this->create_arbitrary_clause({{a.first, a.second}, {c_i.first, c_i.second}, {c_o.first, not c_o.second}});
	// 3)
	this->create_arbitrary_clause({{b.first, b.second}, {c_i.first, c_i.second}, {c_o.first, not c_o.second}});
	// 4)
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, b.second}, {c_o.first, not c_o.second}});
	// 5)
	this->create_arbitrary_clause({{b.first, not b.second}, {c_i.first, not c_i.second}, {c_o.first, c_o.second}});
	// 6)
	this->create_arbitrary_clause({{a.first, not a.second}, {c_i.first, not c_i.second}, {c_o.first, c_o.second}});
}

void mcm::create_half_adder(std::pair<int, bool> a, std::pair<int, bool> b, std::pair<int, bool> sum, std::pair<int, bool> c_o) {
	//this->create_2x1_xor(a, b, sum);
	// 1) a b -sum
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, b.second}, {sum.first, not sum.second}});
	// 2) a -b sum
	this->create_arbitrary_clause({{a.first, a.second}, {b.first, not b.second}, {sum.first, sum.second}});
	// 3) -a b sum
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, b.second}, {sum.first, sum.second}});
	// 4) -a -b -sum
	this->create_arbitrary_clause({{a.first, not a.second}, {b.first, not b.second}, {sum.first, not sum.second}});

	if (c_o.first <= 0) return;
	//this->create_2x1_and(a, b, c_o);
	// 1) -a -b c_o
	this->create_arbitrary_clause({{a.first, not a.second},{b.first, not b.second},{c_o.first, c_o.second}});
	// 2) a -c_o
	this->create_arbitrary_clause({{a.first, a.second},{c_o.first, not c_o.second}});
	// 3) b -c_o
	this->create_arbitrary_clause({{b.first, b.second},{c_o.first, not c_o.second}});
}

int mcm::init_const_one_bit() {
	if (this->const_one_bit < 1) {
		this->const_one_bit = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
		this->force_bit(this->const_one_bit, 1);
	}
	return this->const_one_bit;
}

int mcm::init_const_zero_bit() {
	if (this->const_zero_bit < 1) {
		this->const_zero_bit = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
		this->force_bit(this->const_zero_bit, 0);
	}
	return this->const_zero_bit;
}

std::vector<int> mcm::create_bitheap(const std::vector<std::pair<std::vector<int>, bool>> &x) {
	std::vector<int> result_variables;
	std::map<int, std::vector<std::pair<int, bool>>> y;
	int num_bits = 0;
	for (auto &it : x) {
		auto bits = it.first;
		auto sub = it.second;
		if (bits.size() > num_bits) num_bits = bits.size();
		if (sub) {
			// add 1 for 2k inversion
			y[0].emplace_back(this->init_const_one_bit(), false);
			// add inverted bits
			for (int bit_pos=0; bit_pos<bits.size(); bit_pos++) {
				y[bit_pos].emplace_back(bits[bit_pos], true);
			}
		}
		else {
			// add bits
			for (int bit_pos=0; bit_pos<bits.size(); bit_pos++) {
				y[bit_pos].emplace_back(bits[bit_pos], false);
			}
		}
	}
	int i = 0;
	while (i < num_bits) {
		while (y[i].size() > 1) {
			if (y[i].size() == 2) {
				// half adder
				// create new literals for sum and carry
				auto sum = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
				auto carry = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
				// get bits to add from container
				auto a = y[i].back();
				y[i].pop_back();
				auto b = y[i].back();
				y[i].pop_back();
				// create clauses
				this->create_half_adder(a, b, {sum, false}, {carry, false});
				// add new bits to bitheap
				y[i].emplace_back(sum, false);
				y[i+1].emplace_back(carry, false);
				if (i+2 > num_bits) num_bits = i+2;
			}
			else {
				// full adder
				// create new literals for sum and carry
				auto sum = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
				auto carry = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
				// get bits to add from container
				auto a = y[i].back();
				y[i].pop_back();
				auto b = y[i].back();
				y[i].pop_back();
				auto c = y[i].back();
				y[i].pop_back();
				// create clauses
				this->create_full_adder(a, b, c, {sum, false}, {carry, false});
				y[i].emplace_back(sum, false);
				y[i+1].emplace_back(carry, false);
				if (i+2 > num_bits) num_bits = i+2;
			}
		}
		if (y[i].size() != 1) {
			std::cerr << "Failed compressing bits at position " << i << " -> " << y[i].size() << " bits are left instead of 1" << std::endl;
			throw std::runtime_error("error during bitheap clause generation");
		}
		if (y[i][0].second) {
			std::cerr << "Failed compressing bits at position " << i << " -> the output bit is inverted..." << std::endl;
			throw std::runtime_error("error during bitheap clause generation");
		}
		result_variables.emplace_back(y[i][0].first);
		// advance to next bit position
		i++;
	}
	return result_variables;
}

void mcm::prohibit_current_solution() {
	std::vector<std::pair<int, bool>> clause;
	int v; // variable buffer
	for (int i=1; i<=this->num_adders; i++) {
		auto mux_word_size = this->ceil_log2(i);
		// left input mux
		for (int s=0; s<mux_word_size; s++) {
			v = this->input_select_selection_variables.at({i, mcm::left, s});
			clause.emplace_back(v, this->get_result_value(v));
		}
		// right input mux
		for (int s=0; s<mux_word_size; s++) {
			v = this->input_select_selection_variables.at({i, mcm::right, s});
			clause.emplace_back(v, this->get_result_value(v));
		}
		// pre add shift
		for (int s=0; s<this->shift_word_size; s++) {
			v = this->input_shift_value_variables.at({i, s});
			clause.emplace_back(v, this->get_result_value(v));
		}
		// negate select
		v = this->input_negate_select_variables.at(i);
		clause.emplace_back(v, this->get_result_value(v));
		// negate value
		v = this->input_negate_value_variables.at(i);
		clause.emplace_back(v, this->get_result_value(v));
		// post add shift
		if (this->enable_node_output_shift) {
			for (int s=0; s<this->shift_word_size; s++) {
				v = this->input_post_adder_shift_value_variables.at({i, s});
				clause.emplace_back(v, this->get_result_value(v));
			}
		}
	}
	this->create_arbitrary_clause(clause);
}

void mcm::set_enumerate_all(bool new_enumerate_all) {
	this->enumerate_all = new_enumerate_all;
}

void mcm::solve_enumeration() {
	this->num_FA_opt = true;
	this->num_add_opt = true;
	if (this->verbosity == verbosity_mode::debug_mode) {
		if (this->c_row_size() == 1 and this->c_column_size() == 1) {
			std::cout << "Trying to solve SCM problem for constant: " << this->C[0].front() << std::endl;
		}
		else if (this->c_row_size() == 1){
			std::cout << "Trying to solve MCM problem for following constants: ";
			for (auto &v : this->C) {
			    for (auto &c : v)
				std::cout << "  " << c;
			}
			std::cout << std::endl;
		}
		else{
            std::cout << "Trying to solve CMM problem for following vectors: ";
            for (auto &v : this->C) {
                std::cout << "<";
                for (auto &c : v) {
                    std::cout << "  " << c;
                }
                std::cout << " >";
            }
            std::cout << std::endl;
        }
		std::cout << "with word size " << this->word_size << " and max shift " << this->max_shift << std::endl;
	}
	bool trivial = true;
    //unit vectors are trivial
    for (auto &v : this->C) {
        if (std::accumulate(v.begin(), v.end(), 0) != 1) {
            trivial = false;
            break;
        }
    }
	if (trivial) {
		this->found_solution = true;
		this->ran_into_timeout = false;
		this->output_values[{0,0}] = 1;
		return;
	}
	formulation_mode mode = formulation_mode::reset_all;
	while (!this->found_solution) {
		this->fa_minimization_timeout = this->timeout;
		++this->num_adders;
		this->optimization_loop(mode);
		if (this->ran_into_timeout) {
			// timeout => can't say anything about optimality
			this->num_add_opt = false;
		}
	}
	// now enumerate all solutions for minimum adder count
	mode = formulation_mode::all_FA_clauses;
	while (this->found_solution) {
		this->timeout = this->fa_minimization_timeout;
		// count current # of full adders
		// except for the last node because its output always has a constant number of full adders
		int current_full_adders = 0;
		int MSBs_cut = 0;
		int MSBs_not_cut = 0;
        for(int v=0; v<c_row_size(); v++) {
            for (int idx = 1; idx <= (this->num_adders + idx_input_buffer()); idx++) {
                if (this->output_values.at({idx,v}) == 0) {
                    // more adders allocated than necessary
                    continue;
                }
                int FAs_for_this_node = (int) std::ceil(std::log2(std::abs(this->add_result_values.at({idx,v}))));
                int shifter_input_non_zero_LSBs = 0;
                int shifter_input = 1;
                if (idx > (1 + idx_input_buffer())) {
                    shifter_input = this->input_select_mux_output.at({idx, mcm::left, v});
                }
                while ((shifter_input & 1) == 0) {
                    shifter_input = shifter_input >> 1;
                    shifter_input_non_zero_LSBs++;
                }
                if (this->subtract.at(idx) == 0 or this->negate_select.at(idx) == 0) {
                    // we do not need to use a full adder for the shifted LSBs
                    //   for a + b
                    //   and a - (b << s)
                    FAs_for_this_node -= (this->shift_value.at(idx) + shifter_input_non_zero_LSBs);
                }
                auto can_cut_MSB =
                        (this->output_values.at({idx,v}) >= 0 and this->input_select_mux_output[{idx, mcm::left, v}] >= 0) or
                        (this->output_values.at({idx,v}) >= 0 and this->input_select_mux_output[{idx, mcm::right, v}] >= 0) or
                        (this->output_values.at({idx,v}) < 0 and this->input_select_mux_output[{idx, mcm::left, v}] < 0) or
                        (this->output_values.at({idx,v}) < 0 and this->input_select_mux_output[{idx, mcm::right, v}] < 0);
                if (can_cut_MSB) {
                    MSBs_cut++;
                } else {
                    MSBs_not_cut++;
                }
                if (this->verbosity != verbosity_mode::quiet_mode) {
                    std::cout << "Additional FAs for node " << idx << " = "
                              << (can_cut_MSB ? FAs_for_this_node - 1 : FAs_for_this_node) << std::endl;
                }
                current_full_adders += (FAs_for_this_node - ((int) can_cut_MSB));
            }
        }
		if (current_full_adders > this->max_full_adders and this->max_full_adders != FULL_ADDERS_UNLIMITED) {
			if (this->verbosity != verbosity_mode::quiet_mode) {
				this->print_solution();
			}
			throw std::runtime_error("SAT solver exceeded full adder limit! Limit was "+std::to_string(this->max_full_adders)+" but solver returned solution with "+std::to_string(current_full_adders)+" FAs!");
		}
		else {
			if (this->verbosity != verbosity_mode::quiet_mode) {
				std::cout << "Current solution needs " << current_full_adders << " additional full adders" << std::endl;
				this->print_solution();
			}
		}
		// no full adder limitation since we are interested in ALL solutions
		this->max_full_adders = FULL_ADDERS_UNLIMITED;
		this->prohibit_current_solution(); // ENUMERATE ALL POSSIBLE SOLUTIONS but prohibit the last found one
		this->optimization_loop(mode);
		mode = formulation_mode::only_FA_limit;
		if (this->ran_into_timeout) {
			// timeout => can't say anything about optimality
			this->num_FA_opt = false;
		}
	}
	this->found_solution = true;
}

void mcm::solve_standard() {
	this->num_FA_opt = true;
	this->num_add_opt = true;
    if (this->verbosity == verbosity_mode::debug_mode) {
        if (this->c_row_size() == 1 and this->c_column_size() == 1) {
            std::cout << "Trying to solve SCM problem for constant: " << this->C[0].front() << std::endl;
        }
        else if (this->c_row_size() == 1){
            std::cout << "Trying to solve MCM problem for following constants: ";
            for (auto &v : this->C) {
                for (auto &c : v)
                    std::cout << "  " << c;
            }
            std::cout << std::endl;
        }
        else{
            std::cout << "Trying to solve CMM problem for following vectors: ";
            for (auto &v : this->C) {
                std::cout << "<";
                for (auto &c : v) {
                    std::cout << "  " << c;
                }
                std::cout << " >";
            }
            std::cout << std::endl;
        }
        std::cout << "with word size " << this->word_size << " and max shift " << this->max_shift << std::endl;
    }
    bool trivial = true;
    //unit vectors are trivial
    for (auto &v : this->C) {
        if (std::accumulate(v.begin(), v.end(), 0) != 1) {
            trivial = false;
            break;
        }
    }
	if (trivial) {
		this->found_solution = true;
		this->ran_into_timeout = false;
		this->output_values[{0,0}] = 1;
		return;
	}
	formulation_mode mode = formulation_mode::reset_all;
	while (!this->found_solution) {
		this->fa_minimization_timeout = this->timeout;
		++this->num_adders;
		this->optimization_loop(mode);
		if (this->ran_into_timeout) {
			// timeout => can't say anything about optimality
			this->num_add_opt = false;
		}
	}
	// check if we should even optimize the number of full adders and return if not
	if (!this->minimize_full_adders) {
		this->num_FA_opt = false; // don't know if solution is optimal w.r.t. full adders
		return;
	}

	mode = formulation_mode::all_FA_clauses;
	while (this->found_solution) {
		this->timeout = this->fa_minimization_timeout;
		// count current # of full adders
		// except for the last node because its output always has a constant number of full adders
		int current_full_adders = 0;
		int MSBs_cut = 0;
		int MSBs_not_cut = 0;
        for(int v=0; v<c_row_size(); v++) {
            for (int idx = 1; idx <= (this->num_adders + idx_input_buffer()); idx++) {
                if (this->output_values.at({idx, v}) == 0) {
                    // more adders allocated than necessary
                    continue;
                }
                int FAs_for_this_node = (int) std::ceil(std::log2(std::abs(this->add_result_values.at({idx,v}))));
                int shifter_input_non_zero_LSBs = 0;
                int shifter_input = 1;
                if (idx > (1 + idx_input_buffer())) {
                    shifter_input = this->input_select_mux_output.at({idx, mcm::left, v});
                }
                while ((shifter_input & 1) == 0) {
                    shifter_input = shifter_input >> 1;
                    shifter_input_non_zero_LSBs++;
                }
                if (this->subtract.at(idx) == 0 or this->negate_select.at(idx) == 0) {
                    // we do not need to use a full adder for the shifted LSBs
                    //   for a + b
                    //   and a - (b << s)
                    FAs_for_this_node -= (this->shift_value.at(idx) + shifter_input_non_zero_LSBs);
                }

                auto can_cut_MSB = (this->output_values.at({idx, v}) >= 0 and
                                    this->input_select_mux_output[{idx, mcm::left, v}] >= 0) or
                                   (this->output_values.at({idx, v}) >= 0 and
                                    this->input_select_mux_output[{idx, mcm::right, v}] >= 0) or
                                   (this->output_values.at({idx, v}) < 0 and
                                    this->input_select_mux_output[{idx, mcm::left, v}] < 0) or
                                   (this->output_values.at({idx, v}) < 0 and
                                    this->input_select_mux_output[{idx, mcm::right, v}] < 0);
                if (can_cut_MSB) {
                    MSBs_cut++;
                } else {
                    MSBs_not_cut++;
                }
                if (this->verbosity != verbosity_mode::quiet_mode)
                    std::cout << "Additional FAs for node " << idx << " = " << (can_cut_MSB ? FAs_for_this_node - 1
                                                                                            : FAs_for_this_node)
                              << std::endl;
                current_full_adders += (FAs_for_this_node - ((int) can_cut_MSB));
            }
        }
		if (this->max_full_adders == FULL_ADDERS_UNLIMITED) {
			if (this->verbosity != verbosity_mode::quiet_mode) {
				std::cout << "Initial solution needs " << current_full_adders << " additional full adders" << std::endl;
				this->print_solution();
			}
		}
		else if (current_full_adders > this->max_full_adders) {
			if (this->verbosity != verbosity_mode::quiet_mode) {
				this->print_solution();
			}
			throw std::runtime_error("SAT solver exceeded full adder limit! Limit was "+std::to_string(this->max_full_adders)+" but solver returned solution with "+std::to_string(current_full_adders)+" FAs!");
		}
		else {
			if (this->verbosity != verbosity_mode::quiet_mode) {
				std::cout << "Current solution needs " << current_full_adders << " additional full adders" << std::endl;
				this->print_solution();
			}
		}
		// must add the number of MSBs that could not be cut because the SAT solver allocs an extra LUT for each of them
		this->max_full_adders = current_full_adders - 1;
		if (this->max_full_adders < -(this->num_adders * (this->max_shift+1))) {
			// trivial minimum value reached
			return;
		}
		this->optimization_loop(mode);
		mode = formulation_mode::only_FA_limit;
		if (this->ran_into_timeout) {
			// timeout => can't say anything about optimality
			this->num_FA_opt = false;
		}
	}
	this->found_solution = true;
}

