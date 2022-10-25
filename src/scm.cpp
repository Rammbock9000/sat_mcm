//
// Created by nfiege on 9/26/22.
//

#include "scm.h"
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>

scm::scm(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers, bool write_cnf)
	:	C(C), timeout(timeout), quiet(quiet), threads(threads), write_cnf(write_cnf) {
	// make it even and count shift
	this->calc_twos_complement = allow_negative_numbers;
	for (auto &c : this->C) {
		// ignore 0
		if (c == 0) continue;
		auto original_number = c;
		int shifted_bits = 0;
		// handle negative numbers
		if (c < 0) {
			c = -c;
			this->negative_coeff_requested[c] = true;
		}
		// right shift until odd
		while ((c & 1) == 0) {
			c = c / 2; // do not use shift operation because it is not uniquely defined for negative numbers
			shifted_bits++;
		}
		this->requested_constants[original_number] = {c, shifted_bits};
	}
	// set word sizes & track unique constants
	this->word_size = 1;
	std::set<int> non_one_unique_constants;
	for (auto &c : this->C) {
		if (c != 1 and c != 0) non_one_unique_constants.insert(c);
		auto w = this->ceil_log2(std::abs(c))+1;
		if (w > this->word_size) this->word_size = w;
	}
	this->max_shift = this->word_size-1;
	if (this->calc_twos_complement) {
		// account for sign bit
		this->word_size++;
	}
	this->shift_word_size = this->ceil_log2(this->max_shift+1);
	this->num_adders = (int)non_one_unique_constants.size()-1;
	std::cout << "min num adders = " << this->num_adders+1 << std::endl;
	// set constants vector
	this->C.clear();
	for (auto &c : non_one_unique_constants) {
		this->C.emplace_back(c);
	}
}

void scm::optimization_loop() {
	auto start_time = std::chrono::steady_clock::now();
	if (!this->quiet) std::cout << "  resetting backend now" << std::endl;
	this->reset_backend();
	if (!this->quiet) std::cout << "  constructing problem for " << this->num_adders << " adders" << (this->max_full_adders>=0?" and "+std::to_string(this->max_full_adders)+" full adders":"") << std::endl;
	this->construct_problem();
	if (!this->quiet) std::cout << "  start solving with " << this->variable_counter << " variables and " << this->constraint_counter << " constraints" << std::endl;
	auto [a, b] = this->check();
	auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time).count() / 1000.0;
	this->found_solution = a;
	this->ran_into_timeout = b;
	if (this->found_solution) {
		std::cout << "  found solution for #adders = " << this->num_adders << (this->max_full_adders>=0?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " after " << elapsed_time << " seconds 8-)" << std::endl;
		this->get_solution_from_backend();
		if (this->solution_is_valid()) {
			std::cout << "Solution is verified :-)" << std::endl;
		}
		else {
			throw std::runtime_error("Solution is invalid (found bug) :-(");
		}
	}
	else if (this->ran_into_timeout) {
		std::cout << "  ran into timeout for #adders = " << this->num_adders << (this->max_full_adders>=0?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " after " << elapsed_time << " seconds :-(" << std::endl;
	}
	else {
		std::cout << "  problem for #adders = " << this->num_adders << (this->max_full_adders>=0?" and max. "+std::to_string(this->max_full_adders)+" full adders":"") << " is proven to be infeasible after " << elapsed_time << " seconds... " << (this->max_full_adders>=0?"":"keep trying :-)") << std::endl;
	}
}

void scm::solve() {
	this->num_FA_opt = true;
	this->num_add_opt = true;
	if (!this->quiet) {
		std::cout << "trying to solve SCM problem for following constants: ";
		for (auto &c : this->C) {
			std::cout << "  " << c << std::endl;
		}
		std::cout << "with word size " << this->word_size << " and max shift " << this->max_shift << std::endl;
	}
	bool trivial = true;
	for (auto &c : this->C) {
		if (c != 1) {
			trivial = false;
			break;
		}
	}
	if (trivial) {
		this->found_solution = true;
		this->ran_into_timeout = false;
		this->output_values[0] = 1;
		return;
	}
	while (!this->found_solution) {
		++this->num_adders;
		this->optimization_loop();
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
	while (this->found_solution) {
		// count current # of full adders
		// except for the last node because its output always has a constant number of full adders
		int current_full_adders = 0;
		for (int idx = 1; idx <= this->num_adders; idx++) {
			if (this->output_values.at(idx) == 0) {
				// more adders allocated than necessary
				continue;
			}
			int FAs_for_this_node = (int)std::ceil(std::log2(std::abs(this->add_result_values.at(idx))));
			int shifter_input_non_zero_LSBs = 0;
			int shifter_input = 1;
			if (idx > 1 and this->shift_input_select.at(idx) == 1) {
				// left input
				shifter_input = this->input_select_mux_output.at({idx, scm::left});
			}
			else if (idx > 1 and this->shift_input_select.at(idx) == 0) {
				// right input
				shifter_input = this->input_select_mux_output.at({idx, scm::right});
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
			/*
			else if (this->negate_select.at(idx) == 1) {
				// for subtractions, we do not need to use a full adder
				// for the shifted LSBs if the non-shifted input gets subtracted
				// except for the LSB where we need a HA to account for the carry-in
				// the remaining shifted bits can be handled by the carry chain alone
				FAs_for_this_node -= ((this->shift_value.at(idx) == 0 ? 0 : this->shift_value.at(idx) - 1) + shifter_input_non_zero_LSBs);
			}
			 */
			std::cout << "FAs for node " << idx << " = " << FAs_for_this_node << std::endl;
			current_full_adders += FAs_for_this_node;
		}
		if (this->max_full_adders < 0) {
			std::cout << "Initial solution needs " << current_full_adders << " full adders" << std::endl;
			this->print_solution();
		}
		else if (current_full_adders > this->max_full_adders) {
			this->print_solution();
			throw std::runtime_error("SAT solver exceeded full adder limit! Limit was "+std::to_string(this->max_full_adders)+" but solver returned solution with "+std::to_string(current_full_adders)+" FAs!");
		}
		else {
			std::cout << "Current solution needs " << current_full_adders << " full adders" << std::endl;
			this->print_solution();
		}
		this->max_full_adders = current_full_adders - 1;
		if (this->max_full_adders < 0) {
			return;
		}
		this->optimization_loop();
		if (this->ran_into_timeout) {
			// timeout => can't say anything about optimality
			this->num_FA_opt = false;
		}
	}
	this->found_solution = true;
}

void scm::reset_backend() {
	this->constraint_counter = 0;
	this->variable_counter = 0;
	this->cnf_clauses.str("");
}

void scm::construct_problem() {
	if (!this->quiet) std::cout << "    creating variables now" << std::endl;
	this->create_variables();
	if (!this->quiet) std::cout << "    creating constraints now" << std::endl;
	this->create_constraints();
	if (this->write_cnf) {
		if (!this->quiet) std::cout << "    creating cnf file now" << std::endl;
		this->create_cnf_file();
	}
}

void scm::create_variables() {
	if (!this->quiet) std::cout << "      creating input node variables" << std::endl;
	this->create_input_node_variables();
	for (int i=1; i<=this->num_adders; i++) {
		if (!this->quiet) std::cout << "      creating variables for node " << i << std::endl;
		if (!this->quiet) std::cout << "        create_input_select_mux_variables" << std::endl;
		this->create_input_select_mux_variables(i);
		if (!this->quiet) std::cout << "        create_input_select_selection_variables" << std::endl;
		this->create_input_select_selection_variables(i);
		if (!this->quiet) std::cout << "        create_input_shift_select_variable" << std::endl;
		this->create_input_shift_select_variable(i);
		if (!this->quiet) std::cout << "        create_shift_select_output_variables" << std::endl;
		this->create_shift_select_output_variables(i);
		if (!this->quiet) std::cout << "        create_input_shift_value_variables" << std::endl;
		this->create_input_shift_value_variables(i);
		if (!this->quiet) std::cout << "        create_shift_internal_variables" << std::endl;
		this->create_shift_internal_variables(i);
		if (!this->quiet) std::cout << "        create_input_negate_select_variable" << std::endl;
		this->create_input_negate_select_variable(i);
		if (!this->quiet) std::cout << "        create_negate_select_output_variables" << std::endl;
		this->create_negate_select_output_variables(i);
		if (!this->quiet) std::cout << "        create_input_negate_value_variable" << std::endl;
		this->create_input_negate_value_variable(i);
		if (!this->quiet) std::cout << "        create_xor_output_variables" << std::endl;
		this->create_xor_output_variables(i);
		if (!this->quiet) std::cout << "        create_adder_internal_variables" << std::endl;
		this->create_adder_internal_variables(i);
		if (this->enable_node_output_shift) {
			if (!this->quiet) std::cout << "        create_post_adder_input_shift_value_variables" << std::endl;
			this->create_post_adder_input_shift_value_variables(i);
			if (!this->quiet) std::cout << "        create_post_adder_shift_variables" << std::endl;
			this->create_post_adder_shift_variables(i);
		}
		if (!this->quiet) std::cout << "        create_output_value_variables" << std::endl;
		this->create_output_value_variables(i);
		if (this->C.size() != 1 or (this->calc_twos_complement and this->sign_inversion_allowed[this->C[0]])) {
			if (!this->quiet) std::cout << "        create_mcm_output_variables" << std::endl;
			this->create_mcm_output_variables(i);
		}
		if (this->max_full_adders > 0 and i <= this->num_adders) {
			if (!this->quiet) std::cout << "        create_full_adder_alloc_variables" << std::endl;
			this->create_full_adder_alloc_variables(i);
		}
	}
}

void scm::create_constraints() {
	if (!this->quiet) std::cout << "      create_input_output_constraints" << std::endl;
	this->create_input_output_constraints();
	for (int i=1; i<=this->num_adders; i++) {
		if (!this->quiet) std::cout << "      creating constraints for node " << i << std::endl;
		if (!this->quiet) std::cout << "        create_input_select_constraints" << std::endl;
		this->create_input_select_constraints(i);
		if (!this->quiet) std::cout << "        create_input_select_limitation_constraints" << std::endl;
		this->create_input_select_limitation_constraints(i);
		if (!this->quiet) std::cout << "        create_shift_limitation_constraints" << std::endl;
		this->create_shift_limitation_constraints(i);
		if (!this->quiet) std::cout << "        create_shift_select_constraints" << std::endl;
		this->create_shift_select_constraints(i);
		if (!this->quiet) std::cout << "        create_shift_constraints" << std::endl;
		this->create_shift_constraints(i);
		if (!this->quiet) std::cout << "        create_negate_select_constraints" << std::endl;
		this->create_negate_select_constraints(i);
		if (!this->quiet) std::cout << "        create_xor_constraints" << std::endl;
		this->create_xor_constraints(i);
		if (!this->quiet) std::cout << "        create_adder_constraints" << std::endl;
		this->create_adder_constraints(i);
		if (this->enable_node_output_shift) {
			if (!this->quiet) std::cout << "        create_post_adder_shift_limitation_constraints" << std::endl;
			this->create_post_adder_shift_limitation_constraints(i);
			if (!this->quiet) std::cout << "        create_post_adder_shift_constraints" << std::endl;
			this->create_post_adder_shift_constraints(i);
		}
		if (this->max_full_adders >= 0 and i <= this->num_adders) {
			if (!this->quiet) std::cout << "        create_full_adder_allocation_constraints" << std::endl;
			this->create_full_adder_allocation_constraints(i);
			if (!this->quiet) std::cout << "        create_full_adder_overlap_constraints" << std::endl;
			this->create_full_adder_overlap_constraints(i);
		}
	}
}

void scm::create_input_node_variables() {
	for (int i=0; i<this->word_size; i++) {
		this->output_value_variables[{0, i}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void scm::create_input_select_mux_variables(int idx) {
	if (idx == 1) return;
	auto select_word_size = this->ceil_log2(idx);
	auto num_muxs = (1 << select_word_size) - 1;
	for (auto &dir : input_directions) {
		for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
			for (int w = 0; w < this->word_size; w++) {
				this->input_select_mux_variables[{idx, dir, mux_idx, w}] = ++this->variable_counter;
				this->create_new_variable(this->variable_counter);
			}
		}
	}
}

void scm::create_input_select_selection_variables(int idx) {
	if (idx == 1) return;
	auto select_word_size = this->ceil_log2(idx);
	for (auto &dir : input_directions) {
		for (int w = 0; w < select_word_size; w++) {
			this->input_select_selection_variables[{idx, dir, w}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_input_shift_select_variable(int idx) {
	if (idx == 1) return;
	this->input_shift_select_variables[idx] = ++this->variable_counter;
	this->create_new_variable(this->variable_counter);
}

void scm::create_shift_select_output_variables(int idx) {
	if (idx == 1) return;
	for (auto &dir : input_directions) {
		for (int w = 0; w < this->word_size; w++) {
			this->shift_select_output_variables[{idx, dir, w}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_input_shift_value_variables(int idx) {
	for (int w = 0; w < this->shift_word_size; w++) {
		this->input_shift_value_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void scm::create_shift_internal_variables(int idx) {
	for (int mux_stage = 0; mux_stage < this->shift_word_size; mux_stage++) {
		for (int w = 0; w < this->word_size; w++) {
			this->shift_internal_mux_output_variables[{idx, mux_stage, w}] = ++this->variable_counter;
			if (mux_stage == this->shift_word_size-1) {
				this->shift_output_variables[{idx, w}] = this->variable_counter;
			}
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_post_adder_input_shift_value_variables(int idx) {
	for (int w = 0; w < this->shift_word_size; w++) {
		this->input_post_adder_shift_value_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void scm::create_post_adder_shift_variables(int idx) {
	for (int mux_stage = 0; mux_stage < this->shift_word_size; mux_stage++) {
		for (int w = 0; w < this->word_size; w++) {
			this->post_adder_shift_internal_mux_output_variables[{idx, mux_stage, w}] = ++this->variable_counter;
			if (mux_stage == this->shift_word_size-1) {
				this->post_adder_shift_output_variables[{idx, w}] = this->variable_counter;
			}
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_input_negate_select_variable(int idx) {
	this->input_negate_select_variables[idx] = ++this->variable_counter;
	this->create_new_variable(this->variable_counter);
}

void scm::create_negate_select_output_variables(int idx) {
	for (auto &dir : input_directions) {
		for (int w = 0; w < this->word_size; w++) {
			this->negate_select_output_variables[{idx, dir, w}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_input_negate_value_variable(int idx) {
	this->input_negate_value_variables[idx] = ++this->variable_counter;
	this->create_new_variable(this->variable_counter);
}

void scm::create_xor_output_variables(int idx) {
	for (int w = 0; w < this->word_size; w++) {
		this->xor_output_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void scm::create_adder_internal_variables(int idx) {
	for (int w = 0; w < this->word_size; w++) {
		this->adder_internal_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
		this->adder_output_value_variables[{idx, w}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
	}
}

void scm::create_output_value_variables(int idx) {
	for (int w = 0; w < this->word_size; w++) {
		if (this->enable_node_output_shift) {
			this->output_value_variables[{idx, w}] = this->post_adder_shift_output_variables.at({idx, w});
		}
		else {
			this->output_value_variables[{idx, w}] = this->adder_output_value_variables.at({idx, w});
		}
	}
}

int scm::ceil_log2(int n) {
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

int scm::floor_log2(int n) {
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

void scm::create_new_variable(int idx) {
	(void) idx; // just do nothing -> should be overloaded by backend if a variable must be explicitly created
}

void scm::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	this->constraint_counter++;
	if (!this->write_cnf) return;
	for (const auto &it : a) {
		this->cnf_clauses << (it.second?-it.first:it.first) << " ";
	}
	this->cnf_clauses << " 0" << std::endl;
}

void scm::create_signed_shift_overflow_protection(int sel, int s_a, int a) {
	// 1)
	this->create_arbitrary_clause({{sel, true}, {s_a, true}, {a, false}});
	// 2)
	this->create_arbitrary_clause({{sel, true}, {s_a, false}, {a, true}});
}

void scm::create_signed_add_overflow_protection(int sub, int s_a, int s_b, int s_y) {
	// 1)
	this->create_arbitrary_clause({{sub, false}, {s_a, false}, {s_b, false}, {s_y, true}});
	// 2)
	this->create_arbitrary_clause({{sub, false}, {s_a, true}, {s_b, true}, {s_y, false}});
	// 3)
	this->create_arbitrary_clause({{sub, true}, {s_a, false}, {s_b, true}, {s_y, true}});
	// 4)
	this->create_arbitrary_clause({{sub, true}, {s_a, true}, {s_b, false}, {s_y, false}});
}

void scm::create_or(std::vector<int> &x) {
	std::vector<std::pair<int, bool>> v(x.size());
	for (auto i=0; i<x.size(); i++) {
		auto &val = x[i];
		if (val > 0) v[i] = {val, false};
		else v[i] = {-val, true};
	}
	this->create_arbitrary_clause(v);
}

void scm::create_1x1_implication(int a, int b) {
	this->create_arbitrary_clause({{a, true}, {b, false}});
}

void scm::create_1x1_negated_implication(int a, int b) {
	this->create_arbitrary_clause({{a, true}, {b, true}});
}

void scm::create_1xN_implication(int a, const std::vector<int> &b) {
	std::vector<std::pair<int, bool>> v(b.size()+1);
	for (auto i=0; i<b.size(); i++) {
		v[i] = {b[i], false};
	}
	v[b.size()] = {a, true};
	this->create_arbitrary_clause(v);
}

void scm::create_MxN_implication(const std::vector<int> &a, const std::vector<int> &b) {
	std::vector<std::pair<int, bool>> v(b.size()+a.size());
	for (auto i=0; i<a.size(); i++) {
		v[i] = {a[i], true};
	}
	for (auto i=0; i<b.size(); i++) {
		v[i+a.size()] = {b[i], false};
	}
	this->create_arbitrary_clause(v);
}

void scm::create_1x1_equivalence(int x, int y) {
	// 1)
	this->create_arbitrary_clause({{x, true}, {y, false}});
	// 2)
	this->create_arbitrary_clause({{x, false}, {y, true}});
}

void scm::create_2x1_mux(int a, int b, int s, int y) {
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

void scm::create_2x1_mux_shift_disallowed(int a, int b, int s, int y) {
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

void scm::create_2x1_mux_zero_const(int a, int s, int y) {
	// 1)
	this->create_arbitrary_clause({{s, true}, {y, true}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {y, true}});
	// 3)
	this->create_arbitrary_clause({{a, true}, {s, false}, {y, false}});
}

void scm::create_2x1_xor(int a, int b, int y) {
	// 1)
	this->create_arbitrary_clause({{a, false}, {b, false}, {y, true}});
	// 2)
	this->create_arbitrary_clause({{a, false}, {b, true}, {y, false}});
	// 3)
	this->create_arbitrary_clause({{a, true}, {b, false}, {y, false}});
	// 4)
	this->create_arbitrary_clause({{a, true}, {b, true}, {y, true}});
}

void scm::create_add_sum(int a, int b, int c_i, int s) {
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

void scm::create_add_carry(int a, int b, int c_i, int c_o) {
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

void scm::force_bit(int x, int val) {
	this->create_arbitrary_clause({{x, val != 1}});
}

void scm::forbid_number(const std::vector<int> &x, int val) {
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

void scm::force_number(const std::vector<int> &x, int val) {
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

std::pair<bool, bool> scm::check() {
	throw std::runtime_error("check is impossible in base class");
}

int scm::get_result_value(int var_idx) {
	throw std::runtime_error("get_result_value is impossible in base class");
}

void scm::create_input_output_constraints() {
	std::vector<int> input_bits(this->word_size);
	std::vector<int> output_bits(this->word_size);
	for (auto w=0; w<this->word_size; w++) {
		input_bits[w] = this->output_value_variables.at({0, w});
		output_bits[w] = this->output_value_variables.at({this->num_adders, w});
	}
	// force input to 1 and output to C
	this->force_number(input_bits, 1);
	if (this->C.size() == 1 and (!this->calc_twos_complement or !this->sign_inversion_allowed[this->C[0]])) {
		// SCM
		this->force_number(output_bits, this->C[0]);
	}
	else {
		// MCM
		this->create_mcm_output_constraints();
	}
}

void scm::create_input_select_constraints(int idx) {
	// stage 1 has no input MUX because it can only be connected to the input node with idx=0
	if (idx == 1) return;
	// create constraints for all muxs
	if (!this->quiet) std::cout << "creating input select constraints for node #" << idx << std::endl;
	auto select_word_size = this->ceil_log2(idx);
	for (auto &dir : this->input_directions) {
		int mux_idx = 0;
		for (int mux_stage = 0; mux_stage < select_word_size; mux_stage++) {
			auto num_muxs_per_stage = (1 << mux_stage);
			auto mux_select_var_idx = this->input_select_selection_variables.at({idx, dir, select_word_size-mux_stage-1}); // mux_stage
			for (int mux_idx_in_stage = 0; mux_idx_in_stage < num_muxs_per_stage; mux_idx_in_stage++) {
				if (mux_stage == select_word_size-1) {
					// connect with another node output
					auto zero_input_node_idx = 2 * mux_idx_in_stage;
					auto one_input_node_idx = zero_input_node_idx + 1;
					if (zero_input_node_idx >= idx) zero_input_node_idx = idx-1;
					if (one_input_node_idx >= idx) one_input_node_idx = idx-1;
					for (int w = 0; w < this->word_size; w++) {
						auto mux_output_var_idx = this->input_select_mux_variables.at({idx, dir, mux_idx, w});
						auto zero_input_var_idx = this->output_value_variables.at({zero_input_node_idx, w});
						auto one_input_var_idx = this->output_value_variables.at({one_input_node_idx, w});
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
						auto mux_output_var_idx = this->input_select_mux_variables.at({idx, dir, mux_idx, w});
						auto zero_input_var_idx = this->input_select_mux_variables.at({idx, dir, zero_input_mux_idx, w});
						auto one_input_var_idx = this->input_select_mux_variables.at({idx, dir, one_input_mux_idx, w});
						this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, mux_select_var_idx, mux_output_var_idx);
					}
				}
				// increment current mux idx
				mux_idx++;
			}
		}
	}
}

void scm::create_shift_constraints(int idx) {
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
			auto mux_output_var_idx = this->shift_internal_mux_output_variables.at({idx, stage, w});
			if (stage == 0) {
				// connect shifter inputs
				if (idx == 1) {
					// shifter input is the output of the input node with idx = 0
					zero_input_var_idx = this->output_value_variables.at({0, w});
					zero_input_sign_bit_idx = this->output_value_variables.at({0, this->word_size-1});
					if (!connect_zero_const) {
						one_input_var_idx = this->output_value_variables.at({0, w_prev});
					}
				}
				else {
					// shifter input is the left shift-select-mux
					zero_input_var_idx = this->shift_select_output_variables.at({idx, scm::left, w});
					zero_input_sign_bit_idx = this->shift_select_output_variables.at({idx, scm::left, this->word_size-1});
					if (!connect_zero_const) {
						one_input_var_idx = this->shift_select_output_variables.at({idx, scm::left, w_prev});
					}
				}
			}
			else {
				// connect output of previous stage
				zero_input_var_idx = this->shift_internal_mux_output_variables.at({idx, stage-1, w});
				zero_input_sign_bit_idx = this->shift_internal_mux_output_variables.at({idx, stage-1, this->word_size-1});
				if (!connect_zero_const) {
					one_input_var_idx = this->shift_internal_mux_output_variables.at({idx, stage-1, w_prev});
				}
			}
			if (w >= first_disallowed_shift_bit) {
				if (this->calc_twos_complement) {
					if (connect_zero_const) {
						this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx, mux_output_var_idx);
					}
					else {
						this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx, mux_output_var_idx);
					}
					if (w == this->word_size-1) {
						// these clauses are different for the sign bit
						this->create_signed_shift_overflow_protection(select_input_var_idx, zero_input_sign_bit_idx, one_input_var_idx);
					}
					else {
						this->create_signed_shift_overflow_protection(select_input_var_idx, zero_input_sign_bit_idx, zero_input_var_idx);
					}
				}
				else {
					if (connect_zero_const) {
						this->create_1x1_equivalence(zero_input_var_idx, mux_output_var_idx);
						this->create_1x1_negated_implication(zero_input_var_idx, select_input_var_idx);
						this->create_1x1_negated_implication(mux_output_var_idx, select_input_var_idx);
					}
					else {
						this->create_2x1_mux_shift_disallowed(zero_input_var_idx, one_input_var_idx, select_input_var_idx, mux_output_var_idx);
					}
				}
			}
			else {
				if (connect_zero_const) {
					this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx, mux_output_var_idx);
				}
				else {
					this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx, mux_output_var_idx);
				}
			}
		}
		// sign bits before and after shifting must be identical if calculating in 2's complement
		if (this->calc_twos_complement) {
			int shift_input_sign_bit_idx;
			int shift_output_sign_bit_idx = this->shift_output_variables.at({idx, this->word_size-1});
			if (idx == 1) {
				shift_input_sign_bit_idx = this->output_value_variables.at({0, this->word_size-1});
			}
			else {
				shift_input_sign_bit_idx = this->shift_select_output_variables.at({idx, scm::left, this->word_size-1});
			}
			this->create_1x1_equivalence(shift_input_sign_bit_idx, shift_output_sign_bit_idx);
		}
	}
}

void scm::create_post_adder_shift_constraints(int idx) {
	for (auto stage = 0; stage < this->shift_word_size; stage++) {
		auto shift_width = (1 << stage);
		auto select_input_var_idx = this->input_post_adder_shift_value_variables.at({idx, stage});
		auto last_disallowed_shift_bit = shift_width-1;
		for (auto w = 0; w < this->word_size; w++) {
			auto w_prev = w + shift_width;
			auto connect_zero_const = w_prev >= this->word_size;
			int zero_input_var_idx;
			int zero_input_sign_bit_idx;
			int one_input_var_idx;
			auto mux_output_var_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage, w});
			if (stage == 0) {
				// connect shifter inputs
				// shifter input is the adder output
				zero_input_var_idx = this->adder_output_value_variables.at({idx, w});
				zero_input_sign_bit_idx = this->adder_output_value_variables.at({idx, this->word_size-1});
				if (!connect_zero_const) {
					one_input_var_idx = this->adder_output_value_variables.at({idx, w_prev});
				}
			}
			else {
				// connect output of previous stage
				zero_input_var_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage-1, w});
				zero_input_sign_bit_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage-1, this->word_size-1});
				if (!connect_zero_const) {
					one_input_var_idx = this->post_adder_shift_internal_mux_output_variables.at({idx, stage-1, w_prev});
				}
			}
			if (w <= last_disallowed_shift_bit) {
				// shifting out 1s is not allowed in these places
				if (connect_zero_const) {
					if (this->calc_twos_complement) {
						// connect the sign bit instead of a constant zero
						this->create_2x1_mux_shift_disallowed(zero_input_var_idx, zero_input_sign_bit_idx, select_input_var_idx, mux_output_var_idx);
					}
					else {
						this->create_1x1_equivalence(zero_input_var_idx, mux_output_var_idx);
						this->create_1x1_negated_implication(zero_input_var_idx, select_input_var_idx);
						this->create_1x1_negated_implication(mux_output_var_idx, select_input_var_idx);
					}
				}
				else {
					this->create_2x1_mux_shift_disallowed(zero_input_var_idx, one_input_var_idx, select_input_var_idx, mux_output_var_idx);
				}
			}
			else {
				// we can shift bits around however we like
				if (connect_zero_const) {
					if (this->calc_twos_complement) {
						// connect the sign bit instead of a constant zero
						this->create_2x1_mux(zero_input_var_idx, zero_input_sign_bit_idx, select_input_var_idx, mux_output_var_idx);
					}
					else {
						this->create_2x1_mux_zero_const(zero_input_var_idx, select_input_var_idx, mux_output_var_idx);
					}
				}
				else {
					this->create_2x1_mux(zero_input_var_idx, one_input_var_idx, select_input_var_idx, mux_output_var_idx);
				}
			}
		}
		// sign bits before and after shifting must be identical if calculating in 2's complement
		if (this->calc_twos_complement) {
			//int shift_output_sign_bit_idx = this->shift_output_variables.at({idx, this->word_size-1});
			int shift_output_sign_bit_idx = this->post_adder_shift_output_variables.at({idx, this->word_size-1});
			//int shift_input_sign_bit_idx = this->shift_select_output_variables.at({idx, scm::left, this->word_size-1});
			int shift_input_sign_bit_idx = this->adder_output_value_variables.at({idx, this->word_size-1});
			this->create_1x1_equivalence(shift_input_sign_bit_idx, shift_output_sign_bit_idx);
		}
	}
}

void scm::create_shift_select_constraints(int idx) {
	if (idx == 1) return; // this node doesn't need shift input muxs
	auto select_var = this->input_shift_select_variables.at(idx);
	for (auto &dir : this->input_directions) {
		for (auto w = 0; w < this->word_size; w++) {
			auto mux_output_var_idx = this->shift_select_output_variables.at({idx, dir, w});
			auto left_input_var_idx = this->input_select_mux_variables.at({idx, scm::left, 0, w});
			auto right_input_var_idx = this->input_select_mux_variables.at({idx, scm::right, 0, w});
			if (dir == scm::left) {
				// left mux has left input in b input and right input in a input
				this->create_2x1_mux(right_input_var_idx, left_input_var_idx, select_var, mux_output_var_idx);
			}
			else {
				// right mux has left input in '0' input and right input in '1' input
				this->create_2x1_mux(left_input_var_idx, right_input_var_idx, select_var, mux_output_var_idx);
			}
		}
	}
}

void scm::create_negate_select_constraints(int idx) {
	auto select_var_idx = this->input_negate_select_variables.at(idx);
	for (int w = 0; w < this->word_size; w++) {
		auto left_input_var_idx = this->shift_output_variables.at({idx, w});
		int right_input_var_idx;
		if (idx == 1) {
			// right input is the output of the input node with idx = 0
			right_input_var_idx = this->output_value_variables.at({0, w});
		}
		else {
			// right input is the output of the right shift select mux
			right_input_var_idx = this->shift_select_output_variables.at({idx, scm::right, w});
		}
		for (auto &dir : this->input_directions) {
			auto mux_output_var_idx = this->negate_select_output_variables.at({idx, dir, w});
			if (dir == scm::left) {
				this->create_2x1_mux(right_input_var_idx, left_input_var_idx, select_var_idx, mux_output_var_idx);
			}
			else {
				this->create_2x1_mux(left_input_var_idx, right_input_var_idx, select_var_idx, mux_output_var_idx);
			}
		}
	}
}

void scm::create_xor_constraints(int idx) {
	auto negate_var_idx = this->input_negate_value_variables.at(idx);
	for (int w = 0; w < this->word_size; w++) {
		auto input_var_idx = this->negate_select_output_variables.at({idx, scm::right, w});
		auto output_var_idx = this->xor_output_variables.at({idx, w});
		this->create_2x1_xor(negate_var_idx, input_var_idx, output_var_idx);
	}
}

void scm::create_adder_constraints(int idx) {
	for (int w = 0; w < this->word_size; w++) {
		int c_i;
		if (w == 0) {
			// carry input = input negate value
			c_i = this->input_negate_value_variables.at(idx);
		}
		else {
			// carry input = carry output of last stage
			c_i = this->adder_internal_variables.at({idx, w-1});
		}
		// build sum
		int a = this->negate_select_output_variables.at({idx, scm::left, w});
		int b = this->xor_output_variables.at({idx, w});
		int s = this->adder_output_value_variables.at({idx, w});
		this->create_add_sum(a, b, c_i, s);
		// build carry
		int c_o = this->adder_internal_variables.at({idx, w});
		this->create_add_carry(a, b, c_i, c_o);
	}
	// disallow overflows
	if (this->calc_twos_complement) {
		this->create_signed_add_overflow_protection(this->input_negate_value_variables.at(idx), this->negate_select_output_variables.at({idx, scm::left, this->word_size-1}), this->negate_select_output_variables.at({idx, scm::right, this->word_size-1}), this->output_value_variables.at({idx, this->word_size-1}));
	}
	else {
		this->create_1x1_equivalence(this->adder_internal_variables.at({idx, this->word_size-1}), this->input_negate_value_variables.at(idx));
	}
}

void scm::create_full_adder_allocation_constraints(int idx) {
	// define global FA alloc for each bit with index x
	if (this->calc_twos_complement) {
		// negative fundamentals are allowed
		for (int x = 0; x < this->word_size-1; x++) {
			int container_size = this->max_full_adders + 4;
			std::vector<std::pair<int, bool>> vars(container_size);
			// y: global FA alloc index
			for (int y = 0; y < this->max_full_adders; y++) {
				vars[y] = {this->full_adder_alloc_variables.at({idx, x, y}), false};
			}
			for (int w_1 = x; w_1 < this->word_size-1; w_1++) {
				// for (a << s) - b
				// subtract bit
				vars[this->max_full_adders] = {this->input_negate_value_variables.at(idx), true};
				// negate select bit
				vars[this->max_full_adders + 1] = {this->input_negate_select_variables.at(idx), true};
				// adder output value bit = 1 and adder output value sign bit = 0
				vars[this->max_full_adders + 2] = {this->adder_output_value_variables.at({idx, w_1}), true};
				vars[this->max_full_adders + 3] = {this->adder_output_value_variables.at({idx, this->word_size-1}), false};
				// add to solver
				this->create_arbitrary_clause(vars);
				// adder output value bit = 0 and adder output value sign bit = 1
				vars[this->max_full_adders + 2] = {this->adder_output_value_variables.at({idx, w_1}), false};
				vars[this->max_full_adders + 3] = {this->adder_output_value_variables.at({idx, this->word_size-1}), true};
				// add to solver
				this->create_arbitrary_clause(vars);
				// for a + (b << s) and a - (b << s)
				for (int w_2 = 0; w_2 <= x; w_2++) {
					// subtract bit
					vars[this->max_full_adders] = {this->input_negate_value_variables.at(idx), false};
					// shifted value bit
					vars[this->max_full_adders + 1] = {this->shift_output_variables.at({idx, w_2}), true};
					// adder output value bit = 1 and adder output value sign bit = 0
					vars[this->max_full_adders + 2] = {this->adder_output_value_variables.at({idx, w_1}), true};
					vars[this->max_full_adders + 3] = {this->adder_output_value_variables.at({idx, this->word_size-1}), false};
					// add to solver
					this->create_arbitrary_clause(vars);
					// adder output value bit = 0 and adder output value sign bit = 1
					vars[this->max_full_adders + 2] = {this->adder_output_value_variables.at({idx, w_1}), false};
					vars[this->max_full_adders + 3] = {this->adder_output_value_variables.at({idx, this->word_size-1}), true};
					// add to solver
					this->create_arbitrary_clause(vars);
					// swap subtract with select bit
					vars[this->max_full_adders] = {this->input_negate_select_variables.at(idx), false};
					// add to solver
					this->create_arbitrary_clause(vars);
					// adder output value bit = 1 and adder output value sign bit = 0
					vars[this->max_full_adders + 2] = {this->adder_output_value_variables.at({idx, w_1}), true};
					vars[this->max_full_adders + 3] = {this->adder_output_value_variables.at({idx, this->word_size-1}), false};
					// add to solver
					this->create_arbitrary_clause(vars);
				}
			}
		}
	}
	else {
		// only positive fundamentals are allowed
		for (int x = 0; x < this->word_size; x++) {
			int container_size = this->max_full_adders + 3;
			std::vector<std::pair<int, bool>> vars(container_size);
			// y: global FA alloc index
			for (int y = 0; y < this->max_full_adders; y++) {
				vars[y] = {this->full_adder_alloc_variables.at({idx, x, y}), false};
			}
			for (int w_1 = x; w_1 < this->word_size; w_1++) {
				// for (a << s) - b
				// subtract bit
				vars[this->max_full_adders] = {this->input_negate_value_variables.at(idx), true};
				// adder output value bit
				vars[this->max_full_adders+1] = {this->adder_output_value_variables.at({idx, w_1}), true};
				// negate select bit
				vars[this->max_full_adders+2] = {this->input_negate_select_variables.at(idx), true};
				// add to solver
				this->create_arbitrary_clause(vars);
				// for a + (b << s) and a - (b << s)
				for (int w_2 = 0; w_2 <= x; w_2++) {
					// subtract bit
					vars[this->max_full_adders] = {this->input_negate_value_variables.at(idx), false};
					// adder output value bit
					vars[this->max_full_adders+1] = {this->adder_output_value_variables.at({idx, w_1}), true};
					// shifted value bit
					vars[this->max_full_adders+2] = {this->shift_output_variables.at({idx, w_2}), true};
					// add to solver
					this->create_arbitrary_clause(vars);
					// swap subtract with select bit
					vars[this->max_full_adders] = {this->input_negate_select_variables.at(idx), false};
					// add to solver
					this->create_arbitrary_clause(vars);
				}
			}
		}
	}
}

void scm::create_full_adder_overlap_constraints(int idx_1) {
	if (this->max_full_adders <= 0) return;
	auto w_lim = (this->calc_twos_complement?this->word_size-1:this->word_size);
	for (int w_1 = 0; w_1 < w_lim; w_1++) {
		for (int idx_2 = idx_1; idx_2 <= this->num_adders; idx_2++) {
			for (int w_2 = 0; w_2 < w_lim; w_2++) {
				if (idx_1 == idx_2 and w_1 == w_2) continue; // no overlap with itself
				for (int x = 0; x < this->max_full_adders; x++) {
					this->create_1x1_negated_implication(this->full_adder_alloc_variables.at({idx_1, w_1, x}), this->full_adder_alloc_variables.at({idx_2, w_2, x}));
				}
			}
		}
	}
}

void scm::create_input_select_limitation_constraints(int idx) {
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

void scm::create_shift_limitation_constraints(int idx) {
	int max_representable_shift = (1 << this->shift_word_size) - 1;
	std::vector<int> x(this->shift_word_size);
	for (int w = 0; w < this->shift_word_size; w++) {
		x[w] = this->input_shift_value_variables.at({idx, w});
	}
	for (int forbidden_number = max_representable_shift; forbidden_number > this->max_shift; forbidden_number--) {
		this->forbid_number(x, forbidden_number);
	}
}

void scm::create_post_adder_shift_limitation_constraints(int idx) {
	int max_representable_shift = (1 << this->shift_word_size) - 1;
	std::vector<int> x(this->shift_word_size);
	for (int w = 0; w < this->shift_word_size; w++) {
		x[w] = this->input_post_adder_shift_value_variables.at({idx, w});
	}
	for (int forbidden_number = max_representable_shift; forbidden_number > this->max_shift; forbidden_number--) {
		this->forbid_number(x, forbidden_number);
	}
}

void scm::get_solution_from_backend() {
	// clear containers
	this->input_select.clear();
	this->input_select_mux_output.clear();
	this->shift_input_select.clear();
	this->shift_value.clear();
	this->negate_select.clear();
	this->subtract.clear();
	this->post_adder_shift_value.clear();
	this->add_result_values.clear();
	this->output_values.clear();
	// get solution
	for (int idx = 0; idx <= this->num_adders; idx++) {
		// output_values
		for (int w = 0; w < this->word_size; w++) {
			this->output_values[idx] += (this->get_result_value(this->output_value_variables.at({idx, w})) << w);
		}
		if (this->calc_twos_complement) this->output_values[idx] = sign_extend(this->output_values[idx], this->word_size);
		if (idx > 0) {
			if (idx > 1) {
				// input_select
				for (auto &dir : this->input_directions) {
					auto input_select_width = this->ceil_log2(idx);
					for (auto w = 0; w < input_select_width; w++) {
						this->input_select[{idx, dir}] += (this->get_result_value(this->input_select_selection_variables[{idx, dir, w}]) << w);
					}
				}
				// shift_input_select
				this->shift_input_select[idx] = this->get_result_value(this->input_shift_select_variables[idx]);
			}
			// shift_value
			for (auto w = 0; w < this->shift_word_size; w++) {
				this->shift_value[idx] += (this->get_result_value(this->input_shift_value_variables[{idx, w}]) << w);
			}
			// negate_select
			this->negate_select[idx] = this->get_result_value(this->input_negate_select_variables[idx]);
			// subtract
			this->subtract[idx] = this->get_result_value(this->input_negate_value_variables[idx]);
			// output shift
			if (this->enable_node_output_shift) {
				for (auto w = 0; w < this->shift_word_size; w++) {
					this->post_adder_shift_value[idx] += (this->get_result_value(this->input_post_adder_shift_value_variables[{idx, w}]) << w);
				}
			}
			// add result
			for (auto w = 0; w < this->word_size; w++) {
				this->add_result_values[idx] += (this->get_result_value(this->adder_output_value_variables[{idx, w}]) << w);
			}
			if (this->calc_twos_complement) this->add_result_values[idx] = sign_extend(this->add_result_values[idx], this->word_size);
		}
	}
}

void scm::print_solution() {
	if (this->found_solution) {
		std::cout << "Solution for constants" << std::endl;
		for (auto &c : this->C) {
			std::cout << "  C = " << c << std::endl;
		}
		std::cout << "#adders = " << this->num_adders << ", word size = " << this->word_size << std::endl;
		std::cout << "  node #0 = " << (this->calc_twos_complement?sign_extend(this->output_values[0], this->word_size):this->output_values[0]) << std::endl;
		for (auto idx = 1; idx <= this->num_adders; idx++) {
			std::cout << "  node #" << idx << " = " << (this->calc_twos_complement?sign_extend((int64_t)this->output_values[idx], this->word_size):this->output_values[idx]) << std::endl;
			std::cout << "    left input: node " << this->input_select[{idx, scm::left}] << std::endl;
			std::cout << "    right input: node " << this->input_select[{idx, scm::right}] << std::endl;
			std::cout << "    shift input select: " << this->shift_input_select[idx] << (this->shift_input_select[idx]==1?" (left)":" (right)") << std::endl;
			std::cout << "    shift value: " << this->shift_value[idx] << std::endl;
			std::cout << "    negate select: " << this->negate_select[idx] << (this->negate_select[idx]==1?" (non-shifted)":" (shifted)") << std::endl;
			std::cout << "    subtract: " << this->subtract[idx] << std::endl;
			if (this->enable_node_output_shift) {
				std::cout << "    post adder right shift value: " << this->post_adder_shift_value[idx] << std::endl;
			}
		}
		std::cerr << "Adder graph: " << this->get_adder_graph_description() << std::endl;
	}
	else {
		for (auto &c : this->C) {
			std::cout << "Failed to find solution for constants" << std::endl;
			std::cout << "  C = " << c << std::endl;
		}
	}
}

bool scm::solution_is_valid() {
	bool valid = true;
	for (int idx = 1; idx <= this->num_adders; idx++) {
		// verify node inputs
		int64_t input_node_idx_l = 0;
		int64_t input_node_idx_r = 0;
		int64_t actual_input_value_l = 1;
		int64_t actual_input_value_r = 1;
		if (idx > 1) {
			for (auto &dir : this->input_directions) {
				for (auto w = 0; w < this->word_size; w++) {
					this->input_select_mux_output[{idx, dir}] += (
						this->get_result_value(this->input_select_mux_variables[{idx, dir, 0, w}]) << w);
				}
			}
			input_node_idx_l = this->input_select[{idx, scm::left}];
			input_node_idx_r = this->input_select[{idx, scm::right}];
			actual_input_value_l = this->input_select_mux_output[{idx, scm::left}];
			actual_input_value_r = this->input_select_mux_output[{idx, scm::right}];
			if (this->calc_twos_complement) actual_input_value_l = sign_extend(actual_input_value_l, this->word_size);
			if (this->calc_twos_complement) actual_input_value_r = sign_extend(actual_input_value_r, this->word_size);
		}
		int64_t left_input_value = this->output_values[input_node_idx_l];
		int64_t right_input_value = this->output_values[input_node_idx_r];
		if (!this->quiet) {
			std::cout << "node #" << idx << " left input" << std::endl;
			std::cout << "  input select = " << input_node_idx_l << std::endl;
			std::cout << "  value = " << actual_input_value_l << std::endl;
		}
		if (left_input_value != actual_input_value_l) {
			std::cout << "node #" << idx << " has invalid left input" << std::endl;
			std::cout << "  input select = " << input_node_idx_l << std::endl;
			std::cout << "  expected value " << left_input_value << " but got " << actual_input_value_l << std::endl;
			auto num_muxs = (1 << this->ceil_log2(idx))-1;
			for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
				int64_t mux_output = 0;
				for (auto w = 0; w < this->word_size; w++) {
					mux_output += (this->get_result_value(this->input_select_mux_variables[{idx, scm::left, mux_idx, w}]) << w);
				}
				std::cout << "    mux #" << mux_idx << " output: " << mux_output << std::endl;
			}
			valid = false;
		}
		if (!this->quiet) {
			std::cout << "node #" << idx << " right input" << std::endl;
			std::cout << "  input select = " << input_node_idx_r << std::endl;
			std::cout << "  value = " << actual_input_value_r << std::endl;
		}
		if (right_input_value != actual_input_value_r) {
			std::cout << "node #" << idx << " has invalid right input" << std::endl;
			std::cout << "  input select = " << input_node_idx_r << std::endl;
			std::cout << "  expected value " << right_input_value << " but got " << actual_input_value_r << std::endl;
			int64_t num_muxs = (1 << this->ceil_log2(idx))-1;
			for (int mux_idx = 0; mux_idx < num_muxs; mux_idx++) {
				int64_t mux_output = 0;
				for (auto w = 0; w < this->word_size; w++) {
					mux_output += (this->get_result_value(this->input_select_mux_variables[{idx, scm::right, mux_idx, w}]) << w);
				}
				std::cout << "    mux #" << mux_idx << " output: " << mux_output << std::endl;
			}
			valid = false;
		}
		// verify shift mux outputs
		int64_t shift_mux_output_l = left_input_value;
		int64_t shift_mux_output_r = right_input_value;
		if (idx > 1) {
			if (this->get_result_value(this->input_shift_select_variables[idx]) == 0) {
				shift_mux_output_l = right_input_value;
				shift_mux_output_r = left_input_value;
			}
			std::map<scm::input_direction, int> actual_shift_mux_output;
			for (auto &dir : this->input_directions) {
				for (auto w = 0; w < this->word_size; w++) {
					actual_shift_mux_output[dir] += (this->get_result_value(this->shift_select_output_variables[{idx, dir, w}]) << w);
				}
			}
			if (this->calc_twos_complement) actual_shift_mux_output[scm::left] = sign_extend(actual_shift_mux_output[scm::left], this->word_size);
			if (this->calc_twos_complement) actual_shift_mux_output[scm::right] = sign_extend(actual_shift_mux_output[scm::right], this->word_size);
			if (!this->quiet) {
				std::cout << "node #" << idx << " left shift input select mux output" << std::endl;
				std::cout << "  select = " << this->get_result_value(this->input_shift_select_variables[idx]) << std::endl;
				std::cout << "  output value = " <<  actual_shift_mux_output[scm::left] << std::endl;
			}
			if (shift_mux_output_l != actual_shift_mux_output[scm::left]) {
				std::cout << "node #" << idx << " has invalid left shift input select mux output" << std::endl;
				std::cout << "  select = " << this->get_result_value(this->input_shift_select_variables[idx]) << std::endl;
				std::cout << "  actual value = " <<  actual_shift_mux_output[scm::left] << std::endl;
				std::cout << "  expected value = " << shift_mux_output_l << std::endl;
				valid = false;
			}
			if (!this->quiet) {
				std::cout << "node #" << idx << " right shift input select mux output" << std::endl;
				std::cout << "  select = " << this->get_result_value(this->input_shift_select_variables[idx]) << std::endl;
				std::cout << "  output value = " <<  actual_shift_mux_output[scm::right] << std::endl;
			}
			if (shift_mux_output_r != actual_shift_mux_output[scm::right]) {
				std::cout << "node #" << idx << " has invalid left shift input select mux output" << std::endl;
				std::cout << "  select = " << this->get_result_value(this->input_shift_select_variables[idx]) << std::endl;
				std::cout << "  actual value = " <<  actual_shift_mux_output[scm::right] << std::endl;
				std::cout << "  expected value = " << shift_mux_output_r << std::endl;
				valid = false;
			}
		}
		// verify shifter output
		int64_t expected_shift_output = (((int64_t)shift_mux_output_l) << this->shift_value[idx]);// % (int64_t)(1 << this->word_size);
		if (this->calc_twos_complement) expected_shift_output = sign_extend(expected_shift_output, this->word_size);
		int64_t actual_shift_output = 0;
		for (int w = 0; w < this->word_size; w++) {
			actual_shift_output += (this->get_result_value(this->shift_output_variables[{idx, w}]) << w);
		}
		if (this->calc_twos_complement) actual_shift_output = sign_extend(actual_shift_output, this->word_size);
		if (!this->quiet) {
			std::cout << "node #" << idx << " shift output" << std::endl;
			std::cout << "  input value = " << shift_mux_output_l << std::endl;
			std::cout << "  shift value = " << this->shift_value[idx] << std::endl;
			std::cout << "  output value = " << actual_shift_output << std::endl;
		}
		if (expected_shift_output != actual_shift_output) {
			std::cout << "node #" << idx << " has invalid shift output" << std::endl;
			std::cout << "  input value = " << shift_mux_output_l << std::endl;
			std::cout << "  shift value = " << this->shift_value[idx] << std::endl;
			std::cout << "  expected output value = " << expected_shift_output << std::endl;
			std::cout << "  actual output value = " << actual_shift_output << std::endl;
			valid = false;
		}
		// verify negate mux outputs
		int64_t negate_mux_output_l = actual_shift_output;
		int64_t negate_mux_output_r = shift_mux_output_r;
		if (this->get_result_value(this->input_negate_select_variables[idx]) == 0) {
			negate_mux_output_l = shift_mux_output_r;
			negate_mux_output_r = actual_shift_output;
		}
		std::map<scm::input_direction, int> actual_negate_mux_output;
		for (auto &dir : this->input_directions) {
			for (auto w = 0; w < this->word_size; w++) {
				actual_negate_mux_output[dir] += (this->get_result_value(this->negate_select_output_variables[{idx, dir, w}]) << w);
			}
		}
		if (this->calc_twos_complement) actual_negate_mux_output[scm::left] = sign_extend(actual_negate_mux_output[scm::left], this->word_size);
		if (this->calc_twos_complement) actual_negate_mux_output[scm::right] = sign_extend(actual_negate_mux_output[scm::right], this->word_size);
		if (!this->quiet) {
			std::cout << "node #" << idx << " left negate select mux output" << std::endl;
			std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx]) << std::endl;
			std::cout << "  output value = " <<  actual_negate_mux_output[scm::left] << std::endl;
		}
		if (negate_mux_output_l != actual_negate_mux_output[scm::left]) {
			std::cout << "node #" << idx << " has invalid left negate select mux output" << std::endl;
			std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx]) << std::endl;
			std::cout << "  actual value = " <<  actual_negate_mux_output[scm::left] << std::endl;
			std::cout << "  expected value = " << negate_mux_output_l << std::endl;
			valid = false;
		}
		if (!this->quiet) {
			std::cout << "node #" << idx << " right negate select mux output" << std::endl;
			std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx]) << std::endl;
			std::cout << "  output value = " <<  actual_negate_mux_output[scm::right] << std::endl;
		}
		if (negate_mux_output_r != actual_negate_mux_output[scm::right]) {
			std::cout << "node #" << idx << " has invalid right negate select mux output" << std::endl;
			std::cout << "  select = " << this->get_result_value(this->input_negate_select_variables[idx]) << std::endl;
			std::cout << "  actual value = " <<  actual_negate_mux_output[scm::right] << std::endl;
			std::cout << "  expected value = " << negate_mux_output_r << std::endl;
			valid = false;
		}
		// verify xor output
		int64_t sub = this->get_result_value(this->input_negate_value_variables[idx]);
		int64_t expected_xor_output = sub==1?(~negate_mux_output_r) & ((((int64_t)1) << this->word_size) - 1):negate_mux_output_r;
		if (this->calc_twos_complement) expected_xor_output = sign_extend(expected_xor_output, this->word_size);
		int64_t actual_xor_output = 0;
		for (int w = 0; w < this->word_size; w++) {
			actual_xor_output += (this->get_result_value(this->xor_output_variables[{idx, w}]) << w);
		}
		if (this->calc_twos_complement) actual_xor_output = sign_extend(actual_xor_output, this->word_size);
		if (!this->quiet) {
			std::cout << "node #" << idx << " xor output" << std::endl;
			std::cout << "  sub = " << sub << std::endl;
			std::cout << "  input value = " << negate_mux_output_r << std::endl;
			std::cout << "  output value = " <<  actual_xor_output << std::endl;
		}
		if (expected_xor_output != actual_xor_output) {
			std::cout << "node #" << idx << " has invalid xor output" << std::endl;
			std::cout << "  sub = " << sub << std::endl;
			std::cout << "  input value = " << negate_mux_output_r << std::endl;
			std::cout << "  actual output value = " <<  actual_xor_output << std::endl;
			std::cout << "  expected output value = " << expected_xor_output << std::endl;
			valid = false;
		}
		// verify adder output
		int64_t expected_adder_output = (sub == 1) ? (negate_mux_output_l - negate_mux_output_r) : (negate_mux_output_l + negate_mux_output_r);
		if (this->calc_twos_complement) expected_adder_output = sign_extend(expected_adder_output, this->word_size);
		int64_t actual_adder_output = 0;
		for (int w = 0; w < this->word_size; w++) {
			actual_adder_output += (this->get_result_value(this->adder_output_value_variables[{idx, w}]) << w);
		}
		if (this->calc_twos_complement) actual_adder_output = sign_extend(actual_adder_output, this->word_size);
		if (!this->quiet) {
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
		if (this->enable_node_output_shift) {
			// verify post adder shift output
			int64_t expected_post_adder_shift_output = actual_adder_output >> this->post_adder_shift_value.at(idx);
			if (this->calc_twos_complement) expected_post_adder_shift_output = sign_extend(expected_post_adder_shift_output, this->word_size);
			int64_t actual_post_adder_shift_output = 0;
			for (int w = 0; w < this->word_size; w++) {
				actual_post_adder_shift_output += (this->get_result_value(this->post_adder_shift_output_variables[{idx, w}]) << w);
			}
			if (this->calc_twos_complement) actual_post_adder_shift_output = sign_extend(actual_post_adder_shift_output, this->word_size);
			if (!this->quiet) {
				std::cout << "node #" << idx << " post adder shift output value" << std::endl;
				std::cout << "  shift value = " << this->post_adder_shift_value.at(idx) << std::endl;
				std::cout << "  input value = " << actual_adder_output << std::endl;
				std::cout << "  actual output value = " << actual_post_adder_shift_output << std::endl;
			}
			if (expected_post_adder_shift_output != actual_post_adder_shift_output)  {
				std::cout << "node #" << idx << " has invalid post adder shift output value" << std::endl;
				std::cout << "  shift value = " << this->post_adder_shift_value.at(idx) << std::endl;
				std::cout << "  input value = " << actual_adder_output << std::endl;
				std::cout << "  expected output value = " << expected_post_adder_shift_output << std::endl;
				std::cout << "  actual output value = " << actual_post_adder_shift_output << std::endl;
				valid = false;
			}
		}
	}
	return valid;
}

void scm::create_cnf_file() {
	std::ofstream f;
	std::stringstream constants;
	for (int i=0; i<this->C.size(); i++) {
		if (i != 0) constants << "_";
		constants << this->C[i];
	}
	std::string filename;
	if (this->max_full_adders >= 0) {
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

void scm::create_mcm_output_constraints() {
	for (auto &c : this->C) {
		std::vector<int> or_me;
		for (int idx = 1; idx <= this->num_adders; idx++) {
			or_me.emplace_back(this->mcm_output_variables[{idx, c}]);
			for (int w = 0; w < this->word_size; w++) {
				if (((c >> w) & 1) == 1) {
					this->create_1x1_implication(this->mcm_output_variables[{idx, c}], this->output_value_variables[{idx, w}]);
				}
				else {
					this->create_1x1_negated_implication(this->mcm_output_variables[{idx, c}], this->output_value_variables[{idx, w}]);
				}
			}
		}
		if (this->calc_twos_complement and this->sign_inversion_allowed[c]) {
			// also allow the solver to choose -c instead of c if it's easier to implement
			for (int idx = 1; idx <= this->num_adders; idx++) {
				or_me.emplace_back(this->mcm_output_variables[{idx, -c}]);
				for (int w = 0; w < this->word_size; w++) {
					if ((((-c) >> w) & 1) == 1) {
						this->create_1x1_implication(this->mcm_output_variables[{idx, -c}], this->output_value_variables[{idx, w}]);
					}
					else {
						this->create_1x1_negated_implication(this->mcm_output_variables[{idx, -c}], this->output_value_variables[{idx, w}]);
					}
				}
			}
		}
		this->create_or(or_me);
	}
}

void scm::create_mcm_output_variables(int idx) {
	for (auto &c : this->C) {
		this->mcm_output_variables[{idx, c}] = ++this->variable_counter;
		this->create_new_variable(this->variable_counter);
		if (this->calc_twos_complement and this->sign_inversion_allowed[c]) {
			this->mcm_output_variables[{idx, -c}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

void scm::create_full_adder_alloc_variables(int idx) {
	auto w_lim = (this->calc_twos_complement?this->word_size-1:this->word_size);
	for (int w = 0; w < w_lim; w++) {
		for (int x = 0; x < this->max_full_adders; x++) {
			this->full_adder_alloc_variables[{idx, w, x}] = ++this->variable_counter;
			this->create_new_variable(this->variable_counter);
		}
	}
}

int64_t scm::sign_extend(int64_t x, int w) {
	auto sign_bit = (x >> (w-1)) & 1;
	if (sign_bit == 0) return x; // x >= 0 -> no conversion needed
	auto mask = (1 << w) - 1;
	mask = ~mask;
	x = x | mask;
	return x;
}

std::string scm::get_adder_graph_description() {
	std::stringstream s;
	if (!this->found_solution) return s.str();
	s << "{";
	std::map<int, int> stage;
	stage[0] = 0;
	for (int idx = 1; idx <= this->num_adders; idx++) {
		// insert comma if needed
		if (idx > 1) s << ",";
		// get left and right inputs and their shift
		int left_idx;
		int right_idx;
		int left_input;
		int right_input;
		int left_shift;
		int right_shift;
		if (this->shift_input_select.at(idx) == 1) {
			// do not swap inputs
			left_idx = this->input_select.at({idx, scm::left});
			right_idx = this->input_select.at({idx, scm::right});
		}
		else {
			// swap inputs
			left_idx = this->input_select.at({idx, scm::right});
			right_idx = this->input_select.at({idx, scm::left});
		}
		left_input = this->output_values.at(left_idx);
		right_input = this->output_values.at(right_idx);
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
		int left_stage = stage.at(left_idx);
		int right_stage = stage.at(right_idx);
		int current_stage = std::max(left_stage, right_stage)+1;
		stage[idx] = current_stage;
		// basic node info
		s << "{'A',[" << this->output_values.at(idx) << "]," << current_stage;
		if (this->enable_node_output_shift) {
			s << "," << this->post_adder_shift_value.at(idx);
		}
		// left input
		s << ",[" << left_input << "]," << left_stage << "," << left_shift;
		// right input
		s << ",[" << right_input << "]," << right_stage << "," << right_shift;
		// close bracket
		s << "}";
	}
	s << "}";
	return s.str();
}

void scm::set_min_add(int new_min_add) {
	this->num_adders = std::max(this->num_adders, new_min_add-1);
	this->num_adders = std::max(this->num_adders, 0);
}

void scm::also_minimize_full_adders() {
	this->minimize_full_adders = true;
}

void scm::allow_node_output_shift() {
	this->enable_node_output_shift = true;
}

std::pair<int, int> scm::solution_is_optimal() {
	return {this->num_add_opt, this->num_FA_opt};
}

void scm::ignore_sign(bool only_apply_to_negative_coefficients) {
	for (auto &c : this->C) {
		if (only_apply_to_negative_coefficients) {
			// only allow sign inversion for negative coefficients
			if (this->negative_coeff_requested[c]) {
				this->sign_inversion_allowed[c] = true;
			}
		}
		else {
			// only the solver to choose sign for all coefficients
			this->sign_inversion_allowed[c] = true;
		}
	}
}
