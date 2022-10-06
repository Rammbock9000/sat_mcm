//
// Created by nfiege on 9/26/22.
//

#include "scm_cadical.h"

#ifdef USE_CADICAL

#include <iostream>

scm_cadical::scm_cadical(const std::vector<int> &C, int timeout, bool quiet, int word_size) : scm(C, timeout, quiet, word_size, 1) {}

void scm_cadical::reset_backend() {
	scm::reset_backend();
	// create new solver
	this->solver = std::make_unique<CaDiCaL::Solver>();
	// create and attach new terminator
	this->terminator = cadical_terminator(this->timeout);
	this->solver->connect_terminator(&this->terminator);
}

void scm_cadical::create_1x1_equivalence(int x, int y) {
	scm::create_1x1_equivalence(x, y);
	// 1)
	this->solver->add(-x);
	this->solver->add(y);
	this->solver->add(0);
	// 2)
	this->solver->add(x);
	this->solver->add(-y);
	this->solver->add(0);
}

void scm_cadical::create_2x1_mux(int a, int b, int s, int y) {
	scm::create_2x1_mux(a, b, s, y);
	// 1)
	this->solver->add(-a);
	this->solver->add(s);
	this->solver->add(y);
	this->solver->add(0);
	// 2)
	this->solver->add(-b);
	this->solver->add(-s);
	this->solver->add(y);
	this->solver->add(0);
	// 3)
	this->solver->add(b);
	this->solver->add(-s);
	this->solver->add(-y);
	this->solver->add(0);
	// 4)
	this->solver->add(a);
	this->solver->add(s);
	this->solver->add(-y);
	this->solver->add(0);
	// 5)
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(y);
	this->solver->add(0);
	// 6)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(-y);
	this->solver->add(0);
}

void scm_cadical::create_2x1_mux_shift_disallowed(int a, int b, int s, int y) {
	scm::create_2x1_mux_shift_disallowed(a, b, s, y);
	// 1)
	this->solver->add(-a);
	this->solver->add(y);
	this->solver->add(0);
	// 2)
	this->solver->add(-b);
	this->solver->add(-s);
	this->solver->add(y);
	this->solver->add(0);
	// 3)
	this->solver->add(b);
	this->solver->add(-s);
	this->solver->add(-y);
	this->solver->add(0);
	// 4)
	this->solver->add(a);
	this->solver->add(s);
	this->solver->add(-y);
	this->solver->add(0);
	// 5)
	this->solver->add(-a);
	this->solver->add(-s);
	this->solver->add(0);
	// 6)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(-y);
	this->solver->add(0);
}

void scm_cadical::create_2x1_mux_zero_const(int a, int s, int y) {
	scm::create_2x1_mux_zero_const(a, s, y);
	// 1)
	this->solver->add(-s);
	this->solver->add(-y);
	this->solver->add(0);
	// 2)
	this->solver->add(a);
	this->solver->add(-y);
	this->solver->add(0);
	// 3)
	this->solver->add(-a);
	this->solver->add(s);
	this->solver->add(y);
	this->solver->add(0);
}

void scm_cadical::create_2x1_xor(int a, int b, int y) {
	scm::create_2x1_xor(a, b, y);
	// 1)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(-y);
	this->solver->add(0);
	// 2)
	this->solver->add(a);
	this->solver->add(-b);
	this->solver->add(y);
	this->solver->add(0);
	// 3)
	this->solver->add(-a);
	this->solver->add(b);
	this->solver->add(y);
	this->solver->add(0);
	// 4)
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(-y);
	this->solver->add(0);
}

void scm_cadical::create_add_sum(int a, int b, int c_i, int s) {
	scm::create_add_sum(a, b, c_i, s);
	// 1)
	this->solver->add(a);
	this->solver->add(-b);
	this->solver->add(c_i);
	this->solver->add(s);
	this->solver->add(0);
	// 2)
	this->solver->add(-a);
	this->solver->add(b);
	this->solver->add(c_i);
	this->solver->add(s);
	this->solver->add(0);
	// 3)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(c_i);
	this->solver->add(-s);
	this->solver->add(0);
	// 4)
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(c_i);
	this->solver->add(-s);
	this->solver->add(0);
	// 5)
	this->solver->add(a);
	this->solver->add(-b);
	this->solver->add(-c_i);
	this->solver->add(-s);
	this->solver->add(0);
	// 6)
	this->solver->add(-a);
	this->solver->add(b);
	this->solver->add(-c_i);
	this->solver->add(-s);
	this->solver->add(0);
	// 7)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(-c_i);
	this->solver->add(s);
	this->solver->add(0);
	// 8)
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(-c_i);
	this->solver->add(s);
	this->solver->add(0);
}

void scm_cadical::create_add_carry(int a, int b, int c_i, int c_o) {
	scm::create_add_carry(a, b, c_i, c_o);
	// 1)
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(c_o);
	this->solver->add(0);
	// 2)
	this->solver->add(a);
	this->solver->add(c_i);
	this->solver->add(-c_o);
	this->solver->add(0);
	// 3)
	this->solver->add(b);
	this->solver->add(c_i);
	this->solver->add(-c_o);
	this->solver->add(0);
	// 4)
	this->solver->add(a);
	this->solver->add(b);
	this->solver->add(-c_o);
	this->solver->add(0);
	// 5)
	this->solver->add(-b);
	this->solver->add(-c_i);
	this->solver->add(c_o);
	this->solver->add(0);
	// 6)
	this->solver->add(-a);
	this->solver->add(-c_i);
	this->solver->add(c_o);
	this->solver->add(0);
}

void scm_cadical::force_bit(int x, int val) {
	scm::force_bit(x, val);
	if (val == 1)
		this->solver->add(x);
	else
		this->solver->add(-x);
	this->solver->add(0);
}

void scm_cadical::forbid_number(const std::vector<int> &x, int val) {
	scm::forbid_number(x, val);
	auto num_bits = (int)x.size();
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			this->solver->add(-x[i]);
		}
		else {
			this->solver->add(x[i]);
		}
	}
	this->solver->add(0);
}

void scm_cadical::force_number(const std::vector<int> &x, int val) {
	scm::force_number(x, val);
	auto num_bits = (int)x.size();
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			this->solver->add(x[i]);
		}
		else {
			this->solver->add(-x[i]);
		}
		this->solver->add(0);
	}
}

std::pair<bool, bool> scm_cadical::check() {
	auto stat = this->solver->solve();
	auto sat = stat == CADICAL_SAT;
	auto unsat = stat == CADICAL_UNSAT;
	auto to = !sat and !unsat;
	return {sat, to};
}

int scm_cadical::get_result_value(int var_idx) {
	return this->solver->val(var_idx) > 0 ? 1 : 0;
}

void scm_cadical::create_1x1_implication(int a, int b) {
	scm::create_1x1_implication(a, b);
	this->solver->add(-a);
	this->solver->add(b);
	this->solver->add(0);
}

void scm_cadical::create_1x1_negated_implication(int a, int b) {
	scm::create_1x1_negated_implication(a, b);
	this->solver->add(-a);
	this->solver->add(-b);
	this->solver->add(0);
}

void scm_cadical::create_or(std::vector<int> &x) {
	scm::create_or(x);
	for (auto &i : x) {
		this->solver->add(i);
	}
	this->solver->add(0);
}

void scm_cadical::create_signed_add_overflow_protection(int sub, int s_a, int s_b, int s_y) {
	scm::create_signed_add_overflow_protection(sub, s_a, s_b, s_y);
	// 1)
	this->solver->add(sub);
	this->solver->add(s_a);
	this->solver->add(s_b);
	this->solver->add(-s_y);
	this->solver->add(0);
	// 2)
	this->solver->add(sub);
	this->solver->add(-s_a);
	this->solver->add(-s_b);
	this->solver->add(s_y);
	this->solver->add(0);
	// 3)
	this->solver->add(-sub);
	this->solver->add(s_a);
	this->solver->add(-s_b);
	this->solver->add(-s_y);
	this->solver->add(0);
	// 4)
	this->solver->add(-sub);
	this->solver->add(-s_a);
	this->solver->add(s_b);
	this->solver->add(s_y);
	this->solver->add(0);
}

void scm_cadical::create_signed_shift_overflow_protection(int sel, int s_a, int a) {
	scm::create_signed_shift_overflow_protection(sel, s_a, a);
	// 1)
	this->solver->add(-sel);
	this->solver->add(-s_a);
	this->solver->add(a);
	this->solver->add(0);
	// 2)
	this->solver->add(-sel);
	this->solver->add(s_a);
	this->solver->add(-a);
	this->solver->add(0);
}

cadical_terminator::cadical_terminator(double timeout) : max_time(timeout), timer_start(std::chrono::steady_clock::now()) {}

bool cadical_terminator::terminate() {
	return this->get_elapsed_time() >= this->max_time;
}

void cadical_terminator::reset(double new_timeout) {
	this->timer_start = std::chrono::steady_clock::now();
	this->max_time = new_timeout;
}

double cadical_terminator::get_elapsed_time() const {
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - this->timer_start).count() / 1000.0;
}

#endif //USE_CADICAL