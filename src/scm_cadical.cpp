//
// Created by nfiege on 9/26/22.
//

#include "scm_cadical.h"

#ifdef USE_CADICAL

#include <iostream>

scm_cadical::scm_cadical(int C, int timeout, bool quiet, int word_size) : scm(C, timeout, quiet, word_size) {}

void scm_cadical::reset_backend() {
	// create new solver
	this->solver = std::make_unique<CaDiCaL::Solver>();
	// create and attach new terminator
	this->terminator = cadical_terminator(this->timeout);
	this->solver->connect_terminator(&this->terminator);
}

void scm_cadical::create_1x1_equivalence(int x, int y) {
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

void scm_cadical::create_2x1_xor(int a, int b, int y) {
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
	if (val == 1)
		this->solver->add(x);
	else
		this->solver->add(-x);
	this->solver->add(0);
}

void scm_cadical::forbid_number(const std::vector<int> &x, int val) {
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