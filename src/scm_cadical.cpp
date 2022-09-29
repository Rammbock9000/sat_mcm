//
// Created by nfiege on 9/26/22.
//

#include "scm_cadical.h"
#include <iostream>

scm_cadical::scm_cadical(int C, int timeout, bool quiet) : scm(C, timeout, quiet) {}

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

void scm_cadical::create_2x1_xor(int a, int b, int y) {
	/*!
	 * force y = a XOR b
	 * clauses are:
	 *   1)  a  b -y
	 *   2)  a -b  y
	 *   3) -a  b  y
	 *   4) -a -b -y
	 */
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
	/*!
	 * force s to be the sum output of a full adder
	 * clauses are:
	 *   1)  a -b  c_i  s
	 *   2) -a  b  c_i  s
	 *   3)  a  b  c_i -s
	 *   4) -a -b  c_i -s
	 *   5)  a -b -c_i -s
	 *   6) -a  b -c_i -s
	 *   7)  a  b -c_i  s
	 *   8) -a -b -c_i  s
	 */
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
	/*!
	 * force c_o to be the carry output of a full adder
	 * clauses are:
	 *   1) -a -b       c_o
	 *   2)  a     c_i -c_o
	 *   3)     b  c_i -c_o
	 *   4)  a  b      -c_o
	 *   5)    -b -c_i  c_o
	 *   6) -a    -c_i  c_o
	 */
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
