//
// Created by nfiege on 9/30/22.
//

#include "scm_z3.h"

#ifdef USE_Z3

scm_z3::scm_z3(int C, int timeout, bool quiet, int word_size)
	:	scm(C, timeout, quiet, word_size), solver(this->context) {}

std::pair<bool, bool> scm_z3::check() {
	if (this->timeout > 0) {
		this->solver.set("timeout", (unsigned int)this->timeout*1000);
	}
	auto stat = this->solver.check();
	auto sat = stat == z3::sat;
	auto unsat = stat == z3::unsat;
	auto to = !sat and !unsat;
	return {sat, to};
}

void scm_z3::reset_backend() {
	this->solver.reset();
	this->variables.clear();
	this->variables.emplace_back(this->context.bool_const("dummy")); // reset variables and add a dummy expression because indices start at 1
}

int scm_z3::get_result_value(int var_idx) {
	return this->solver.get_model().eval(this->variables.at(var_idx)).is_true()?1:0;
}

void scm_z3::create_1x1_equivalence(int x, int y) {
	// 1)
	this->solver.add((not this->variables.at(x)) or (this->variables.at(y)));
	// 2)
	this->solver.add((this->variables.at(x)) or (not this->variables.at(y)));
}

void scm_z3::create_2x1_mux(int a, int b, int s, int y) {
	// 1)
	this->solver.add((not this->variables.at(a)) or (this->variables.at(s)) or (this->variables.at(y)));
	// 2)
	this->solver.add((not this->variables.at(b)) or (not this->variables.at(s)) or (this->variables.at(y)));
	// 3)
	this->solver.add((this->variables.at(b)) or (not this->variables.at(s)) or (not this->variables.at(y)));
	// 4)
	this->solver.add((this->variables.at(a)) or (this->variables.at(s)) or (not this->variables.at(y)));
	// 5)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(b)) or (this->variables.at(y)));
	// 6)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(y)));
}

void scm_z3::create_2x1_mux_shift_disallowed(int a, int b, int s, int y) {
	// 1)
	this->solver.add((not this->variables.at(a)) or (this->variables.at(y)));
	// 2)
	this->solver.add((not this->variables.at(b)) or (not this->variables.at(s)) or (this->variables.at(y)));
	// 3)
	this->solver.add((this->variables.at(b)) or (not this->variables.at(s)) or (not this->variables.at(y)));
	// 4)
	this->solver.add((this->variables.at(a)) or (this->variables.at(s)) or (not this->variables.at(y)));
	// 5)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(s)));
	// 6)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(y)));
}

void scm_z3::create_2x1_xor(int a, int b, int y) {
	// 1)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(y)));
	// 2)
	this->solver.add((this->variables.at(a)) or (not this->variables.at(b)) or (this->variables.at(y)));
	// 3)
	this->solver.add((not this->variables.at(a)) or (this->variables.at(b)) or (this->variables.at(y)));
	// 4)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(b)) or (not this->variables.at(y)));
}

void scm_z3::create_add_sum(int a, int b, int c_i, int s) {
	// 1)
	this->solver.add((this->variables.at(a)) or (not this->variables.at(b)) or (this->variables.at(c_i)) or (this->variables.at(s)));
	// 2)
	this->solver.add((not this->variables.at(a)) or (this->variables.at(b)) or (this->variables.at(c_i)) or (this->variables.at(s)));
	// 3)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (this->variables.at(c_i)) or (not this->variables.at(s)));
	// 4)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(b)) or (this->variables.at(c_i)) or (not this->variables.at(s)));
	// 5)
	this->solver.add((this->variables.at(a)) or (not this->variables.at(b)) or (not this->variables.at(c_i)) or (not this->variables.at(s)));
	// 6)
	this->solver.add((not this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(c_i)) or (not this->variables.at(s)));
	// 7)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(c_i)) or (this->variables.at(s)));
	// 8)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(b)) or (not this->variables.at(c_i)) or (this->variables.at(s)));
}

void scm_z3::create_add_carry(int a, int b, int c_i, int c_o) {
	// 1)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(b)) or (this->variables.at(c_o)));
	// 2)
	this->solver.add((this->variables.at(a)) or (this->variables.at(c_i)) or (not this->variables.at(c_o)));
	// 3)
	this->solver.add((this->variables.at(b)) or (this->variables.at(c_i)) or (not this->variables.at(c_o)));
	// 4)
	this->solver.add((this->variables.at(a)) or (this->variables.at(b)) or (not this->variables.at(c_o)));
	// 5)
	this->solver.add((not this->variables.at(b)) or (not this->variables.at(c_i)) or (this->variables.at(c_o)));
	// 6)
	this->solver.add((not this->variables.at(a)) or (not this->variables.at(c_i)) or (this->variables.at(c_o)));
}

void scm_z3::force_bit(int x, int val) {
	if (val == 1)
		this->solver.add(this->variables.at(x));
	else
		this->solver.add(not this->variables.at(x));
}

void scm_z3::forbid_number(const std::vector<int> &x, int val) {
	auto num_bits = (int)x.size();
	z3::expr e(this->context);
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			if (i == 0)
				e = not this->variables.at(x[i]);
			else
				e = e or (not this->variables.at(x[i]));
		}
		else {
			if (i == 0)
				e = this->variables.at(x[i]);
			else
				e = e or this->variables.at(x[i]);
		}
	}
	this->solver.add(e);
}

void scm_z3::force_number(const std::vector<int> &x, int val) {
	auto num_bits = (int)x.size();
	for (int i=0; i<num_bits; i++) {
		auto bit = (val >> i) & 1;
		if (bit == 1) {
			this->solver.add(this->variables.at(x[i]));
		}
		else {
			this->solver.add(not this->variables.at(x[i]));
		}
	}
}

void scm_z3::create_new_variable(int idx) {
	auto name = std::to_string(idx);
	this->variables.emplace_back(this->context.bool_const(name.c_str()));
}

#endif //USE_Z3