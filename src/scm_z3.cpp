//
// Created by nfiege on 9/30/22.
//

#include "scm_z3.h"

#ifdef USE_Z3

scm_z3::scm_z3(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers, bool write_cnf)
	:	scm(C, timeout, quiet, threads, allow_negative_numbers, write_cnf), solver(this->context) {}

std::pair<bool, bool> scm_z3::check() {
	if (this->timeout > 0) {
		this->solver.set("timeout", (unsigned int)this->timeout*1000);
	}
	this->solver.set("threads", (unsigned int)this->threads);
	auto stat = this->solver.check();
	auto sat = stat == z3::sat;
	auto unsat = stat == z3::unsat;
	auto to = !sat and !unsat;
	return {sat, to};
}

void scm_z3::reset_backend(bool incremental) {
	scm::reset_backend(incremental);
	if (incremental) return;
	this->solver.reset();
	this->variables.clear();
	this->variables.emplace_back(this->context.bool_const("dummy")); // reset variables and add a dummy expression because indices start at 1
}

int scm_z3::get_result_value(int var_idx) {
	return this->solver.get_model().eval(this->variables.at(var_idx)).is_true()?1:0;
}

void scm_z3::create_new_variable(int idx) {
	auto name = std::to_string(idx);
	this->variables.emplace_back(this->context.bool_const(name.c_str()));
}

void scm_z3::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	scm::create_arbitrary_clause(a);
	z3::expr e(this->context);
	if (a.at(0).second) {
		e = not this->variables.at(a.at(0).first);
	}
	else {
		e = this->variables.at(a.at(0).first);
	}
	for (int i=1; i<a.size(); i++) {
		auto &it = a.at(i);
		if (it.second) {
			e = e or (not this->variables.at(it.first));
		}
		else {
			e = e or this->variables.at(it.first);
		}
	}
	this->solver.add(e);
}

#endif //USE_Z3