//
// Created by nfiege on 9/26/22.
//

#include "mcm_cadical.h"

#ifdef USE_CADICAL

#include <iostream>

mcm_cadical::mcm_cadical(const std::vector<std::vector<int>> &C, int timeout, verbosity_mode verbosity, bool allow_negative_numbers, bool write_cnf)
	: mcm(C, timeout, verbosity, 1, allow_negative_numbers, write_cnf) {}

void mcm_cadical::reset_backend(formulation_mode mode) {
	mcm::reset_backend(mode);
	if (mode != formulation_mode::reset_all) return;
	// create new solver
	this->solver = std::make_unique<CaDiCaL::Solver>();
	// create and attach new terminator
	this->terminator = cadical_terminator(this->timeout);
	this->solver->connect_terminator(&this->terminator);
}

std::pair<bool, bool> mcm_cadical::check() {
	auto stat = this->solver->solve();
	auto sat = stat == CADICAL_SAT;
	auto unsat = stat == CADICAL_UNSAT;
	auto to = !sat and !unsat;
	return {sat, to};
}

int mcm_cadical::get_result_value(int var_idx) {
	return this->solver->val(var_idx) > 0 ? 1 : 0;
}

void mcm_cadical::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	mcm::create_arbitrary_clause(a);
	for (auto &it : a) {
		if (it.second) {
			this->solver->add(-it.first);
		}
		else {
			this->solver->add(it.first);
		}
	}
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