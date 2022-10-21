//
// Created by nfiege on 9/26/22.
//

#include "scm_cadical.h"

#ifdef USE_CADICAL

#include <iostream>

scm_cadical::scm_cadical(const std::vector<int> &C, int timeout, bool quiet, bool allow_negative_numbers)
	: scm(C, timeout, quiet, 1, allow_negative_numbers) {}

void scm_cadical::reset_backend() {
	scm::reset_backend();
	// create new solver
	this->solver = std::make_unique<CaDiCaL::Solver>();
	// create and attach new terminator
	this->terminator = cadical_terminator(this->timeout);
	this->solver->connect_terminator(&this->terminator);
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

void scm_cadical::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	scm::create_arbitrary_clause(a);
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