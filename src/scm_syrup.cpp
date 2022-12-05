//
// Created by nfiege on 10/21/22.
//

#include "scm_syrup.h"

#ifdef USE_SYRUP
#include <iostream>
#include <cstdlib>
#include <mtl/Vec.h>

scm_syrup::scm_syrup(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers, bool write_cnf)
	: scm(C, timeout, quiet, threads, allow_negative_numbers, write_cnf) {}

void *scm_syrup::timeout_thread(std::pair<int, pthread_t*>* p) {
	sleep(p->first); // first argument: timeout
	pthread_cancel(*p->second); // second argument: thread to cancel
	return nullptr;
}

void *scm_syrup::worker_thread(std::pair<scm_syrup*, pthread_t*>* p) {
	auto result = p->first->solver->solve();
	p->first->solved_instance = result == l_True;
	p->first->ran_into_timeout = result == l_Undef;
	pthread_cancel(*p->second); // second argument: thread to cancel
	return nullptr;
}

std::pair<bool, bool> scm_syrup::check() {
	// use pthread to handle timeout (dirty but ok I guess)
	std::pair <int, pthread_t*> timeout_args = {this->timeout, &this->worker_thread_id};
	std::pair <scm_syrup*, pthread_t*> worker_args = {this, &this->timeout_thread_id};
	pthread_create(&this->timeout_thread_id, nullptr, reinterpret_cast<void *(*)(void *)>(&scm_syrup::timeout_thread), &timeout_args);
	pthread_create(&this->worker_thread_id, nullptr, reinterpret_cast<void *(*)(void *)>(&scm_syrup::worker_thread), &worker_args);
	pthread_join(this->worker_thread_id, nullptr);
	return {this->solved_instance, this->ran_into_timeout};
}

void scm_syrup::reset_backend(formulation_mode mode) {
	scm::reset_backend(mode);
	if (mode != formulation_mode::reset_all) return;
	this->glucoseVariableCounter = 0;
	this->solver = std::make_unique<Glucose::MultiSolvers>(this->threads);
	if (this->quiet) this->solver->setVerbosity(0);
	else this->solver->setVerbosity(2);
}

int scm_syrup::get_result_value(int var_idx) {
	auto val = this->solver->model[var_idx-1]; // internally, glucose starts with idx 0 for the first variable
	if (val == l_True) return 1;
	return 0;
}

void scm_syrup::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	scm::create_arbitrary_clause(a);
	Glucose::vec<Glucose::Lit> literals;
	for (auto &it : a) {
		auto var = it.first-1;
		if (it.second) literals.push(~Glucose::mkLit(var)); // negate variable
		else literals.push(Glucose::mkLit(var)); // do not negate
	}
	this->solver->addClause(literals);
}

void scm_syrup::create_new_variable(int idx) {
	scm::create_new_variable(idx);
	while (idx > this->glucoseVariableCounter) {
		this->glucoseVariableCounter++;
		this->solver->newVar();
	}
}


#endif //USE_SYRUP