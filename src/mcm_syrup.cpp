//
// Created by nfiege on 10/21/22.
//

#include "mcm_syrup.h"

#ifdef USE_SYRUP
#include <iostream>
#include <cstdlib>
#include <mtl/Vec.h>

mcm_syrup::mcm_syrup(const std::vector<std::vector<int>> &C, int timeout, verbosity_mode verbosity, int threads, bool allow_negative_numbers, bool write_cnf)
	: mcm(C, timeout, verbosity, threads, allow_negative_numbers, write_cnf) {}

void *mcm_syrup::timeout_thread(std::pair<int, pthread_t*>* p) {
	sleep(p->first); // first argument: timeout
	pthread_cancel(*p->second); // second argument: thread to cancel
	return nullptr;
}

void *mcm_syrup::worker_thread(std::pair<mcm_syrup*, pthread_t*>* p) {
	auto result = p->first->solver->solve();
	p->first->solved_instance = result == l_True;
	p->first->ran_into_timeout = result == l_Undef;
	pthread_cancel(*p->second); // second argument: thread to cancel
	return nullptr;
}

std::pair<bool, bool> mcm_syrup::check() {
	// use pthread to handle timeout (dirty but ok I guess)
	std::pair <int, pthread_t*> timeout_args = {this->timeout, &this->worker_thread_id};
	std::pair <mcm_syrup*, pthread_t*> worker_args = {this, &this->timeout_thread_id};
	pthread_create(&this->timeout_thread_id, nullptr, reinterpret_cast<void *(*)(void *)>(&mcm_syrup::timeout_thread), &timeout_args);
	pthread_create(&this->worker_thread_id, nullptr, reinterpret_cast<void *(*)(void *)>(&mcm_syrup::worker_thread), &worker_args);
	pthread_join(this->worker_thread_id, nullptr);
	return {this->solved_instance, this->ran_into_timeout};
}

void mcm_syrup::reset_backend(formulation_mode mode) {
	mcm::reset_backend(mode);
	if (mode != formulation_mode::reset_all) return;
	this->glucoseVariableCounter = 0;
	this->solver = std::make_unique<Glucose::MultiSolvers>(this->threads);
	if (this->verbosity != mcm::verbosity_mode::debug_mode) this->solver->setVerbosity(0);
	else this->solver->setVerbosity(2);
}

int mcm_syrup::get_result_value(int var_idx) {
	auto val = this->solver->model[var_idx-1]; // internally, glucose starts with idx 0 for the first variable
	if (val == l_True) return 1;
	return 0;
}

void mcm_syrup::create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) {
	mcm::create_arbitrary_clause(a);
	Glucose::vec<Glucose::Lit> literals;
	for (auto &it : a) {
		auto var = it.first-1;
		if (it.second) literals.push(~Glucose::mkLit(var)); // negate variable
		else literals.push(Glucose::mkLit(var)); // do not negate
	}
	this->solver->addClause(literals);
}

void mcm_syrup::create_new_variable(int idx) {
	mcm::create_new_variable(idx);
	while (idx > this->glucoseVariableCounter) {
		this->glucoseVariableCounter++;
		this->solver->newVar();
	}
}


#endif //USE_SYRUP