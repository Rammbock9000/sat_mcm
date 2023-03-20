//
// Created by nfiege on 10/21/22.
//

#ifndef SATMCM_MCM_SYRUP_H
#define SATMCM_MCM_SYRUP_H

#ifdef USE_SYRUP

#include <mcm.h>
#include <parallel/MultiSolvers.h>
#include <core/SolverTypes.h>
#include <chrono>
#include <memory>
#include <utility>
#include <vector>
#include <pthread.h>


class mcm_syrup : public mcm {
public:
	mcm_syrup(const std::vector<int> &C, int timeout, verbosity_mode verbosity, int threads, bool allow_negative_numbers, bool write_cnf);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend(formulation_mode mode) override;
	int get_result_value(int var_idx) override;
	void create_new_variable(int idx) override;

	void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) override;

private:
	std::unique_ptr<Glucose::MultiSolvers> solver;
	int glucoseVariableCounter = -1;

	static void* timeout_thread(std::pair<int, pthread_t*>* p);
	static void* worker_thread(std::pair<mcm_syrup*, pthread_t*>* p);
	pthread_t worker_thread_id = 0;
	pthread_t timeout_thread_id = 0;

	bool ran_into_timeout = false;
	bool solved_instance = false;
};

#endif //USE_SYRUP
#endif //SATMCM_MCM_SYRUP_H
