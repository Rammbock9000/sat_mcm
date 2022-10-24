//
// Created by nfiege on 10/21/22.
//

#ifndef SATSCM_SCM_SYRUP_H
#define SATSCM_SCM_SYRUP_H

#ifdef USE_SYRUP

#include <scm.h>
#include <parallel/MultiSolvers.h>
#include <core/SolverTypes.h>
#include <chrono>
#include <memory>
#include <utility>
#include <vector>
#include <pthread.h>


class scm_syrup : public scm {
public:
	scm_syrup(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers, bool write_cnf);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend() override;
	int get_result_value(int var_idx) override;
	void create_new_variable(int idx) override;

	void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) override;

private:
	std::unique_ptr<Glucose::MultiSolvers> solver;
	int glucoseVariableCounter = -1;

	static void* timeout_thread(std::pair<int, pthread_t*>* p);
	static void* worker_thread(std::pair<scm_syrup*, pthread_t*>* p);
	pthread_t worker_thread_id = 0;
	pthread_t timeout_thread_id = 0;

	bool ran_into_timeout = false;
	bool solved_instance = false;
};

#endif //USE_SYRUP
#endif //SATSCM_SCM_SYRUP_H
