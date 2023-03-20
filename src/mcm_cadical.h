//
// Created by nfiege on 9/26/22.
//

#ifndef SATMCM_SCM_CADICAL_H
#define SATMCM_SCM_CADICAL_H

#ifdef USE_CADICAL

#include <mcm.h>
#include <cadical.hpp>
#include <chrono>
#include <memory>
#include <utility>

class cadical_terminator : public CaDiCaL::Terminator {
public:
	explicit cadical_terminator(double timeout = 0.0);
	bool terminate () override;
	void reset(double newTimeout);
	double get_elapsed_time() const;
private:
	double max_time;
	std::chrono::steady_clock::time_point timer_start;
};

class mcm_cadical : public mcm {

#define CADICAL_SAT 10
#define CADICAL_UNSAT 20

public:
	mcm_cadical(const std::vector<int> &C, int timeout, verbosity_mode verbosity, bool allow_negative_numbers, bool write_cnf);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend(formulation_mode mode) override;
	int get_result_value(int var_idx) override;

	void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) override;

private:
	std::unique_ptr<CaDiCaL::Solver> solver;
	cadical_terminator terminator;
};

#endif //USE_CADICAL

#endif //SATMCM_SCM_CADICAL_H
