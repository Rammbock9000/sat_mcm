//
// Created by nfiege on 9/30/22.
//

#ifndef SATSCM_SAT_Z3_H
#define SATSCM_SAT_Z3_H

#ifdef USE_Z3

#include <scm.h>
#include <z3++.h>
#include <chrono>
#include <memory>
#include <utility>
#include <vector>


class scm_z3 : public scm {
public:
	scm_z3(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend() override;
	int get_result_value(int var_idx) override;
	void create_new_variable(int idx) override;

	void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) override;

private:
	z3::context context;
	z3::solver solver;
	std::vector<z3::expr> variables;
};

#endif //USE_Z3

#endif //SATSCM_SAT_Z3_H
