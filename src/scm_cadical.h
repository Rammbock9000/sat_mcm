//
// Created by nfiege on 9/26/22.
//

#ifndef SATSCM_SCM_CADICAL_H
#define SATSCM_SCM_CADICAL_H

#include <scm.h>
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

class scm_cadical : public scm {

#define CADICAL_SAT 10
#define CADICAL_UNSAT 20

public:
	scm_cadical(int C, int timeout, bool quiet);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend() override;
	int get_result_value(int var_idx) override;

	void create_1x1_equivalence(int x, int y) override;
	void create_2x1_mux(int a, int b, int s, int y) override;
	void create_2x1_xor(int a, int b, int y) override;
	void create_add_sum(int a, int b, int c_i, int s) override;
	void create_add_carry(int a, int b, int c_i, int c_o) override;
	void force_bit(int x, int val) override;
	void forbid_number(const std::vector<int> &x, int val) override;
	void force_number(const std::vector<int> &x, int val) override;

private:
	std::unique_ptr<CaDiCaL::Solver> solver;
	cadical_terminator terminator;
};

#endif //SATSCM_SCM_CADICAL_H
