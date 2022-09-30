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
	scm_z3(int C, int timeout, bool quiet, int word_size, int threads);

protected:
	std::pair<bool, bool> check() override;
	void reset_backend() override;
	int get_result_value(int var_idx) override;
	void create_new_variable(int idx) override;

	void create_1x1_equivalence(int x, int y) override;
	void create_2x1_mux(int a, int b, int s, int y) override;
	void create_2x1_mux_shift_disallowed(int a, int b, int s, int y) override;
	void create_2x1_xor(int a, int b, int y) override;
	void create_add_sum(int a, int b, int c_i, int s) override;
	void create_add_carry(int a, int b, int c_i, int c_o) override;
	void force_bit(int x, int val) override;
	void forbid_number(const std::vector<int> &x, int val) override;
	void force_number(const std::vector<int> &x, int val) override;

private:
	z3::context context;
	z3::solver solver;
	std::vector<z3::expr> variables;
};

#endif //USE_Z3

#endif //SATSCM_SAT_Z3_H
