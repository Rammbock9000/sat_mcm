//
// Created by nfiege on 9/26/22.
//

#ifndef SATSCM_SCM_H
#define SATSCM_SCM_H

#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <cstdint>

#define SHIFT_SELECT_OLD 0

class scm {
public:
	enum input_direction {
		left, right
	};
	const std::set<input_direction> input_directions = {left, right};
	/*!
	 * constructor
	 * @param C the constant we want to compute
	 * @param timeout in seconds
	 * @param quiet true/false
	 */
	scm(const std::vector<int> &C, int timeout, bool quiet, int threads, bool allow_negative_numbers, bool write_cnf);
	/*!
	 * define the minimum number of needed adders to help the algorithm converge faster
	 * @param new_min_add value
	 */
	void set_min_add(int new_min_add);
	/*!
	 * also minimize full adders for the optimal number of adder nodes during this->solve()
	 */
	void also_minimize_full_adders();
	/*!
	 * also allow a shift at each node's output during this->solve()
	 * this reduces the number of needed adders in very few corner cases by 1 at the cost of an increased runtime
	 */
	void allow_node_output_shift();
	/*!
	 * allow the solver to choose whether to implement a coefficient C_i or -C_i
	 * depending on the implementation costs
	 * @param only_apply_to_negative_coefficients only allow this option for negative requested coefficients
	 */
	void ignore_sign(bool only_apply_to_negative_coefficients);
	/*!
	 * solve the problem
	 */
	void solve();
	/*!
	 * print solution values
	 */
	void print_solution();
	/*!
	 * sign extend x and return that badboy
	 * @param x 2's complement number with w bits
	 * @param w word size of x
	 * @return sign extended x
	 */
	static int64_t sign_extend(int64_t x, int w);
	/*!
	 * visit 'https://gitlab.com/kumm/pagsuite' for info
	 * @return a string that uniquely describes the resulting adder graph
	 */
	std::string get_adder_graph_description();
	/*!
	 * @return whether the computed solution is optimal
	 *   -> pair.first: optimal w.r.t. number of add operations
	 *   -> pair.second: optimal w.r.t. full adder count (low level costs)
	 */
	std::pair<int, int> solution_is_optimal();

protected:
	/*!
	 * check feasibility after constructing the problem
	 * @return whether the problem is feasible
	 */
	virtual std::pair<bool, bool> check();
	/*!
	 * get result value for the variable with index "var_idx" if a solution was found
	 * @param var_idx
	 * @return the value
	 */
	virtual int get_result_value(int var_idx);
	/*!
	 * reset the solver backend
	 */
	virtual void reset_backend();

	/*!
	 * create new variable (if backend needs it)
	 * @param idx variable index (=name)
	 */
	virtual void create_new_variable(int idx);
	/*!
	 * helper function to create an arbitrary clause:
	 * @param a < variable idx, negate >
	 *   -> negate the variable if negate == true
	 */
	virtual void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a);

	///////////////////////////////////////////////////
	//// create clauses for the following circuits ////
	///////////////////////////////////////////////////
	/*!
	 * disallow shifting bits that are not equal to the sign bit
	 * clauses are:
	 *   1) -sel -s_a  a
	 *   2) -sel  s_a -a
	 * @param sel
	 * @param s_a
	 * @param a
	 */
	virtual void create_signed_shift_overflow_protection(int sel, int s_a, int a);
	/*!
	 * disallow overflows for signed additions/subtractions (sub=1: subtraction, sub=0: addition)
	 * clauses are:
	 *   1)  sub  s_a  s_b -s_y
	 *   2)  sub -s_a -s_b  s_y
	 *   3) -sub  s_a -s_b -s_y
	 *   4) -sub -s_a  s_b  s_y
	 * @param sub
	 * @param s_a
	 * @param s_b
	 * @param s_y
	 */
	virtual void create_signed_add_overflow_protection(int sub, int s_a, int s_b, int s_y);
	/*!
	 * force x_0 or x_1 or ... x_n = 1
	 * clauses are:
	 *   1) x_0 x_1 x_2 ...
	 * @param a
	 * @param b
	 */
	virtual void create_or(std::vector<int> &x);
	/*!
	 * force a -> b
	 * clauses are:
	 *   1) -a  b
	 * @param a
	 * @param b
	 */
	virtual void create_1x1_implication(int a, int b);
	/*!
	 * force a -> !b
	 * clauses are:
	 *   1) -a  -b
	 * @param a
	 * @param b
	 */
	virtual void create_1x1_negated_implication(int a, int b);
	/*!
	 * force a -> (b_0 or b_1 or ...)
	 *   1) -a  b_0  b_1 ...
	 * @param a
	 * @param b
	 */
	virtual void create_1xN_implication(int a, const std::vector<int> &b);
	/*!
	 * force (a_0 and a_1 and ...) -> (b_0 or b_1 or ...)
	 *   1) -a_0 -a_1 ...  b_0  b_1 ...
	 * @param a
	 * @param b
	 */
	virtual void create_MxN_implication(const std::vector<int> &a, const std::vector<int> &b);
	/*!
	 * force y = x
	 * clauses are:
	 *   1) -x  y
	 *   2)  x -y
	 * @param x
	 * @param y
	 */
	virtual void create_1x1_equivalence(int x, int y);
	/*!
	 * force y = s ? a : b
	 * clauses are:
	 *   1) -a     s  y
	 *   2)    -b -s  y
	 *   3)     b -s -y
	 *   4)  a     s -y
	 *   5) -a -b     y
	 *   6)  a  b    -y
	 * @param a
	 * @param b
	 * @param s
	 * @param y
	 */
	virtual void create_2x1_mux(int a, int b, int s, int y);
	/*!
	 * force y = s ? a : b AND not (s and a)
	 * clauses are:
	 *   1) -a        y
	 *   2)    -b -s  y
	 *   3)     b -s -y
	 *   4)  a     s -y
	 *   5) -a    -s
	 *   6)  a  b    -y
	 * @param a
	 * @param b
	 * @param s
	 * @param y
	 */
	virtual void create_2x1_mux_shift_disallowed(int a, int b, int s, int y);
	/*!
	 * force y = s ? a : b where b = 0
	 * => y = !s and a
	 * clauses are:
	 *   1)    -s -y
	 *   2)  a    -y
	 *   3) -a  s  y
	 * @param a
	 * @param s
	 * @param y
	 */
	virtual void create_2x1_mux_zero_const(int a, int s, int y);
	/*!
	 * force y = a XOR b
	 * clauses are:
	 *   1)  a  b -y
	 *   2)  a -b  y
	 *   3) -a  b  y
	 *   4) -a -b -y
	 * @param a
	 * @param b
	 * @param y
	 */
	virtual void create_2x1_xor(int a, int b, int y);
	/*!
	 * force s to be the sum output of a full adder
	 * clauses are:
	 *   1)  a -b  c_i  s
	 *   2) -a  b  c_i  s
	 *   3)  a  b  c_i -s
	 *   4) -a -b  c_i -s
	 *   5)  a -b -c_i -s
	 *   6) -a  b -c_i -s
	 *   7)  a  b -c_i  s
	 *   8) -a -b -c_i  s
	 * @param a
	 * @param b
	 * @param c_i
	 * @param s
	 */
	virtual void create_add_sum(int a, int b, int c_i, int s);
	/*!
	 * force c_o to be the carry output of a full adder
	 * clauses are:
	 *   1) -a -b       c_o
	 *   2)  a     c_i -c_o
	 *   3)     b  c_i -c_o
	 *   4)  a  b      -c_o
	 *   5)    -b -c_i  c_o
	 *   6) -a    -c_i  c_o
	 * @param a
	 * @param b
	 * @param c_i
	 * @param c_o
	 */
	virtual void create_add_carry(int a, int b, int c_i, int c_o);
	/*!
	 * these clauses are not needed but increase performance during unit propagation
	 * clauses are:
	 *   1)  a         -s -c_o
	 *   2)     b      -s -c_o
	 *   3)        c_i -s -c_o
	 *   4) -a          s  c_o
	 *   5)    -b       s  c_o
	 *   6)       -c_i  s  c_o
	 * @param a
	 * @param c_i
	 * @param s
	 * @param c_o
	 */
	virtual void create_add_redundant(int a, int b, int c_i, int s, int c_o);
	/*!
	 * set x = val
	 * @param x
	 * @param val must be 0 or 1
	 */
	virtual void force_bit(int x, int val);
	/*!
	 * force x != num
	 * @param x vector that contains all bits
	 * @param num
	 */
	virtual void forbid_number(const std::vector<int> &x, int val);
	/*!
	 * force x == num
	 * @param x vector that contains all bits
	 * @param num
	 */
	virtual void force_number(const std::vector<int> &x, int val);

	/*!
	 * @param n
	 * @return ceil(log2(n))
	 */
	int ceil_log2(int n);
	/*!
	 * @param n
	 * @return floor(log2(n))
	 */
	int floor_log2(int n);

	/*!
	 * count #variables
	 */
	int variable_counter = 0;
	/*!
	 * count #constraints
	 */
	int constraint_counter = 0;

	/*!
	 * the constant by which we want to multiply
	 */
	std::vector<int> C;
	/*!
	 * store info whether or not the negative version of a coefficient was requested by the user
	 */
	std::map<int, bool> negative_coeff_requested;
	/*!
	 * store info whether the solver can decide to implement C or -C
	 * this is only relevant when ...
	 *   ... minimizing full adders
	 *   ... the solver is allowed to use negative numbers
	 *   ... doing MCM
	 */
	std::map<int, bool> sign_inversion_allowed;
	/*!
	 * word size of all operations
	 */
	int word_size;
	/*!
	 * maximum allowed shift
	 */
	int max_shift;
	/*!
	 * word size of the corresponding node input
	 */
	int shift_word_size;
	/*!
	 * current number of adders
	 */
	int num_adders = 0;
	/*!
	 * keep track how to get the requested constants from the computed nodes (relevant if constants are negative or even)
	 * requested constant -> < adder node output, number of shifted bits >
	 * e.g. "18 -> < 9, 1 >" because 18 is computed from 9 left-shifted by 1 bit
	 */
	std::map<int, std::pair<int, int>> requested_constants;
	/*!
	 * if we found a solution, yet
	 */
	bool found_solution = false;
	/*!
	 * if we ran into a timeout during solving
	 */
	bool ran_into_timeout = false;
	/*!
	 * the solution has the optimal number of adders
	 */
	bool num_add_opt = false;
	/*!
	 * the solution has the optimal number of full adders
	 */
	bool num_FA_opt = false;
	/*!
	 * solver timeout
	 */
	int timeout;
	/*!
	 * suppress debug outputs if quiet = true
	 */
	bool quiet;
	/*!
	 * the number of CPU threads the backend is allowed to use
	 */
	int threads;
	/*!
	 * also write the corresponding cnf files for all solving attempts
	 * e.g. 521_2.cnf for C = 521 and #adders = 2
	 */
	bool write_cnf;
	/*!
	 * whether we are performing all computation in 2's complement
	 * i.e. at least one coefficient is negative
	 */
	bool calc_twos_complement;
	/*!
	 * whether we also minimize the number of full adders for the minimum number of adders
	 */
	bool minimize_full_adders = false;
	/*!
	 * whether we also allow a shift at each node's output
	 */
	bool enable_node_output_shift = false;

private:
	/*!
	 * limit on the number of full adders used
	 * -1 = no limit
	 */
	int max_full_adders = -1;
	/*!
	 * store all cnf clauses for cnf file generation
	 */
	std::stringstream cnf_clauses;
	/*!
	 * creates a .cnf file for the current SAT problem
	 */
	void create_cnf_file();
	/*!
	 * get solution from backend and store result in containers below
	 */
	void get_solution_from_backend();
	/*!
	 * verify whether the found solution is valid
	 * @return if it is valid
	 */
	bool solution_is_valid();
	/*!
	 * < node idx, left/right > -> int value
	 */
	std::map<std::pair<int, input_direction>, int> input_select;
	/*!
	 * < node idx, left/right > -> int value
	 */
	std::map<std::pair<int, input_direction>, int> input_select_mux_output;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> shift_value;
	/*!
	 * node idx -> 1/0
	 */
	std::map<int, int> negate_select;
	/*!
	 * node idx -> 1/0
	 */
	std::map<int, int> subtract;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> add_result_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> post_adder_shift_value;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> output_values;
	/*!
	 * create backend solver variables and keep track of their indices
	 */
	void create_variables();
	/*!
	 * create solver constraints
	 */
	void create_constraints();
	/*!
	 * construct the problem (everything needed for solving)
	 */
	void construct_problem();
	/*!
	 * optimize #adders or #full_adders within this loop
	 */
	void optimization_loop();
	/*!
	 * cache values for ceil(log2(n))
	 */
	std::map<int, int> ceil_log2_cache;
	/*!
	 * cache values for floor(log2(n))
	 */
	std::map<int, int> floor_log2_cache;

	//////////////////////////////
	//// CREATE ALL VARIABLES ////
	//////////////////////////////
	void create_input_node_variables(); // idx = 0 is the input node that has a constant value 0 as output
	void create_input_select_mux_variables(int idx);
	void create_input_select_selection_variables(int idx);
	void create_input_shift_value_variables(int idx);
	void create_shift_internal_variables(int idx);
	void create_input_negate_select_variable(int idx);
	void create_negate_select_output_variables(int idx);
	void create_input_negate_value_variable(int idx);
	void create_xor_output_variables(int idx);
	void create_adder_internal_variables(int idx);
	void create_post_adder_input_shift_value_variables(int idx);
	void create_post_adder_shift_variables(int idx);
	void create_output_value_variables(int idx);
	void create_mcm_output_variables(int idx);
	void create_full_adder_alloc_variables(int idx);

	////////////////////////////////
	//// CREATE ALL CONSTRAINTS ////
	////////////////////////////////
	void create_input_output_constraints();
	void create_input_select_constraints(int idx);
	void create_input_select_limitation_constraints(int idx);
	void create_shift_limitation_constraints(int idx);
	void create_shift_constraints(int idx);
	void create_negate_select_constraints(int idx);
	void create_xor_constraints(int idx);
	void create_adder_constraints(int idx);
	void create_post_adder_shift_limitation_constraints(int idx);
	void create_post_adder_shift_constraints(int idx);
	void create_full_adder_allocation_constraints(int idx);
	void create_full_adder_overlap_constraints(int idx_1);
	void create_odd_fundamentals_constraints(int idx);
	void create_mcm_output_constraints();

	///////////////////////////////////
	//// INDICES FOR ALL VARIABLES ////
	//////////////////////////////////////////////////
	//// THE FIRST VARIABLE INDEX STARTS WITH 1!! ////
	//////////////////////////////////////////////////
	/*!
	 * < node idx, left/right, mux idx, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int, int>, int> input_select_mux_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_select_mux_output_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_select_selection_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_value_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> input_shift_value_variables;
	/*!
	 * < node idx, mux stage, bit > -> variable idx
	 */
	std::map<std::tuple<int, int, int>, int> shift_internal_mux_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! identical to the last MUX stage of shift_internal_mux_output_variables
	 */
	std::map<std::pair<int, int>, int> shift_output_variables;
	/*!
	 * node idx -> variable idx
	 */
	std::map<int, int> input_negate_select_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> negate_select_output_variables;
	/*!
	 * node idx -> variable idx
	 */
	std::map<int, int> input_negate_value_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> xor_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> adder_carry_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * used for optimized adder clauses
	 */
	std::map<std::pair<int, int>, int> adder_XOR_internal_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! node idx = 0 is the input node with constant value 0
	 */
	std::map<std::pair<int, int>, int> adder_output_value_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> input_post_adder_shift_value_variables;
	/*!
	 * < node idx, mux stage, bit > -> variable idx
	 */
	std::map<std::tuple<int, int, int>, int> post_adder_shift_internal_mux_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! identical to the last MUX stage of shift_internal_mux_output_variables
	 */
	std::map<std::pair<int, int>, int> post_adder_shift_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! node idx = 0 is the input node with constant value 0
	 */
	std::map<std::pair<int, int>, int> output_value_variables;
	/*!
	 * < node idx, mcm constant > -> variable idx
	 */
	std::map<std::pair<int, int>, int> mcm_output_variables;
	/*!
	 * < node idx, bit, full adder idx > -> variable idx
	 */
	std::map<std::tuple<int, int, int>, int> full_adder_alloc_variables;
};


#endif //SATSCM_SCM_H
