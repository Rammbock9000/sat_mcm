//
// Created by nfiege on 9/26/22.
//

#ifndef SATSCM_SCM_H
#define SATSCM_SCM_H

#include <map>
#include <set>
#include <vector>

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
	scm(int C, int timeout, bool quiet);
	/*!
	 * solve the problem
	 */
	void solve();
	/*!
	 * print solution values
	 */
	void print_solution();

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

	///////////////////////////////////////////////////
	//// create clauses for the following circuits ////
	///////////////////////////////////////////////////
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
	virtual void forbid_number(const std::vector<int> &x, int num);
	/*!
	 * force x == num
	 * @param x vector that contains all bits
	 * @param num
	 */
	virtual void force_number(const std::vector<int> &x, int num);

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
	int variable_counter;
	/*!
	 * count #constraints
	 */
	int constraint_counter;

	/*!
	 * the constant by which we want to multiply
	 */
	int C;
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
	 * number of shifted bits of the output to re-create the original constant
	 */
	int output_shift;
	/*!
	 * if we found a solution, yet
	 */
	bool found_solution = false;
	/*!
	 * if we ran into a timeout during solving
	 */
	bool ran_into_timeout = false;
	/*!
	 * solver timeout
	 */
	int timeout;
	/*!
	 * suppress debug outputs if quiet = true
	 */
	bool quiet;

private:
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
	 * node idx -> 1/0
	 */
	std::map<int, int> shift_input_select;
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
	void create_constant_zero_variable();
	void create_input_select_mux_variables(int idx);
	void create_input_select_selection_variables(int idx);
	void create_input_shift_select_variable(int idx);
	void create_shift_select_output_variables(int idx);
	void create_input_shift_value_variables(int idx);
	void create_shift_internal_variables(int idx);
	void create_input_negate_select_variable(int idx);
	void create_negate_select_output_variables(int idx);
	void create_input_negate_value_variable(int idx);
	void create_xor_output_variables(int idx);
	void create_adder_internal_variables(int idx);
	void create_output_value_variables(int idx);

	////////////////////////////////
	//// CREATE ALL CONSTRAINTS ////
	////////////////////////////////
	void create_input_output_constraints();
	void create_input_select_constraints(int idx);
	void create_input_select_limitation_constraints(int idx);
	void create_shift_limitation_constraints(int idx);
	void create_shift_select_constraints(int idx);
	void create_shift_constraints(int idx);
	void create_negate_select_constraints(int idx);
	void create_xor_constraints(int idx);
	void create_adder_constraints(int idx);

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
	std::map<std::tuple<int, input_direction, int>, int> input_select_selection_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_value_variables;
	/*!
	 * node idx -> variable idx
	 */
	std::map<int, int> input_shift_select_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 */
	std::map<std::tuple<int, input_direction, int>, int> shift_select_output_variables;
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
	std::map<std::pair<int, int>, int> adder_internal_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! node idx = 0 is the input node with constant value 0
	 */
	std::map<std::pair<int, int>, int> output_value_variables;
	/*!
	 * a variable that is equal to a constant zero (needed for shifter)
	 * variable idx
	 */
	int constant_zero_variable;
};


#endif //SATSCM_SCM_H
