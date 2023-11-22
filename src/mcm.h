//
// Created by nfiege on 9/26/22.
//

#ifndef SATMCM_MCM_H
#define SATMCM_MCM_H

#include <map>
#include <set>
#include <vector>
#include <sstream>
#include <cstdint>
#include <limits>

#define SHIFT_SELECT_OLD 0
#define FULL_ADDERS_UNLIMITED std::numeric_limits<long>::min()

class mcm {
public:
	enum input_direction {
		left, right
	};
	const std::set<input_direction> input_directions = {left, right};
	enum formulation_mode {
		reset_all, all_FA_clauses, only_FA_limit
	};
	enum verbosity_mode {
		debug_mode, normal_mode, quiet_mode
	};
	/*!
	 * constructor
	 * @param C the vector of constants we want to compute
	 * @param timeout in seconds
	 * @param quiet true/false
	 */
	mcm(const std::vector<std::vector<int>> &C, int timeout, verbosity_mode verbosity, int threads, bool allow_negative_numbers, bool write_cnf);
	/*!
	 * define the minimum number of needed adders to help the algorithm converge faster
	 * @param new_min_add value
	 */
	void set_min_add(int new_min_add);
	/*!
	 * define whether all solutions shall be enumerated or only the optimal solution is of interest
	 * @param new_enumerate_all value
	 */
	void set_enumerate_all(bool new_enumerate_all);
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
	 * only permit solutions with minimal adder depth by calling this function prior to solving
	 */
	void minimize_adder_depth();
	/*!
	 * assume that the adder graph will be implemented with a pipeline stage after each adder stage
	 * => optimize for that
	 * => minimize the total number of registered operations
	 * => i.e., an adder is "as expensive" as a register
	 */
	void enable_pipelining();
	/*!
	 * force all outputs into the same pipeline stage
	 */
	void equalize_output_stages();
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
	virtual void reset_backend(formulation_mode mode);

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
	 * clauses for a full adder (i.e., 3:2 compressor)
	 * @param a < variable, whether the bit should be negated at the adder input >
	 * @param b < variable, whether the bit should be negated at the adder input >
	 * @param c_i < variable, whether the bit should be negated at the adder input >
	 * @param sum < variable, whether the bit should be negated at the adder output >
	 * @param c_o < variable, whether the bit should be negated at the adder output >
	 */
	virtual void create_full_adder(std::pair<int, bool> a, std::pair<int, bool> b, std::pair<int, bool> c_i, std::pair<int, bool> sum, std::pair<int, bool> c_o={-1, false});
	/*!
	 * clauses for a half adder (i.e., 2:2 compressor, a full adder with c_i=0)
	 * @param a < variable, whether the bit should be negated at the adder input >
	 * @param b < variable, whether the bit should be negated at the adder input >
	 * @param sum < variable, whether the bit should be negated at the adder output >
	 * @param c_o < variable, whether the bit should be negated at the adder output >
	 */
	virtual void create_half_adder(std::pair<int, bool> a, std::pair<int, bool> b, std::pair<int, bool> sum, std::pair<int, bool> c_o={-1, false});
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
	 * force !a -> b
	 * clauses are:
	 *   1)  a   b
	 * @param a
	 * @param b
	 */
	virtual void create_1x1_reversed_negated_implication(int a, int b);
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
	 * force y = not (a XOR b)
	 * clauses are:
	 *   1)  a  b  y
	 *   2)  a -b -y
	 *   3) -a  b -y
	 *   4) -a -b  y
	 * @param a
	 * @param b
	 * @param y
	 */
	virtual void create_2x1_equiv(int a, int b, int y);
	/*!
	 * force y = a OR b
	 * clauses are:
	 *   1)  a  b -y
	 *   2) -a     y
	 *   3)    -b  y
	 * @param a
	 * @param b
	 * @param y
	 */
	virtual void create_2x1_or(int a, int b, int y);
	/*!
	 * force y = a AND b
	 * clauses are:
	 *   1) -a -b  y
	 *   2)  a    -y
	 *   3)     b -y
	 * @param a
	 * @param b
	 * @param y
	 */
	virtual void create_2x1_and(int a, int b, int y);
	/*!
	 * force y = a AND (not b)
	 * clauses are:
	 *   1) -a  b  y
	 *   2)  a    -y
	 *   3)    -b -y
	 * @param a
	 * @param b
	 * @param y
	 */
	virtual void create_2x1_and_b_inv(int a, int b, int y);
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
	 * create clauses for one bit of a "ripple-carry" version of the computation "c = max(a,b)", where a,b,c are integers
	 * @param a_i current bit for a
	 * @param b_i current bit for b
	 * @param x_i current bit that determines whether "a >= b" was already determined in a previous stage
	 * @param y_i current bit that determines whether x_i is valid (i.e., could be determined, yet)
	 * @param c_o current output bit for the result
	 * @param x_o serves as "x_i" in the next stage
	 * @param y_o serves as "y_i" in the next stage
	 */
	virtual void create_max_cell(int a_i, int b_i, int x_i, int y_i, int c_o, int x_o, int y_o);

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
	 * the constant(s) by which we want to multiply
	 */
	std::vector<std::vector<int>> C;
	/*!
	 * store info whether or not the negative version of a coefficient was requested by the user
	 */
	std::map<std::vector<int>, bool> inverted_coeff_requested;
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
	 * word size for adder depth computation (relevant for min adder depth constraint and for pipelining)
	 */
	int adder_depth_word_size;
	/*!
	 * whether to permit only solutions with minimal adder depth
	 */
	bool force_min_adder_depth = false;
	/*!
	 * value for the optimum adder depth
	 */
	int opt_adder_depth = 0;
	/*!
	 * used to compute this->opt_adder_depth before solving
	 */
	void compute_opt_adder_depth_value();
	/*!
	 * preprocess the requested constants
	 */
	void preprocess_constants();
	/*!
	 * whether to assume a fully pipelined adder graph
	 * -> i.e., minimize the number of registered operations (#reg_add + #reg)
	 */
	bool pipelining_enabled = false;
	/*!
	 * force all outputs into the same pipeline stage
	 */
	bool force_output_stages_equal = false;
	/*!
	 * current number of adders
	 */
	int num_adders = 0;
	/*!
	 * keep track how to get the requested constants from the computed nodes (relevant if constants are negative or even)
	 * requested constant -> < adder node output, number of shifted bits >
	 * e.g. "18 -> < 9, 1 >" because 18 is computed from 9 left-shifted by 1 bit
	 */
	std::map<std::vector<int>, std::pair<std::vector<int>, int>> requested_vectors;
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
	 * timeout for FA minimization (duh!)
	 */
	int fa_minimization_timeout;
	/*!
	 * define cout level
	 */
	verbosity_mode verbosity;
	/*!
	 * the number of CPU threads the backend is allowed to use
	 */
	int threads;
	/*!
	 * also write the corresponding cnf files for all solving attempts
	 * e.g. 521_2.cnf for C = 521 and #adders = 2
	 * e.g. 412_532_2.cnf for V = <412,532> and #adders  = 2
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
    /*!
    * whether the solver supports incremental solving
    * @return default = true (overload derived class if this is not the case)
    */
    virtual bool supports_incremental_solving() { return true; }

private:
	/*!
	 * whether we want to enumerate all solutions for minimum adder count
	 */
	bool enumerate_all = false;
	/*!
	 * solving method for this->enumerate_all == true
	 */
	void solve_enumeration();
	/*!
	 * solving method for this->enumerate_all == false (standard case)
	 */
	void solve_standard();
	/*!
	 * limit on the number of full adders used
	 * FULL_ADDERS_UNLIMITED = no limit
	 */
	long int max_full_adders = FULL_ADDERS_UNLIMITED;
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
	std::map<std::tuple<int, input_direction>, int> input_select;
	/*!
	 * < node idx, left/right > -> int value
	 * //CMM Dimension
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_select_mux_output;
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
	 * CMM Dimension
	 */
    std::map<std::pair<int, int>, int> add_result_values;
	/*!
	 * node idx -> int value
	 */
    std::map<int, int> post_adder_shift_value;
	/*!
	 * node idx -> int value
	 * CMM Dimension
	 */
	std::map<std::pair<int, int>, int> output_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> coeff_word_size_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> can_cut_msb_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> coeff_word_size_sum_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> shift_gain_values;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> shift_sum_values;
	/*!
	 * node idx -> 1/0
	 */
	std::map<int, int> is_register;
	/*!
	 * node idx -> int value
	 */
	std::map<int, int> pipeline_stage;
	/*!
	 *  int value
	 */
	int num_FAs_value;
	/*!
	 * create backend solver variables and keep track of their indices
	 */
	void create_variables();
	/*!
	 * create solver constraints
	 */
	void create_constraints(formulation_mode mode);
	/*!
	 * construct the problem (everything needed for solving)
	 */
	void construct_problem(formulation_mode mode);
	/*!
	 * optimize #adders or #full_adders within this loop
	 */
	void optimization_loop(formulation_mode mode);
    /*!
     * check if vector has only positv values after they were flipped in the constructor
     * the solver can then allow sign inversion for this vector
     * @param v current vector to check
     * @return true when vector is all positiv else false
     */
    bool vector_all_positive(const std::vector<int> v);
    /*!
     * @return C[0].size()
     */
    inline int c_row_size(){return C[0].size();}
    /*!
     * @return C.size()
     */
    inline int c_column_size(){return C.size();}
    /*!
     * C[0].size() idx are reserved by the inputs
     * @return idx input buffer caused by multiple inputs in SOP and CMM
     * not sure yet if the -1 is correct
     */
    inline int idx_input_buffer(){return C[0].size()-1;}
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
	void create_input_node_depth_variables(); // idx = 0 is the input node that has a constant value 0 as output
	void create_input_select_mux_variables(int idx); // select input variables
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
	void create_adder_depth_variables(int idx);
	void create_pipelining_variables(int idx);
	void create_output_stage_eq_variables();


	////////////////////////////////
	//// CREATE ALL CONSTRAINTS ////
	////////////////////////////////
	void create_input_output_constraints(formulation_mode mode);
	void create_input_select_constraints(int idx, formulation_mode mode);
	void create_input_select_limitation_constraints(int idx, formulation_mode mode);
	void create_shift_limitation_constraints(int idx, formulation_mode mode);
	void create_shift_constraints(int idx, formulation_mode mode);
	void create_negate_select_constraints(int idx, formulation_mode mode);
	void create_xor_constraints(int idx, formulation_mode mode);
	void create_adder_constraints(int idx, formulation_mode mode);
	void create_post_adder_shift_limitation_constraints(int idx, formulation_mode mode);
	void create_post_adder_shift_constraints(int idx, formulation_mode mode);
	void create_odd_fundamentals_constraints(int idx, formulation_mode mode);
	void create_mcm_output_constraints(formulation_mode mode);
	void create_mcm_input_constraints(formulation_mode mode);
	void create_full_adder_coeff_word_size_constraints(int idx, formulation_mode mode);
	void create_full_adder_msb_constraints(int idx, formulation_mode mode);
	void create_full_adder_coeff_word_size_sum_constraints(int idx, formulation_mode mode);
	void create_full_adder_shift_gain_constraints(int idx, formulation_mode mode);
	void create_full_adder_shift_sum_constraints(int idx, formulation_mode mode);
	void create_full_adder_msb_sum_constraints(formulation_mode mode);
	void create_full_adder_add_subtract_inputs_constraints(formulation_mode mode);
	void create_full_adder_cpa_constraints(formulation_mode mode);
	void create_full_adder_result_constraints();
	void create_adder_depth_computation_select_constraints(int idx, formulation_mode mode);
	void create_adder_depth_computation_max_constraints(int idx, formulation_mode mode);
	void create_adder_depth_computation_add_constraints(int idx, formulation_mode mode);
	void create_adder_depth_computation_limit_constraints(int idx, formulation_mode mode);
	void create_pipelining_register_constraints(int idx, formulation_mode mode);
	void create_pipelining_input_stage_equality_constraints(int idx, formulation_mode mode);
	void create_pipelining_output_stage_equality_constraints(int idx, formulation_mode mode);
	void prohibit_current_solution();

	/*!
	 * helper function to create a bitheap in SAT
	 * @param x with size [a][b_a] were a is the bit position and b_a is the bitheap height at position a
	 * @return the resulting bitvector after compression
	 */
	std::vector<int> create_bitheap(const std::vector<std::pair<std::vector<int>, bool>> &x);

	///////////////////////////////////
	//// INDICES FOR ALL VARIABLES ////
	//////////////////////////////////////////////////
	//// THE FIRST VARIABLE INDEX STARTS WITH 1!! ////
	//////////////////////////////////////////////////
	/*!
	 * < node idx, left/right, mux idx, bit > -> variable idx
	 * CMM Dimension
	 */
	std::map<std::tuple<int, input_direction, int, int, int>, int> input_select_mux_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 * /l(i) and r(i)
	 * CMM Dimension
	 */
	std::map<std::tuple<int, input_direction, int, int>, int> input_select_mux_output_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 * /alpha(i)
	 */
	std::map<std::tuple<int, input_direction, int>, int> input_select_selection_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * /gamma(i)
	 */
	std::map<std::pair<int, int>, int> input_shift_value_variables;
	/*!
	 * < node idx, mux stage, bit > -> variable idx
	 * /s(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int, int>, int> shift_internal_mux_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! identical to the last MUX stage of shift_internal_mux_output_variables
	 * /s(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int>, int> shift_output_variables;
	/*!
	 * node idx -> variable idx
	 * /delta(i)
	 */
	std::map<int, int> input_negate_select_variables;
	/*!
	 * < node idx, left/right, bit > -> variable idx
	 * /x(i) and w(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, input_direction, int, int>, int> negate_select_output_variables;
	/*!
	 * node idx -> variable idx
	 * /epsilon(i)
	 */
	std::map<int, int> input_negate_value_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * /y(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int>, int> xor_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * /C_in
	 * CMM Dimension
	 */
	std::map<std::tuple<int, int, int>, int> adder_carry_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * used for optimized adder clauses
	 * /no use
	 */
	std::map<std::pair<int, int>, int> adder_XOR_internal_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! node idx = 0 is the input node with constant value 0
	 * /z(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int>, int> adder_output_value_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * /zeta(i)
	 */
	std::map<std::pair<int, int>, int> input_post_adder_shift_value_variables;
	/*!
	 * < node idx, mux stage, bit > -> variable idx
	 * /c(i) mux
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int, int>, int> post_adder_shift_internal_mux_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! identical to the last MUX stage of shift_internal_mux_output_variables
	 * /c(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int>, int> post_adder_shift_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 * !!! node idx = 0 is the input node with constant value 0
	 * /z(i) or c(i)
	 * CMM dimension
	 */
	std::map<std::tuple<int, int, int>, int> output_value_variables;
	/*!
	 * < node idx, mcm constant > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> mcm_output_variables;
	/*!
	 * < node idx, bit > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> full_adder_coeff_word_size_variables;
	/*!
	 * < idx, stage, bit > -> variable idx
	 */
	std::map<std::tuple<int, int, int>, int> full_adder_coeff_word_size_internal_variables;
	/*!
	 * < idx, stage, bit > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> full_adder_coeff_word_size_internal_carry_input_variables;
	/*!
	 * < node idx > -> variable idx
	 */
	std::map<int, int> full_adder_msb_variables;
	/*!
	 * < idx, bit > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> full_adder_word_size_sum_variables;
	/*!
	 * < idx, bit > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> full_adder_shift_gain_variables;
	/*!
	 * < idx, bit > -> variable idx
	 */
	std::map<std::tuple<int, int>, int> full_adder_shift_sum_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_msb_sum_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_add_subtract_inputs_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_cpa_internal_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_result_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_comparator_ok_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> full_adder_comparator_carry_variables;
	/*!
	 * < idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> adder_depth_variables;
	/*!
	 * < idx, direction, bit > -> variable idx
	 * /eta and /theta
	 */
	std::map<std::tuple<int, input_direction, int>, int> adder_depth_computation_input_variables;
	/*!
	 * < idx, direction, mux idx, bit > -> variable idx
	 * /eta and /theta
	 */
	std::map<std::tuple<int, input_direction, int, int>, int> adder_depth_computation_input_mux_variables;
	/*!
	 * < idx, bit > -> variable idx
	 */
	std::map<std::pair<int, int>, int> adder_depth_computation_max_variables;
	/*!
	 * < idx > -> variable idx
	 */
	std::map<int, int> register_variables;
	/*!
	 * < idx > -> variable idx
	 */
	std::map<int, int> input_stages_equal_variables;
	/*!
	 * < bit > -> variable idx
	 */
	std::map<int, int> output_stage_eq_variables;
	/*!
	 * a variable that is forced to 1
	 */
	int const_one_bit = -1;
	/*!
	 * a variable that is forced to 0
	 */
	int const_zero_bit = -1;
	/*!
	 * init this->const_one_bit if not yet initialized
	 * @return this->const_one_bit after init
	 */
	int init_const_one_bit();
	/*!
	 * init this->const_one_bit if not yet initialized
	 * @return this->const_one_bit after init
	 */
	int init_const_zero_bit();
};


#endif //SATMCM_MCM_H
