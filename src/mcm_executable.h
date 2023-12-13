//
// Created by nfiege on 12/7/23.
//

#ifndef SATSCM_MCM_EXECUTABLE_H
#define SATSCM_MCM_EXECUTABLE_H

#include <mcm.h>
#include <cstdlib>

class mcm_executable : public mcm {
public:
    mcm_executable(const std::vector<std::vector<int>> &C, const std::string &new_path_to_executable, const std::string &new_pre_cnf_file_args, const std::string &new_post_cnf_file_args, const std::string &solver_log_filename, const std::string &solver_err_filename, verbosity_mode verbosity, bool allow_negative_numbers, bool write_cnf);

protected:
    std::pair<bool, bool> check() override;
    void reset_backend(formulation_mode mode) override;
    int get_result_value(int var_idx) override;
    bool needs_cnf_generation() const override { return true; }

private:
    std::string path_to_executable;
    std::string pre_cnf_file_args;
    std::string post_cnf_file_args;
    std::string solver_log_filename;
    std::string solver_err_filename;
    std::map<int, int> result_values;
};


#endif //SATSCM_MCM_EXECUTABLE_H
