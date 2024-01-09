//
// Created by nfiege on 05/19/2023.
//

#ifndef SATSCM_MCM_KISSAT_H
#define SATSCM_MCM_KISSAT_H

#ifdef USE_KISSAT

#include <mcm.h>
#include <kissatpp.hpp>
#include <memory>

class mcm_kissat : public mcm {

#define KISSAT_SAT 10
#define KISSAT_UNSAT 20

public:
    mcm_kissat(const std::vector<std::vector<int>> &C, int timeout, verbosity_mode verbosity, bool allow_negative_numbers, bool write_cnf);

protected:
    std::pair<bool, bool> check() override;
    void reset_backend(formulation_mode mode) override;
    int get_result_value(int var_idx) override;

    void create_arbitrary_clause(const std::vector<std::pair<int, bool>> &a) override;
	bool supports_incremental_solving() const override { return false; }

private:
    std::unique_ptr<kissatpp::kissatpp> solver;
};

#endif //USE_KISSAT

#endif //SATSCM_MCM_KISSAT_H
