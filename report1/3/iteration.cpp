#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cmdline.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>


int main(int argc, char* argv[])
{
    cmdline::parser parser;
    parser.add<int>("rows", 'm');
    parser.add<int>("cols", 'n');
    parser.add<double>("lambda");

    if (!parser.parse(argc, argv)) {
        std::cerr << parser.error_full() << std::endl
                  << parser.usage() << std::endl;
        return EXIT_FAILURE;
    }

    const int cols = parser.get<int>("cols"), rows = parser.get<int>("rows");
    const double lambda = parser.get<double>("lambda");

    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(rows, cols);
    const Eigen::VectorXd b = Eigen::VectorXd::Random(rows),
                          w0 = Eigen::VectorXd::Random(cols) * 8192.0;


    const Eigen::MatrixXd AT_A = A.transpose() * A;
    const Eigen::VectorXd AT_b = A.transpose() * b;
    const double AT_A_max_eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>{AT_A}.eigenvalues().maxCoeff();


    const auto min = b.squaredNorm()
                     - (AT_A + lambda * Eigen::MatrixXd::Identity(cols, cols))
                           .fullPivHouseholderQr()
                           .solve(AT_b)
                           .dot(AT_b);


    constexpr auto NumIter = 500;

    std::vector<double> log_fixed_step;
    {
        log_fixed_step.reserve(NumIter);

        const double L = 2.0 * (AT_A_max_eigenvalue + lambda),
                     L_inv = 1.0 / L;

        auto w = w0;
        for (auto i = decltype(NumIter){0}; i < NumIter; i++) {
            const Eigen::VectorXd df_w = 2.0 * (AT_A * w + lambda * w - AT_b);
            w += -L_inv * df_w;

            const auto f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
            log_fixed_step.push_back(f_w);
        }
    }

    std::vector<double> log_Nesterov;
    {
        log_Nesterov.reserve(NumIter);

        const double L = 2.0 * (AT_A_max_eigenvalue + lambda),
                     L_inv = 1.0 / L;

        auto w = w0, v = w;

        for (auto i = decltype(NumIter){0}; i < NumIter; i++) {
            const auto prev_w = w;

            const Eigen::VectorXd df_v = 2.0 * (AT_A * v + lambda * v - AT_b);
            w = v - L_inv * df_v;
            v = w + (prev_w - w) * i / (i + 3);

            const double f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
            log_Nesterov.push_back(f_w);
        }
    }


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (auto iter_idx = decltype(NumIter){0}; iter_idx < NumIter + 1; iter_idx++) {
        std::cout << iter_idx + 1 << '\t'
                  << log_fixed_step[iter_idx] - min << '\t'
                  << log_Nesterov[iter_idx] - min << std::endl;
    }

    return 0;
}
