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

    parser.add<double>("step0");
    parser.add<double>("step_ratio");
    parser.add<double>("threshold");

    if (!parser.parse(argc, argv)) {
        std::cerr << parser.error_full() << std::endl
                  << parser.usage() << std::endl;
        return EXIT_FAILURE;
    }

    const int cols = parser.get<int>("cols"), rows = parser.get<int>("rows");
    const double lambda = parser.get<double>("lambda");

    const double step0 = parser.get<double>("step0"), step_ratio = parser.get<double>("step_ratio"), threshold = parser.get<double>("threshold");

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


    constexpr auto NumIter = 200;

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

    std::vector<double> log_Armijo;
    {
        log_Armijo.reserve(NumIter);

        auto w = w0;
        double f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();

        for (auto i = decltype(NumIter){0}; i < NumIter; i++) {
            const Eigen::VectorXd df_w = 2.0 * (AT_A * w + lambda * w - AT_b);

            double threshold_offset = -threshold * step0 * df_w.squaredNorm();
            Eigen::VectorXd step = -step0 * df_w;

            Eigen::VectorXd tmp_w = w + step;
            while ((A * tmp_w - b).squaredNorm() + lambda * tmp_w.squaredNorm() > f_w + threshold_offset) {
                threshold_offset *= step_ratio;
                step *= step_ratio;
                tmp_w = w + step;
            }
            w = tmp_w;

            f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
            log_Armijo.push_back(f_w);
        }
    }


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (auto iter_idx = decltype(NumIter){0}; iter_idx < NumIter + 1; iter_idx++) {
        std::cout << iter_idx + 1 << '\t'
                  << log_fixed_step[iter_idx] - min << '\t'
                  << log_Armijo[iter_idx] - min << std::endl;
    }

    return 0;
}
