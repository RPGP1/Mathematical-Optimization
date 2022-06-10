#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <cmdline.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>

template <class T, size_t N, class Func>
std::array<std::invoke_result_t<Func, T>, N> transform_array(std::array<T, N> input, Func&& f)
{
    std::array<std::invoke_result_t<Func, T>, N> tmp;
    std::transform(input.cbegin(), input.cend(), tmp.begin(), std::forward<Func>(f));
    return tmp;
}

int main(int argc, char* argv[])
{
    cmdline::parser parser;
    parser.add<int>("rows", 'm');
    parser.add<int>("cols", 'n');

    if (!parser.parse(argc, argv)) {
        std::cout << parser.error_full() << std::endl
                  << parser.usage() << std::endl;
        return EXIT_FAILURE;
    }

    const int cols = parser.get<int>("cols"), rows = parser.get<int>("rows");

    const Eigen::MatrixXd A = Eigen::MatrixXd::Random(rows, cols);
    const Eigen::VectorXd b = Eigen::VectorXd::Random(rows),
                          w0 = Eigen::VectorXd::Random(cols) * 8192.0;


    const Eigen::MatrixXd AT_A = A.transpose() * A;
    const Eigen::VectorXd AT_b = A.transpose() * b;
    const double AT_A_max_eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>{AT_A}.eigenvalues().maxCoeff();


    const std::array<double, 3> lambda_list{0.0, 1.0, 10.0};

    const auto min_list = transform_array(lambda_list, [&](double lambda) {
        const Eigen::VectorXd argmin = (AT_A + lambda * Eigen::MatrixXd::Identity(cols, cols)).fullPivHouseholderQr().solve(AT_b);
        return b.squaredNorm() - argmin.dot(AT_b);
    });


    constexpr auto NumIter = 200;

    const auto logs = transform_array(lambda_list, [&](double lambda) {
        std::vector<double> log;
        log.reserve(NumIter);

        const double L = 2.0 * (AT_A_max_eigenvalue + lambda),
                     L_inv = 1.0 / L;

        auto w = w0;
        for (auto i = decltype(NumIter){0}; i < NumIter; i++) {
            const Eigen::VectorXd df_w = 2.0 * (AT_A * w + lambda * w - AT_b);
            w += -L_inv * df_w;

            const auto f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
            log.push_back(f_w);
        }

        return log;
    });


    std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    for (auto iter_idx = decltype(NumIter){0}; iter_idx < NumIter + 1; iter_idx++) {
        std::cout << iter_idx + 1;
        for (unsigned log_idx = 0; log_idx < lambda_list.size(); log_idx++) {
            std::cout << '\t' << logs[log_idx][iter_idx] - min_list[log_idx];
        }
        std::cout << std::endl;
    }

    return 0;
}
