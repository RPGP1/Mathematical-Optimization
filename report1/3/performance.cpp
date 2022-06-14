#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>


int main()
{
    constexpr int rows = 3, cols = 10;
    constexpr double lambda = 1;

    const auto A = Eigen::Matrix<double, rows, cols>::Random().eval();
    const auto b = Eigen::Vector<double, rows>::Random().eval();
    const auto w0 = (Eigen::Vector<double, cols>::Random() * 8192.0).eval();


    constexpr auto NumIter = 1 << 16;
    std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1);


    {
        const auto fixed_step_benchmark = [&] {
            const auto AT_A = (A.transpose() * A).eval();
            const auto AT_b = (A.transpose() * b).eval();
            const double AT_A_max_eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, cols, cols>>{AT_A}.eigenvalues().maxCoeff();

            const double L = 2.0 * (AT_A_max_eigenvalue + lambda),
                         L_inv = 1.0 / L;

            auto w = w0;
            double prev_f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
            do {
                for (auto i = 0; i < 10; i++) {
                    w += (-L_inv * 2.0 * (AT_A * w + lambda * w - AT_b)).eval();
                }

                const auto f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
                if (prev_f_w - f_w < 1e-13) {
                    break;
                }
                prev_f_w = f_w;
            } while (true);
        };

        for (auto i = 0; i < 3; i++) {
            fixed_step_benchmark();
        }
        std::chrono::high_resolution_clock::duration duration{};
        for (auto i = 0; i < NumIter; i++) {
            const auto begin_time = std::chrono::high_resolution_clock::now();
            fixed_step_benchmark();
            duration += std::chrono::high_resolution_clock::now() - begin_time;
        }
        std::cout << "fixed_step\t" << std::chrono::duration_cast<std::chrono::duration<long double, std::nano>>(duration).count() / NumIter << std::endl;
    }


    {
        const auto benchmark = [&] {
            const auto AT_A = (A.transpose() * A).eval();
            const auto AT_b = (A.transpose() * b).eval();
            const double AT_A_max_eigenvalue = Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, cols, cols>>{AT_A}.eigenvalues().maxCoeff();

            const double L = 2.0 * (AT_A_max_eigenvalue + lambda),
                         L_inv = 1.0 / L;

            auto w = w0, v = w;
            int iter = 0;
            double prev_f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();

            do {
                for (auto i = 0; i < 10; i++) {
                    const auto prev_w = w;

                    w = v - L_inv * 2.0 * (AT_A * v + lambda * v - AT_b);
                    v = w + (prev_w - w) * iter / (iter + 3);

                    iter++;
                }

                const auto f_w = (A * w - b).squaredNorm() + lambda * w.squaredNorm();
                if (prev_f_w - f_w < 1e-13) {
                    break;
                }
                iter++;
                prev_f_w = f_w;
            } while (true);
        };

        for (auto i = 0; i < 3; i++) {
            benchmark();
        }
        std::chrono::high_resolution_clock::duration duration{};
        for (auto i = 0; i < NumIter; i++) {
            const auto begin_time = std::chrono::high_resolution_clock::now();
            benchmark();
            duration += std::chrono::high_resolution_clock::now() - begin_time;
        }
        std::cout << "Nesterov_acc\t" << std::chrono::duration_cast<std::chrono::duration<long double, std::nano>>(duration).count() / NumIter << std::endl;
    }


    return 0;
}
