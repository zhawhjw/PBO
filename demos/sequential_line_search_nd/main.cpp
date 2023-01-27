#include "core.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <mathtoolbox/data-normalization.hpp>
#include <sequential-line-search/preferential-bayesian-optimizer.hpp>

//#include <iostream>
//#include <sequential-line-search/sequential-line-search.hpp>
//#include <sequential-line-search/utils.hpp>
#ifdef PARALLEL
#include <parallel-util.hpp>
#endif
#include <timer.hpp>
#include <fstream>

using namespace std;

using namespace sequential_line_search;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Core::Core()
{
    reset();
}

void Core::reset()
{
    // Define a generator that generates the initial options
    auto generator = [](const int num_dims, const int) {
        return std::vector<VectorXd>{VectorXd::Constant(num_dims, 0.3), VectorXd::Constant(num_dims, 0.6)};
    };



    constexpr bool   use_map      = false;
    constexpr double signal_var   = 0.5;
    constexpr double length_scale = 0.5;
    constexpr double noise_level  = 0.0;
    constexpr double prior        = 0.0;

    // Note: This app is for a one-dimensional problem.
    constexpr int num_dims = 6;

    // Note: This app uses pairwise comparison.
    constexpr int num_options = 2;

    m_optimizer = std::make_shared<PreferentialBayesianOptimizer>(num_dims,
                                                                  use_map,
                                                                  KernelType::ArdMatern52Kernel,
                                                                  AcquisitionFuncType::ExpectedImprovement,
                                                                  generator,
                                                                  CurrentBestSelectionStrategy::LastSelection,
                                                                  num_options);
    m_optimizer->SetHyperparams(signal_var, length_scale, noise_level, prior);

    m_y = VectorXd::Zero(6);

    m_x_max = VectorXd::Zero(6);
    m_y_max = NAN;
}

void Core::proceedOptimization(vector<vector<double>> t, vector<vector<int>> p, int iteration)
{

    vector<double> PBD_On_1_Off_0 = t[0];
    vector<double> alpha = t[1];
    vector<double> w_smooth_acceleration = t[2];
    vector<double> w_col = t[3];
    vector<double> w_comfort = t[4];
    vector<double> Cheby_1_Or_Euclidean_2 = t[5];

    vector<double> liked = t[6];

    vector<int> cp = p[iteration];

    VectorXd best;
    vector<VectorXd> others;
    others.clear();

    int best_option_index = -1;

    double max = 0.0;
    for (int k =0; k<cp.size();k++){

        int index = cp[k];

        if(liked[index] > max){
            max = liked[index] ;
            best_option_index = index;
        }

    }

    for (int k =0; k<cp.size();k++){

        int index = cp[k];

        VectorXd content(6);
        content << PBD_On_1_Off_0[index], alpha[index], w_smooth_acceleration[index], w_col[index], w_comfort[index], Cheby_1_Or_Euclidean_2[index];

        if(index == best_option_index){
            best = content;
        }
        else{
            others.push_back(content);
        }

    }






    m_optimizer->SubmitCustomFeedbackData(best, others);

    // Update internal data according to the new preference model
    const auto data_points     = m_optimizer->GetRawDataPoints();
    cout << data_points << endl;

    const int  num_data_points = data_points.cols();

    VectorXd f(num_data_points);
    for (int i = 0; i < num_data_points; ++i)
    {
        f(i) = m_optimizer->GetPreferenceValueMean(data_points.col(i));
    }

    m_normalizer = std::make_shared<mathtoolbox::DataNormalizer>(f.transpose());
    m_y          = VectorXd::Constant(f.size(), 1.0) + m_normalizer->GetNormalizedDataPoints().transpose();

    int best_index;
    m_y_max = m_y.maxCoeff(&best_index);
    cout << "Best Index = " << best_index << endl;
    cout << "All Y values = " << m_y << endl;


    m_x_max = data_points.col(best_index);

    // Determine the next pairwise comparison query
    m_optimizer->DetermineNextQuery();
}

double Core::evaluateObjectiveFunction(const Eigen::VectorXd& x) const
{
    return 1.0 - 1.5 * x(0) * std::sin(x(0) * 13.0);
}

int main(int argc, char* argv[]){

    std::ifstream csvFile;
    std::string strPathCSVFile = "../../../csv/Corridor_Scenario.csv";
#ifdef _WIN32
    strPathCSVFile = "../../../../../csv/Corridor_Scenario.csv";
#endif

    std::cout << strPathCSVFile << std::endl;

    csvFile.open(strPathCSVFile.c_str());

    if (!csvFile.is_open())
    {
        std::cout << "Path Wrong!!!!" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<double> PBD_On_1_Off_0;
    std::vector<double> alpha;
    std::vector<double> w_smooth_acceleration;
    std::vector<double> w_col;

    std::vector<double> w_comfort;
    std::vector<double> Cheby_1_Or_Euclidean_2;
    std::vector<double> liked;

    std::vector<vector<double>> total;

    vector<int> first = {0, 1} ;
    vector<int> second = {2, 3, 4};
    vector<int> third = {5, 6, 7};
    vector<int> forth = {8, 9, 10};
    vector<int> fifth = {11, 12, 13};
    vector<int> sixth = {14, 15};

    vector<vector<int>> pairs;

    pairs.push_back(first);
    pairs.push_back(second);
    pairs.push_back(third);
    pairs.push_back(forth);
    pairs.push_back(fifth);
    pairs.push_back(sixth);


    std::string line;
    std::vector <std::string > vec;
    getline(csvFile, line); // skip the 1st line

    while (getline(csvFile,line))
    {
        if (line.empty()) // skip empty lines:
        {
            //cout << "empty line!" << endl;
            continue;
        }

        std::istringstream iss(line);
        std::string lineStream;
        std::string::size_type sz;

        std::vector <double> row;

        while (getline(iss, lineStream, ','))
        {
            row.push_back(stod(lineStream,&sz)); // convert to double
        }

        PBD_On_1_Off_0.push_back(row[0]);

        alpha.push_back(row[1]);
        w_smooth_acceleration.push_back(row[2]);
        w_col.push_back(row[3]);

        w_comfort.push_back(row[4]);
        Cheby_1_Or_Euclidean_2.push_back(row[5]);
        liked.push_back(row[6]);
    }

    total.push_back(PBD_On_1_Off_0);
    total.push_back(alpha);
    total.push_back(w_smooth_acceleration);
    total.push_back(w_col);
    total.push_back(w_comfort);
    total.push_back(Cheby_1_Or_Euclidean_2);

    total.push_back(liked);

    Core& core = Core::getInstance();
    core.reset();
    constexpr unsigned n_iterations = 6;

//    std::vector<VectorXd> xx = std::vector<VectorXd>{VectorXd::Constant(6, 0.4), VectorXd::Constant(6, 0.6)};
//    std::cout << xx[0].size() << std::endl;
    //    for (VectorXd content: xx)
//
//        std::cout << content << std::endl;

    int trials = 1;
    for (int t =0 ; t< trials; t++){

                std::cout << "========================" << std::endl;
                std::cout << "Trial " << t + 1 << std::endl;
                std::cout << "========================" << std::endl;

        for (unsigned i = 0; i < n_iterations; ++i)
        {
            core.proceedOptimization(total, pairs, i);

            std::vector<Eigen::VectorXd> tt = core.m_optimizer->GetCurrentOptions();
            std::cout << "m_x_max = " << core.m_x_max << std::endl;
            std::cout << "m_y_max = " << setprecision(12) << core.m_y_max << std::endl;
            std::cout << "---------------------------------------------" << std::endl;



        }
    }

    return 0;
}







//
//#include <iostream>
//#include <sequential-line-search/sequential-line-search.hpp>
//#include <sequential-line-search/utils.hpp>
//#ifdef PARALLEL
//#include <parallel-util.hpp>
//#endif
//#include <timer.hpp>
//#include <stdio.h>
//#include <fstream>
//
////Found maximizer: 0.626 0.836 0.308 0.221 0.198 0.612
////    Found maximum: 0.692362
//
////Found maximizer:  0.527  0.848   0.25  0.256 0.0732  0.514
////    Found maximum: 0.684029
//
////Found maximizer:   0.599       1   0.144   0.187 5.6e-08   0.625
////    Found maximum: 0.486117
//
//using namespace std;
//
//namespace
//{
//    constexpr double a         = 0.500;
//    constexpr double r         = 0.500;
//    constexpr double b         = 0.001;
//    constexpr double variance  = 0.100;
//    constexpr double btl_scale = 0.010;
//
//    constexpr unsigned n_trials     = 8;
//    constexpr unsigned n_iterations = 16;
//
//    constexpr bool use_slider_enlargement = true;
//    constexpr bool use_MAP_hyperparams    = true;
//
//    constexpr unsigned test_dimension = 6;
//
//    // Define a test function
//    double evaluateObjectiveFunction(const Eigen::VectorXd& x)
//    {
//        assert(x.rows() == test_dimension);
//
//
////        std::cout << "Here is the vector v:\n" << x << std::endl;
////        std::cout << "Here is the data v:\n" << *x.data() << std::endl;
//
//
//        auto lambda = [](const Eigen::VectorXd& x, const Eigen::VectorXd& mu, const double sigma)
//        { return std::exp(-(x - mu).squaredNorm() / (sigma * sigma)); };
//
//        return lambda(x, Eigen::VectorXd::Constant(test_dimension, 0.4), 1.0);
//    }
//
//    // Define an analytic solution of the test function
//    const Eigen::VectorXd analytic_solution = Eigen::VectorXd::Constant(test_dimension, 0.4);
//} // namespace
//
//
//void readCSV(const string &strPath2Dataset, const string &csvName)
//{
//    std::ifstream csvFile;
//    std::string strPathCSVFile = strPath2Dataset + csvName;
//    std::cout << strPathCSVFile << std::endl;
//
//    csvFile.open(strPathCSVFile.c_str());
//
//    if (!csvFile.is_open())
//    {
//        std::cout << "Path Wrong!!!!" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//    std::vector<long double> PBD_On_1_Off_0;
//    std::vector<long double> alpha;
//    std::vector<long double> w_smooth_acceleration;
//    std::vector<long double> w_col;
//
//    std::vector<long double> w_comfort;
//    std::vector<long double> Cheby_1_Or_Euclidean_2;
//    std::vector<long double> liked;
//
//    std::string line;
//    std::vector <std::string > vec;
//    getline(csvFile, line); // skip the 1st line
//
//    while (getline(csvFile,line))
//    {
//        if (line.empty()) // skip empty lines:
//        {
//            //cout << "empty line!" << endl;
//            continue;
//        }
//
//        std::istringstream iss(line);
//        std::string lineStream;
//        std::string::size_type sz;
//
//        std::vector <long double> row;
//
//        while (getline(iss, lineStream, ','))
//        {
//            row.push_back(stold(lineStream,&sz)); // convert to double
//        }
//
//        PBD_On_1_Off_0.push_back(row[0]);
//
//        alpha.push_back(row[1]);
//        w_smooth_acceleration.push_back(row[2]);
//        w_col.push_back(row[3]);
//
//        w_comfort.push_back(row[4]);
//        Cheby_1_Or_Euclidean_2.push_back(row[5]);
//        liked.push_back(row[6]);
//    }
//
//    //cout << "size ts = " << timeStampIMU.size() << endl;
//    for (size_t i = 0; i < PBD_On_1_Off_0.size(); i++)
//    {
//        std::cout << "PBD_On_1_Off_0 = " << std::setprecision(12) << PBD_On_1_Off_0[i] << std::endl;
//
//        std::cout << "alpha = " << std::setprecision(12) << alpha[i] << std::endl;
//        std::cout << "w_smooth_acceleration = " << std::setprecision(12) << w_smooth_acceleration[i] << std::endl;
//        std::cout << "w_col = " << std::setprecision(12) << w_col[i] << std::endl;
//
//        std::cout << "w_comfort = " << std::setprecision(12) << w_comfort[i] << std::endl;
//        std::cout << "Cheby_1_Or_Euclidean_2 = " << std::setprecision(12) << Cheby_1_Or_Euclidean_2[i] << std::endl;
//        std::cout << "liked = " << std::setprecision(12) << liked[i] << std::endl;
//        std::cout << "--------------------------------" << std::endl;
//    }
//}
//
//
//int main(int argc, char* argv[])
//{
//    // Storage for performance reports
//    Eigen::MatrixXd objective_values(n_iterations, n_trials);
//    Eigen::MatrixXd residual_norms(n_iterations, n_trials);
//    Eigen::MatrixXd elapsed_times(n_iterations, n_trials);
//
////    readCSV("../../../csv/", "Corridor_Scenario.csv");
//
//    std::ifstream csvFile;
//    std::string strPathCSVFile = "../../../csv/Hallway_Scenario.csv";
//    std::cout << strPathCSVFile << std::endl;
//
//    csvFile.open(strPathCSVFile.c_str());
//
//    if (!csvFile.is_open())
//    {
//        std::cout << "Path Wrong!!!!" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//    std::vector<double> PBD_On_1_Off_0;
//    std::vector<double> alpha;
//    std::vector<double> w_smooth_acceleration;
//    std::vector<double> w_col;
//
//    std::vector<double> w_comfort;
//    std::vector<double> Cheby_1_Or_Euclidean_2;
//    std::vector<double> liked;
//
//    std::string line;
//    std::vector <std::string > vec;
//    getline(csvFile, line); // skip the 1st line
//
//    while (getline(csvFile,line))
//    {
//        if (line.empty()) // skip empty lines:
//        {
//            //cout << "empty line!" << endl;
//            continue;
//        }
//
//        std::istringstream iss(line);
//        std::string lineStream;
//        std::string::size_type sz;
//
//        std::vector <double> row;
//
//        while (getline(iss, lineStream, ','))
//        {
//            row.push_back(stod(lineStream,&sz)); // convert to double
//        }
//
//        PBD_On_1_Off_0.push_back(row[0]);
//
//        alpha.push_back(row[1]);
//        w_smooth_acceleration.push_back(row[2]);
//        w_col.push_back(row[3]);
//
//        w_comfort.push_back(row[4]);
//        Cheby_1_Or_Euclidean_2.push_back(row[5]);
//        liked.push_back(row[6]);
//    }
//
//    const auto perform_test = [&](const int trial_index)
//    {
//        sequential_line_search::SequentialLineSearchOptimizer optimizer(
//            test_dimension, use_slider_enlargement, use_MAP_hyperparams);
//
//        optimizer.SetHyperparams(a, r, b, variance, btl_scale);
//
//        std::cout << "========================" << std::endl;
//        std::cout << "Trial " << trial_index + 1 << std::endl;
//        std::cout << "========================" << std::endl;
//
//
//
//        // Iterate optimization steps
//        for (unsigned i = 0; i < n_iterations; ++i)
//        {
//            std::cout << "---- Iteration " << i + 1 << " ----" << std::endl;
//
//            constexpr double search_epsilon = 1e-05;
//
//            // Search the best position in the current slider space
//            double max_slider_position = 0.0;
//            double max_y               = -1e+10;
//            for (double slider_position = 0.0; slider_position <= 1.0; slider_position += search_epsilon)
//            {
//                Eigen::VectorXd query = optimizer.CalcPointFromSliderPosition(slider_position);
//
////                query[0] = PBD_On_1_Off_0[11+i];
////                query[1] = alpha[11+i];
////                query[2] = w_smooth_acceleration[11+i];
////                query[3] = w_col[11+i];
////
////                query[4] = w_comfort[11+i];
////                query[5] = Cheby_1_Or_Euclidean_2[11+i];
////
//                query[0] = PBD_On_1_Off_0[i];
//                query[1] = alpha[i];
//                query[2] = w_smooth_acceleration[i];
//                query[3] = w_col[i];
//
//                query[4] = w_comfort[i];
//                query[5] = Cheby_1_Or_Euclidean_2[i];
//
//
//
//                const double y = evaluateObjectiveFunction(query);
//                if (y > max_y)
//                {
//                    max_y               = y;
//                    max_slider_position = slider_position;
//                }
//            }
//
//            const Eigen::VectorXd max_x = optimizer.CalcPointFromSliderPosition(max_slider_position);
//
//            std::cout << "x: " << max_x.transpose().format(Eigen::IOFormat(3)) << std::endl;
//            std::cout << "y: " << max_y << std::endl;
//
//            objective_values(i, trial_index) = max_y;
//            residual_norms(i, trial_index)   = (max_x - analytic_solution).norm();
//
//            timer::Timer t;
//
//            // Perform Bayesian optimization
//            optimizer.SubmitFeedbackData(max_slider_position);
//
//            elapsed_times(i, trial_index) = t.get_elapsed_time_in_milliseconds();
//        }
//
//        const Eigen::VectorXd x_star = optimizer.GetMaximizer();
//        const double          y_star = evaluateObjectiveFunction(x_star);
//
//        std::cout << std::endl;
//        std::cout << "Found maximizer: " << x_star.transpose().format(Eigen::IOFormat(3)) << std::endl;
//        std::cout << "Found maximum: " << y_star << std::endl << std::endl;
//    };
//
//#ifdef PARALLEL
//    parallelutil::queue_based_parallel_for(n_trials, perform_test);
//#else
//    for (unsigned trial_index = 0; trial_index < n_trials; ++trial_index)
//    {
//        perform_test(trial_index);
//    }
//#endif
//
//    // Export a report as a CSV file
//    sequential_line_search::utils::ExportMatrixToCsv("objective_values.csv", objective_values);
//    sequential_line_search::utils::ExportMatrixToCsv("residual_norms.csv", residual_norms);
//    sequential_line_search::utils::ExportMatrixToCsv("elapsed_times.csv", elapsed_times);
//
//    return 0;
//}
