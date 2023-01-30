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
    cout << "\n" << endl;


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
    cout << "\n" << endl;
    cout << "All Scores = \n" << m_y << endl;
    cout << "\n" << endl;


    m_x_max = data_points.col(best_index);

    // Determine the next pairwise comparison query
    m_optimizer->DetermineNextQuery();
}

double Core::evaluateObjectiveFunction(const Eigen::VectorXd& x) const
{
    return 1.0 - 1.5 * x(0) * std::sin(x(0) * 13.0);
}

std::string parameter(int argc, char* argv[]){

    std::string scenario;

    if (argc == 1){

        scenario = "../../../csv/Corridor_Scenario.csv";
#ifdef _WIN32
        scenario = "../../../../../csv/Corridor_Scenario.csv";
#endif

    } else if (argc > 2){

        std::string csv_arg = argv[1];

        if (!csv_arg.compare("--csv")){

            scenario = argv[2];

        } else{
            std::cout << "Please use '--csv' to specify the absolute csv file path.\n" << std::endl;
            exit( 1 );
        }



    } else{
        std::cout << "NO Parameter!\n" << std::endl;
        exit( 1 );

    }


    return scenario;

}

int main(int argc, char* argv[]){


    int trials = 1;

    std::ifstream csvFile;
    std::string strPathCSVFile = parameter(argc, argv);

    std::cout << "The parameters in scenario:\n" << strPathCSVFile << "\n is chosen to be evaluated." << std::endl;

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
    unsigned n_iterations = (unsigned)pairs.size();



    for (int t =0 ; t< trials; t++){

                std::cout << "========================" << std::endl;
                std::cout << "Trial " << t + 1 << std::endl;
                std::cout << "========================" << std::endl;

        for (unsigned i = 0; i < n_iterations; ++i)
        {
            std::cout << "-------------All Combinations Candidates---------" << std::endl;

            core.proceedOptimization(total, pairs, i);

            std::vector<Eigen::VectorXd> tt = core.m_optimizer->GetCurrentOptions();


            std::cout << "The current best hyper parameter combination is: \n" << core.m_x_max << std::endl;
            std::cout << "\n" << std::endl;

            std::cout << "The score for current best hyper parameter combination is: \n" << setprecision(12) << core.m_y_max << std::endl;
            std::cout << "\n" << std::endl;


        }
    }

    return 0;
}