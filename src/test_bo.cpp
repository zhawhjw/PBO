//
// Created by fantasticapple on 12/13/22.
//
#include <iostream>
#include <sequential-line-search/sequential-line-search.hpp>

using Eigen::VectorXd;
using sequential_line_search::SequentialLineSearchOptimizer;

double AskHumanForSliderManipulation(const std::pair<VectorXd, VectorXd>& slider_ends)
{
    // ...
    // ...

    return (double)0.01;
}

int main()
{
    // Instantiate an optimizer
    constexpr int dim = 6;
    SequentialLineSearchOptimizer optimizer(dim);

    // Iterate optimization steps
    constexpr int num_iters = 1;
    for (int i = 0; i < num_iters; ++i)
    {
        // Retrieve a slider space
        const std::pair<VectorXd, VectorXd> slider_ends = optimizer.GetSliderEnds();

        // Query slider manipulation
        const double slider_position = AskHumanForSliderManipulation(slider_ends);

        // Feed the slider manipulation result to the optimizer
//        optimizer.SubmitLineSearchResult(slider_position);
    }

    // Display the found solution
    std::cout << optimizer.GetMaximizer() << std::endl;

    return 0;
}