/*
 * MIT License
 *
 * Copyright (c) 2024 Björnborg Nguyen
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "gnuplot-iostream.h"

#include "reaching.hpp"
#include "tinypso.hpp"

#include <format>
#include <fstream>


TEST_CASE("PSO Curve fitting on normal distributions")
{

  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution d{0.0, 0.1};

  Eigen::MatrixXd solution(3, 3);
  solution << 2.0, 4.0, -3.0, 0.20, 0.30, 0.7, 0.00159572, 0.00159706, 0.01;

  Eigen::ArrayXd const x = Eigen::ArrayXd::LinSpaced(1001, 0.0, 1.0);
  Eigen::ArrayXd noise = x;

  for (auto &&n : noise) {
    n = d(gen);
  }
  Eigen::ArrayXd const yClean = superpositionGaussians(x, solution);
  Eigen::ArrayXd const yNoisy = yClean + noise;

  auto evaluate = [&x, &yNoisy](Particle const &a_particle) -> double {
    Eigen::MatrixXd posMat{a_particle.getPosition()};
    if (posMat.rows() < 3 && posMat.cols() > 0) {
      std::cerr << "Dimensionality of particle init is incorrect." << std::endl;
      return 0.0;
    }
    // Reconstructed fit
    auto const y = superpositionGaussians(x, posMat);

    // Fitness, the larger the better
    // double const fitness{1.0 / (std::sqrt(posMat.cols() * (yRef -
    // y).square().mean()))};
    double const fitness{
        1.0 / (posMat.cols() * std::sqrt((yNoisy - y).square().mean()))};

    return fitness;
  };
  Eigen::Index optimumCounter;
  yNoisy.abs().matrix().maxCoeff(&optimumCounter);

  std::shared_ptr<PsoSettings> psoSettings = std::make_shared<PsoSettings>();
  psoSettings->objectiveFunction = evaluate;
  psoSettings->particleRow = 3; //  dimensionality
  psoSettings->particleCol = 1; // initial number of bump guesses
  psoSettings->swarmSize = 10000;
  psoSettings->posMin =
      Eigen::ArrayXd{{yNoisy.minCoeff(), x.minCoeff(), 0.00001}};
  psoSettings->posMax = Eigen::ArrayXd{{yNoisy.maxCoeff(), x.maxCoeff(), 0.4}};
  // psoSettings->initialGuess = Eigen::ArrayXd{{yRef(optimumCounter),
  // x(optimumCounter), 0.2}}; // not used
  psoSettings->alpha = 1.0;
  psoSettings->dt = 1.0;
  psoSettings->speedMax = 1.1;
  psoSettings->maxInertia = 1.4; // exploration vs exploitation
  psoSettings->minInertia = 0.3;
  psoSettings->cognitiveFactor = 3;
  psoSettings->inertiaDecay = 0.995; // transforms exploration -> exploitation
  psoSettings->socialFactor = 4.0 - psoSettings->cognitiveFactor;
  psoSettings->numSteps = 200;
  psoSettings->threads = 16;
  psoSettings->convergeIterations = 30;

  ParticleSwarmOptimization pso(psoSettings);

  Gnuplot gnuplot("");
  gnuplot << "set term wxt 1 noraise\n";
  gnuplot << "set title \"Pso curvefitting\"\n";
  gnuplot << "set grid\n";
  std::tuple<Eigen::ArrayXd, Eigen::ArrayXd> gnudataNoisy =
      std::forward_as_tuple(x, yNoisy);
  std::tuple<Eigen::ArrayXd, Eigen::ArrayXd> gnudataClean =
      std::forward_as_tuple(x, yClean);
  gnuplot << "plot"
          << gnuplot.binFile1d(gnudataNoisy, "record", "noisyreference.dat")
          << "with lines title 'Noisy ref'"
          << ","
          << gnuplot.binFile1d(gnudataClean, "record", "clearnreference.dat")
          << "with lines title 'Clean ref'"
          << "," << gnuplot.binFile1d(gnudataNoisy, "record", "fit.dat")
          << "with lines title 'Fit'" << std::endl;


  auto plot{[&gnuplot, &gnudataNoisy, &x](Particle const &a_particle) -> bool {
    auto const fitness = a_particle.getBestFitness();

    std::get<1>(gnudataNoisy) = superpositionGaussians(x, std::get<1>(fitness));
    gnuplot.binFile1d(gnudataNoisy, "record", "fit.dat");
    gnuplot << "replot" << std::endl;
    return true;
  }};
  psoSettings->plotFunction = plot;
  // Run the optimization and get the results
  Particle const particle = pso.runOptimization();
  auto const fitness = particle.getBestFitness();

  plot(particle);

  std::clog << "Best Particle fitness  = " << std::get<0>(fitness) << "\n"
            << std::get<1>(fitness) << std::endl;
  // Solution
  // 2023-09-04 10:37:14 bb | Individual data,
  // 0: amplitude
  // 1: timeshift (statistical mean)
  // 2: variance

  Eigen::MatrixXd proposal = std::get<1>(fitness);
  Eigen::ArrayXd const fit = superpositionGaussians(x, std::get<1>(fitness));

  double const error{(solution - proposal).cwiseAbs().sum()};
  // Overfitting
  REQUIRE(solution.cols() == 3);
  // Correct solution
  REQUIRE(error < 0.06);


  // std::fstream file("data.csv", std::ios::out);
  // file << "#time;signal;noisysignal;fit;\n";
  // for (int64_t i = 0; i < x.size(); i++) {
  //   file << std::format("{};{};{};{};\n", x(i), yClean(i), yNoisy(i),
  //   fit(i));
  // }
  // file.close();
}
