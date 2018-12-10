/*
 * MIT License
 *
 * Copyright (c) 2018 Ola Benderius
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
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

#include "tinysoa.hpp"

TEST_CASE("Genetic algorithm")
{
  double const min_x = 2.0;
  double const min_y = 5.0;

  auto evaluate = [&min_x, &min_y](tinysoa::Individual const &individual, uint32_t const) -> double
  {
    double x = 20.0 * individual[0] - 10.0;
    double y = 20.0 * individual[1] - 10.0;
    double z = 1 + sqrt((x - min_x) * (x - min_x) + (y - min_y) * (y - min_y));
    return 1.0 / z;
  };

  tinysoa::GeneticAlgorithm ga(evaluate, tinysoa::CrossoverMethod::Split, 2, 1, 
      30, 2, 0.5, 0.05, 0.5);

  for (uint32_t i = 0; i < 1000; i++) {
    ga.NextGeneration(4);
    auto ind = ga.GetBestIndividual();

    double x = 20.0 * ind[0] - 10.0;
    double y = 20.0 * ind[1] - 10.0;
    std::cout << i << ": Best f(" << x << ", " << y << ") = " << 1 / evaluate(ind, 0) << std::endl;
  }

  auto bestInd = ga.GetBestIndividual();
  double x = 20.0 * bestInd[0] - 10.0;
  double y = 20.0 * bestInd[1] - 10.0;
  REQUIRE(std::abs(x - min_x) < 0.1);
  REQUIRE(std::abs(y - min_y) < 0.1);
}
