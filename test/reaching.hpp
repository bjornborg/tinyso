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

#pragma once

#include <eigen3/Eigen/Dense>

Eigen::ArrayXd gaussianMembershipFunction(
    Eigen::ArrayXd const x, Eigen::ArrayXd const &a_position)
{
  // 2023-09-04 10:37:14 bb | Individual data,
  // 0: amplitude
  // 1: timeshift (statistical mean)
  // 2: variance

  double const amp = a_position.x();
  double const mu = a_position.y();
  double const variance = a_position.z();
  Eigen::ArrayXd val{amp * Eigen::exp(-(x - mu).square() / (2 * variance))};

  return val;
}

Eigen::ArrayXd superpositionGaussians(
    Eigen::ArrayXd const x, Eigen::MatrixXd const &a_position)
{
  Eigen::ArrayXd y;
  //  = gaussianMembershipFunction(x, a_particle);
  y.resize(x.size());
  y.setZero();
  for (auto const pos : a_position.colwise()) {
    y += gaussianMembershipFunction(x, pos);
  }
  return y;
}


Eigen::ArrayXd gaussianMembershipHeavyTailFunction(
    Eigen::ArrayXd const x, Eigen::ArrayXd const &a_position)
{

  // 2024-11-22 10:05:31 bb | Assuming onset at 0, offset at 1

  // 2023-09-04 10:37:14 bb | Individual data,
  // 0: amplitude
  // 1: timeshift (statistical mean)
  // 2: variance

  double const amp = a_position.x();
  double const mu = a_position.y();
  double const variance = a_position.z();
  double const onset = a_position(3);
  double const offset = a_position(4);
  auto t = (x - onset) / (offset - onset);
  Eigen::ArrayXd val{
      (amp / t) * Eigen::exp(-(Eigen::log(t) - mu).square() / (2 * variance))};

  return val;
}

Eigen::ArrayXd superpositionGaussiansHeavyTail(
    Eigen::ArrayXd const x, Eigen::MatrixXd const &a_position)
{
  Eigen::ArrayXd y;
  //  = gaussianMembershipFunction(x, a_particle);
  y.resize(x.size());
  y.setZero();
  for (auto const pos : a_position.colwise()) {
    y += gaussianMembershipHeavyTailFunction(x, pos);
  }
  return y;
}
