## tinysoa: a minimalistic header-only C++ library for stochastic optimization algorithms, such as a genetic algorithm

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

tinysoa is a small and efficient library written in modern C++ to provide functionality for stochastic optimization.

tinysoa is available as single-file, header-only library. Insert [tinysoa.hpp](https://raw.githubusercontent.com/olbender/tinysoa/master/tinysoa.hpp) into your project, `#include "tinysoa.hpp"`, and compile your project with a modern C++ compiler (C++11 or newer).


## Table of Contents
* [Features](#features)
* [Dependencies](#dependencies)
* [Contributing](#contributing)
* [License](#license)


## Features
* Written in highly portable and high quality C++11
* **Available as header-only, single-file distribution. Insert [tinysoa.hpp](https://raw.githubusercontent.com/olbender/tinysoa/master/tinysoa.hpp) into your project, `#include "tinysoa.hpp"`, and compile your project with a modern C++ compiler (C++11 or newer)**
* Genetic algorithm with creep mutation and tournament selection. Each individual is a list of doubles.


## Dependencies
No dependencies! All you need is a C++11-compliant compiler as the project ships the following dependencies as part of the source distribution:

* [Unit Test Framework Catch2](https://github.com/catchorg/Catch2/releases/tag/v2.1.1) - [![License: Boost Software License v1.0](https://img.shields.io/badge/License-Boost%20v1-blue.svg)](http://www.boost.org/LICENSE_1_0.txt) - [Source](https://github.com/olbender/tinysoa/blob/master/test/catch.hpp)


## Installation
### Installation as single-file, header-only library
tinysoa is provided as header-only, single-file library as well. Insert [tinysoa.hpp](https://raw.githubusercontent.com/olbender/tinysoa/master/tinysoa.hpp) into your project, `#include "tinysoa.hpp"`, and compile your project with a modern C++ compiler (C++11 or newer)

## License
* This project is released under the terms of the MIT License - [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
