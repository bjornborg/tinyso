# MIT License
#
# Copyright (c) 2018 Ola Benderius
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.10)

project(tinyso)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} \
      -D_FORTIFY_SOURCE=2 -O3 -s")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} \
      -g -fsanitize=thread")
  # -fsanitize=address
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} \
      -D_XOPEN_SOURCE=700 -fstack-protector \
      -fomit-frame-pointer -pipe -pedantic -pedantic-errors -Werror -Weffc++ \
      -Wall -Wextra -Wshadow -Wdeprecated -Wdiv-by-zero -Wfloat-equal \
      -Wfloat-conversion -Wsign-compare -Wpointer-arith -Wuninitialized \
      -Wunreachable-code -Wunused -Wunused-function -Wunused-label \
      -Wunused-parameter -Wunused-but-set-parameter -Wunused-but-set-variable \
      -Wunused-value -Wunused-variable -Wunused-result \
      -Wmissing-field-initializers -Wmissing-format-attribute \
      -Wmissing-include-dirs -Wmissing-noreturn")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty)

# ##############################################################################
# gnuplot-iostream lib
find_package(Boost REQUIRED COMPONENTS iostreams system filesystem)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/gnuplot-iostream)
set(GnuPlotIostream_BuildTests
    OFF
    CACHE INTERNAL "build gnuplot iostream tests" FORCE)
set(GnuPlotIostream_BuildExamples
    OFF
    CACHE INTERNAL "build gnuplot iostream examples" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/gnuplot-iostream)
set(LIBRARIES ${LIBRARIES} gnuplot_iostream Boost::iostreams Boost::system
              Boost::filesystem)

set(LIBRARIES ${LIBRARIES} Threads::Threads)

enable_testing()
add_executable(
  ${PROJECT_NAME}-ga-runner ${CMAKE_CURRENT_SOURCE_DIR}/test/test-ga.cpp
  # $<TARGET_OBJECTS:tinyso-core>
)
target_link_libraries(${PROJECT_NAME}-ga-runner ${LIBRARIES})
add_test(NAME ${PROJECT_NAME}-ga-runner COMMAND ${PROJECT_NAME}-ga-runner)

add_executable(
  ${PROJECT_NAME}-pso-runner ${CMAKE_CURRENT_SOURCE_DIR}/test/test-pso.cpp
  # $<TARGET_OBJECTS:tinyso-core>
)
target_link_libraries(${PROJECT_NAME}-pso-runner ${LIBRARIES})
add_test(NAME ${PROJECT_NAME}-pso-runner COMMAND ${PROJECT_NAME}-pso-runner)
