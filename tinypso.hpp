#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <eigen3/Eigen/Dense>
#include <future>
#include <vector>

class Particle; // 2023-09-25 13:04:01 bb |  Probably ok?
typedef std::vector<Particle> Swarm;

struct PsoSettings
{
  std::function<double(Particle const &)> objectiveFunction;
  std::function<bool(Particle const &)> plotFunction;
  std::function<void(Particle const &, std::shared_ptr<PsoSettings> const &)> printResultFunction;
  uint32_t particleRow{0}; // dimensionality
  uint32_t particleCol{0}; // number of bumps, may be increased after a while
  uint32_t swarmSize{0};
  Eigen::VectorXd posMin;
  Eigen::VectorXd posMax;
  Eigen::VectorXd initialGuess;
  double alpha{1.0};
  double dt{1.0};
  double speedMax{0.0};
  double maxInertia{1.4}; // exploration vs exploitation, modified throughout the optimization
  double minInertia{0.3};
  double inertia{0.0};
  double inertiaDecay{0.995}; // transforms exploration -> exploitation
  double cognitiveFactor{2.0};
  double socialFactor{2.0};
  double m_probCrazy{0.001}; // actually makes it worse?

  uint32_t numSteps{1000};
  uint32_t threads{16};
  uint32_t convergeIterations{100};
  double minumumFitness{1.0};
  uint32_t maxCol{5};
};

class Particle
{
public:
  Particle(std::shared_ptr<PsoSettings> const &a_settings);
  ~Particle();
  Particle(Particle const &) = default;
  Particle &operator=(Particle const &) = default;

  bool step();
  bool evaluateCurrentFitness();
  Eigen::MatrixXd getPosition() const;
  bool setPosition(Eigen::MatrixXd const &);
  bool setGuess(Eigen::MatrixXd const &);
  Eigen::MatrixXd getVelocity() const;
  bool setVelocity(Eigen::MatrixXd const &);
  bool updateVelocity(Particle const &a_particle);
  bool resize();

  std::tuple<double, Eigen::MatrixXd> getBestFitness() const;
  std::tuple<double, Eigen::MatrixXd> getCurrentFitness() const;

private:
  Eigen::MatrixXd initPosition();
  Eigen::MatrixXd initVelocity();

  std::string m_name;
  std::shared_ptr<PsoSettings> m_settings;
  Eigen::MatrixXd m_position;
  Eigen::MatrixXd m_velocity;
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> m_permutation;
  double m_currentFitness;
  double m_bestFitness;
  Eigen::MatrixXd m_bestPosition;
};

inline Particle::Particle(std::shared_ptr<PsoSettings> const &a_settings)
    : m_name{"[Pso particle] "},
      m_settings{a_settings},
      m_position{initPosition()},
      m_velocity{initVelocity()},
      m_permutation{Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>(m_settings->particleCol)},
      m_currentFitness{0.0},
      m_bestFitness{0.0},
      m_bestPosition{m_position}
{
  if (!m_settings)
  {
    throw std::runtime_error(m_name + " Pso settings is a null pointer");
  }
  // resize();
  m_permutation.setIdentity();
}

inline Particle::~Particle() {}

inline bool Particle::step()
{
  return setPosition(m_position + m_velocity * m_settings->dt);
}

// Return position in uniform random in range [posMin, posMax]
inline bool Particle::evaluateCurrentFitness()
{
  m_currentFitness = m_settings->objectiveFunction(*this);
  if (m_currentFitness > m_bestFitness)
  {
    m_bestFitness = m_currentFitness;
    m_bestPosition = m_position;
  }
  return true;
}

inline Eigen::MatrixXd Particle::getPosition() const
{
  return m_position;
}

inline bool Particle::setPosition(Eigen::MatrixXd const &a_position)
{
  // 2023-09-21 16:42:45 bb |  Clamp the position within defined box set
  m_position = a_position.cwiseMin(m_settings->posMax.replicate(1, m_settings->particleCol)).cwiseMax(m_settings->posMin.replicate(1, m_settings->particleCol));

  // 2023-09-21 16:42:47 bb | Order the time dimension to truncate the search space a bit
  std::sort(m_permutation.indices().begin(), m_permutation.indices().end(),
            [data = m_position.row(1)](int32_t const &a, int32_t const &b)
            { return data[a] < data[b]; });

  m_position *= m_permutation;
  m_velocity *= m_permutation;
  m_permutation.setIdentity();
  return true;
}

inline bool Particle::setGuess(Eigen::MatrixXd const &a_guess)
{

  Eigen::MatrixXd guess = initPosition();
  if (a_guess.cols() > guess.cols())
    return false;

  guess.leftCols(a_guess.cols()) = a_guess;
  return setPosition(guess);
}

inline Eigen::MatrixXd Particle::getVelocity() const
{
  return m_velocity;
}

inline bool Particle::setVelocity(Eigen::MatrixXd const &a_velocity)
{
  // Clamp the velocity within the speedMax limits
  // Using the euclidian norm might be computational taxing. If it's too maxing, might consider to switch over box clamping instead.
  m_velocity = (a_velocity.norm() > m_settings->speedMax)
                   ? a_velocity.normalized() * m_settings->speedMax
                   : a_velocity;
  return true;
}

inline bool Particle::updateVelocity(Particle const &a_particleBest)
{
  // Eigen::MatrixXd newVelocity;
  Eigen::Vector3d r{Eigen::Vector3d::Random().cwiseAbs()}; // stochastic varibles [0,1)
  // crazy check
  // if (r.z() < m_settings->m_probCrazy)
  // {
  //   newVelocity = initVelocity();
  // }
  // else
  // {
  //   // inertia
  //   newVelocity = m_settings->inertia * m_velocity
  //                 // cognitive contribution
  //                 + m_settings->cognitiveFactor * r.x() * (m_bestPosition - m_position) / m_settings->dt
  //                 // social contribution
  //                 + m_settings->socialFactor * r.y() * (std::get<1>(a_particleBest.getBestFitness()) - m_position) / m_settings->dt;
  // }

  setVelocity(r.z() < m_settings->m_probCrazy
                  ? initVelocity() // crazy
                  : m_settings->inertia * m_velocity
                        // cognitive contribution
                        + m_settings->cognitiveFactor * r.x() * (m_bestPosition - m_position) / m_settings->dt
                        // social contribution
                        + m_settings->socialFactor * r.y() * (std::get<1>(a_particleBest.getBestFitness()) - m_position) / m_settings->dt);
  return true;
}

inline bool Particle::resize()
{
  m_permutation.resize(m_settings->particleCol);
  m_permutation.setIdentity();

  auto newVelocity = initVelocity();
  if (m_velocity.cols() <= newVelocity.cols())
    newVelocity.leftCols(m_velocity.cols()) = m_velocity;
  else
    newVelocity = m_velocity.leftCols(newVelocity.cols());

  setVelocity(newVelocity);

  auto newPosition = initPosition();
  if (m_position.cols() <= newPosition.cols())
    newPosition.leftCols(m_position.cols()) = m_position;
  else
    newPosition = m_position.leftCols(newPosition.cols());

  setPosition(newPosition);

  if (m_bestPosition.cols() <= newPosition.cols())
    newPosition.leftCols(m_bestPosition.cols()) = m_bestPosition;
  else
    newPosition = m_bestPosition.leftCols(newPosition.cols());

  std::sort(m_permutation.indices().begin(), m_permutation.indices().end(),
            [data = newPosition.row(1)](int32_t const &a, int32_t const &b)
            { return data[a] < data[b]; });

  m_bestPosition = newPosition * m_permutation;
  m_permutation.setIdentity();

  m_bestFitness = m_settings->objectiveFunction(*this);

  return true;
}

inline std::tuple<double, Eigen::MatrixXd> Particle::getBestFitness() const
{
  return std::make_tuple(m_bestFitness, m_bestPosition);
}

inline std::tuple<double, Eigen::MatrixXd> Particle::getCurrentFitness() const
{
  return std::make_tuple(m_currentFitness, m_position);
}

inline Eigen::MatrixXd Particle::initPosition()
{
  Eigen::MatrixXd position =
      m_settings->posMin.replicate(1, m_settings->particleCol).array() +
      // random in elementwise interval [-1,1], abs -> [0,1]
      Eigen::MatrixXd::Random(m_settings->particleRow, m_settings->particleCol).cwiseAbs().array() * (m_settings->posMax - m_settings->posMin).replicate(1, m_settings->particleCol).array();

  return position;
}

inline Eigen::MatrixXd Particle::initVelocity()
{
  Eigen::MatrixXd velocity =
      m_settings->alpha / m_settings->dt +
      ((-(m_settings->posMax - m_settings->posMin).replicate(1, m_settings->particleCol).array() / 2) +
       Eigen::MatrixXd::Random(m_settings->particleRow, m_settings->particleCol).cwiseAbs().array() * (m_settings->posMax - m_settings->posMin).replicate(1, m_settings->particleCol).array());
  return velocity;
}

class ParticleSwarmOptimization
{
public:
  ParticleSwarmOptimization(
      std::shared_ptr<PsoSettings> const a_psoSettings);
  ~ParticleSwarmOptimization();

  int32_t getStepCounter() const;
  bool step();
  Particle runOptimization();
  Swarm getSwarm() const;

private:
  ParticleSwarmOptimization(ParticleSwarmOptimization const &);
  ParticleSwarmOptimization &operator=(ParticleSwarmOptimization const &);

  double getBestFitness() const;
  std::vector<double> getFitnesses() const;
  Particle getBestParticle() const;
  bool evaluateSwarm(Swarm &, std::vector<double> &, size_t &a_bestParticleIndex, uint32_t &a_convergeCounter, uint32_t const) const;
  bool stepSwarm(Swarm &, std::vector<double> &, size_t &a_bestParticleIndex, uint32_t const) const;
  Particle *findBestParticle(std::vector<Particle> &a_swarm, std::vector<double> const &a_fitnesses) const;
  Swarm generatePopulation(std::shared_ptr<PsoSettings> const &a_psoSettings);

  double updateInertia(double const, double const) const;
  Particle trimTruncate(Particle const &a_particle) const;

  std::string const m_name;
  std::shared_ptr<PsoSettings> const m_psoSettings;

  std::vector<double> m_fitnesses;
  double m_currentOverallFitness;
  uint32_t m_stepCounter;
  uint32_t m_convergeCounter;
  Swarm m_swarm;
  // Particle *m_bestParticle;
  size_t m_bestParticleIndex;
};

inline ParticleSwarmOptimization::ParticleSwarmOptimization(
    std::shared_ptr<PsoSettings> const a_psoSettings)
    : m_name{"[Pso] "},
      m_psoSettings{a_psoSettings},
      m_fitnesses{std::vector<double>(m_psoSettings->swarmSize)},
      m_currentOverallFitness{0.0},
      m_stepCounter{0},
      m_convergeCounter{0},
      m_swarm{generatePopulation(m_psoSettings)},
      // m_bestParticle{&m_swarm.front()}
      m_bestParticleIndex{0}
{
  if (!m_psoSettings)
  {
    throw std::runtime_error(m_name + " settings pointer is a nullpointer");
  }
}
inline ParticleSwarmOptimization::~ParticleSwarmOptimization() {}
inline double ParticleSwarmOptimization::getBestFitness() const
{
  return std::get<double>(m_swarm[m_bestParticleIndex].getBestFitness());
}
inline Particle ParticleSwarmOptimization::getBestParticle() const
{
  return m_swarm[m_bestParticleIndex];
}
inline int32_t ParticleSwarmOptimization::getStepCounter() const
{
  return m_stepCounter;
}
inline bool ParticleSwarmOptimization::step()
{
  m_stepCounter++;
  m_convergeCounter++;
  // Update inertia
  m_psoSettings->inertia = updateInertia(m_psoSettings->inertia, m_psoSettings->inertiaDecay);

  // check and resize
  if (m_swarm.front().getPosition().cols() != m_psoSettings->particleCol)
  {
    // std::clog << "Resizing to " << std::to_string(m_psoSettings->particleCol) << " cols." << std::endl;
    for (auto &particle : m_swarm)
    {
      particle.resize();
    }
    m_bestParticleIndex = 0;
  }
  // 2023-09-21 16:17:09 bb | Evaluate the swarm
  evaluateSwarm(m_swarm, m_fitnesses, m_bestParticleIndex, m_convergeCounter, m_psoSettings->threads);
  stepSwarm(m_swarm, m_fitnesses, m_bestParticleIndex, m_psoSettings->threads);

  return true;
}
inline Particle ParticleSwarmOptimization::runOptimization()
{
  bool converged{false};
  Particle bestParticle = getBestParticle();

  while (!converged)
  {
    // Reset counters
    m_psoSettings->inertia = m_psoSettings->maxInertia;
    m_convergeCounter = 0;
    m_stepCounter = 0;
    for (uint32_t i = 0; i < m_psoSettings->numSteps && m_convergeCounter < m_psoSettings->convergeIterations; i++)
    {
      step();
      if (m_psoSettings->plotFunction)
        m_psoSettings->plotFunction(getBestParticle());
    }
    Particle candidate = getBestParticle();
    // auto const fitness = candidate.getBestFitness();

    if ((std::get<0>(candidate.getBestFitness()) < std::get<0>(bestParticle.getBestFitness()) &&
         std::get<0>(bestParticle.getBestFitness()) > m_psoSettings->minumumFitness) ||
        m_psoSettings->particleCol > m_psoSettings->maxCol)
    {
      // 2023-10-04 13:57:18 bb | Found a good solution with minimal improvements and a minimal requirement of some base score
      converged = true;
    }
    else
    {
      // Keep expanding and look for new solution
      bestParticle = candidate;
      m_psoSettings->particleCol++;
    }
  }
  // truncate and remove trivial solutions in the particle
  bestParticle = trimTruncate(bestParticle);

  if (m_psoSettings->printResultFunction)
    m_psoSettings->printResultFunction(bestParticle, m_psoSettings);

  return bestParticle;
}
inline Swarm ParticleSwarmOptimization::getSwarm() const
{
  return m_swarm;
}
inline std::vector<double> ParticleSwarmOptimization::getFitnesses() const
{
  return m_fitnesses;
}
inline bool ParticleSwarmOptimization::evaluateSwarm(Swarm &a_swarm, std::vector<double> &a_fitnesses, size_t &a_bestParticleIndex, uint32_t &a_convergeCounter, uint32_t const a_threads) const
{

  std::mutex fitnessWriteMutex;
  std::mutex evaluatedWriteMutex;
  std::mutex bestParticleMutex;
  uint32_t evaluated = 0;
  auto worker{
      [this, &a_swarm, &a_fitnesses, &evaluated, &fitnessWriteMutex, &evaluatedWriteMutex, &a_bestParticleIndex, &bestParticleMutex, &a_convergeCounter]() mutable
      {
        while (evaluated < a_fitnesses.size())
        {
          uint32_t index;
          {
            std::lock_guard<std::mutex> lock(evaluatedWriteMutex);
            if (evaluated >= a_fitnesses.size())
            {
              break;
            }
            index = evaluated;
            evaluated++;
          }
          a_swarm[index].evaluateCurrentFitness();
          // You want the  all time best fitness.
          // double fitness = std::get<double>(a_swarm[index].getCurrentFitness());
          double const fitness = std::get<double>(a_swarm[index].getBestFitness());
          {
            std::lock_guard<std::mutex> lock(fitnessWriteMutex);
            a_fitnesses[index] = fitness;
          }
          {
            std::lock_guard<std::mutex> lock(bestParticleMutex);
            if (fitness > std::get<double>(a_swarm[a_bestParticleIndex].getBestFitness()))
            {
              a_bestParticleIndex = index;
              a_convergeCounter = 0;
            }
          }
        }
      }};
  std::vector<std::thread> threads;
  for (uint32_t i{0}; i < a_threads; i++)
  {
    threads.emplace_back(worker);
  }
  for (auto &t : threads)
  {
    t.join();
  }

  return true;
}
inline bool ParticleSwarmOptimization::stepSwarm(Swarm &a_swarm, std::vector<double> &, size_t &a_bestParticleIndex, uint32_t const) const
{
  std::mutex evaluatedWriteMutex;
  uint32_t evaluated{0};
  auto worker{
      [&a_swarm,
       &a_bestParticleIndex,
       &evaluated,
       &evaluatedWriteMutex]() -> void
      {
        while (evaluated < a_swarm.size())
        {
          uint32_t index;
          // Get the next index
          {
            std::lock_guard<std::mutex> lock(evaluatedWriteMutex);
            if (evaluated >= a_swarm.size())
            {
              break;
            }
            index = evaluated;
            evaluated++;
          }
          // Update the velocities
          a_swarm[index].updateVelocity(a_swarm[a_bestParticleIndex]);
          // Update the positons
          a_swarm[index].step();
        }
      }};

  std::vector<std::thread> threads;
  // Start the jobs
  for (uint32_t i{0}; i < m_psoSettings->threads; i++)
  {
    threads.push_back(std::thread(worker));
  }
  // Merge them
  for (auto &t : threads)
  {
    t.join();
  }
  return true;
}
inline Particle *ParticleSwarmOptimization::findBestParticle(std::vector<Particle> &a_swarm, std::vector<double> const &a_fitnesses) const
{
  auto const it = std::max_element(a_fitnesses.cbegin(), a_fitnesses.cend());
  size_t index = std::distance(a_fitnesses.begin(), it);

  return &a_swarm[index];
}
inline Swarm ParticleSwarmOptimization::generatePopulation(
    std::shared_ptr<PsoSettings> const &a_psoSettings)
{
  Swarm swarm;
  // swarm.reserve(a_psoSettings->swarmSize);
  for (size_t n{0}; n < a_psoSettings->swarmSize; n++)
  {
    swarm.emplace_back(m_psoSettings);
    // swarm.back().setPosition(guess);
  }
  // swarm.front().setPosition(a_psoSettings->initialGuess);
  return swarm;
}

inline double ParticleSwarmOptimization::updateInertia(double const a_inertia, double const a_beta) const
{
  return std::clamp(a_inertia * m_psoSettings->inertiaDecay, m_psoSettings->minInertia, m_psoSettings->maxInertia);
}

inline Particle ParticleSwarmOptimization::trimTruncate(Particle const &a_particle) const
{
  // To do if trimming and truncation is needed
  return a_particle;
}
