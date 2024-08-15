#ifndef PARALLELS_ACO_H
#define PARALLELS_ACO_H

#include <algorithm>
#include <iostream>
#include <mutex>
#include <random>
#include <set>
#include <stack>
#include <thread>
#include <vector>

#include "../matrix/matrix.h"

namespace Parallels {
class Ant {
 public:
  Ant(int position, const Matrix &distances, const Matrix &pheromones,
      double alpha = 1.0, double beta = 1.0)
      : started_from_(position),
        position_(position),
        distances_(distances),
        pheromones_(pheromones),
        alpha_(alpha),
        beta_(beta) {
    for (int idx = 0; idx < distances_.GetRows(); ++idx) {
      to_visit_.insert(idx);
    }
    to_visit_.erase(position);
    route_.push_back(position);
  }

  std::vector<std::pair<int, double>> GetAvailableDirections() {
    std::vector<std::pair<int, double>> available_directions;

    VVDouble distances = distances_.GetMatrix();
    for (int idx = 0; idx < distances_.GetRows(); ++idx) {
      double direction = distances[position_][idx];
      if (direction != 0.) {
        available_directions.push_back({idx, direction});
      }
    }
    return available_directions;
  }

  void Travel(int next_vertex) {
    position_ = next_vertex;
    to_visit_.erase(next_vertex);
    route_.push_back(next_vertex);
  }

  bool IsJourneyOver() { return to_visit_.empty(); }

  std::pair<int, double> FirstDirection() {
    std::vector<std::pair<int, double>> available_directions =
        GetAvailableDirections();
    return available_directions[rand() % available_directions.size()];
  }

  double GetDesire(double distance, double pheromone) {
    return pow(1.0 / distance, alpha_) * pow(pheromone, beta_);
  }

  double DesireToProbability(double desire, double sum_of_desires) {
    return desire / sum_of_desires;
  }

  double PheromoneForDirection(int direction) {
    return pheromones_.GetMatrix()[position_][direction];
  }

  std::vector<std::pair<int, double>> FilterNotVisited(
      const std::vector<std::pair<int, double>> &available_directions) {
    std::vector<std::pair<int, double>> filtered_directions;
    for (const auto &direction : available_directions) {
      if (to_visit_.find(direction.first) != to_visit_.end()) {
        filtered_directions.push_back(direction);
      }
    }
    return filtered_directions;
  }

  std::vector<std::pair<int, std::pair<double, double>>> PrepareRoullete(
      const std::vector<std::pair<int, double>> &directions) {
    std::vector<std::pair<int, std::pair<double, double>>> roulette;
    std::vector<std::pair<int, double>> desires_for_directions;

    for (const auto &direction : directions) {
      double pheromone = PheromoneForDirection(direction.first);
      double desire = GetDesire(direction.second, pheromone);
      desires_for_directions.push_back({direction.first, desire});
    }

    double sum_of_desires = 0.0;
    for (const auto &d : desires_for_directions) {
      sum_of_desires += d.second;
    }

    double mark = 0.0;
    for (const auto &desire : desires_for_directions) {
      double probability = DesireToProbability(desire.second, sum_of_desires);
      roulette.push_back({desire.first, {mark, mark + probability}});
      mark += probability;
    }

    return roulette;
  }

  int MakeDecision(
      const std::vector<std::pair<int, std::pair<double, double>>> &roulette,
      double dice_roll) {
    for (const auto &entry : roulette) {
      if (dice_roll >= entry.second.first && dice_roll < entry.second.second) {
        return entry.first;
      }
    }
    return -1;
  }

  int NextDirection(bool not_visited = true) {
    std::vector<std::pair<int, double>> available_directions =
        GetAvailableDirections();
    std::vector<std::pair<int, double>> choose_from = available_directions;

    if (not_visited) {
      choose_from = FilterNotVisited(available_directions);
    }

    if (!choose_from.empty()) {
      std::vector<std::pair<int, std::pair<double, double>>> roulette =
          PrepareRoullete(choose_from);
      double dice_roll = static_cast<double>(rand()) / RAND_MAX;
      return MakeDecision(roulette, dice_roll);
    } else {
      return -1;
    }
  }

  void OneWayTour() {
    while (!IsJourneyOver()) {
      int next_vertex = NextDirection(true);
      if (next_vertex != -1) {
        Travel(next_vertex);
      } else {
        next_vertex = NextDirection(false);
        Travel(next_vertex);
      }
    }
  }

  std::vector<int> FullTour() {
    OneWayTour();
    to_visit_.insert(started_from_);
    OneWayTour();
    return route_;
  }

  double GetRouteLength() {
    VVDouble distances = distances_.GetMatrix();
    if (route_.size() > 1) {
      double route_length = 0.0;
      for (int idx = 0; idx < static_cast<int>(route_.size() - 1); ++idx) {
        route_length += distances[route_[idx]][route_[idx + 1]];
      }
      return route_length;
    } else {
      return 0.0;
    }
  }

  std::vector<std::vector<double>> GetPheromoneAdditive() {
    int n_vertices = distances_.GetCols();
    double Q = static_cast<double>(n_vertices);
    double route_length = GetRouteLength();
    double additive = Q / route_length;

    VVDouble additives(n_vertices, std::vector<double>(n_vertices, 0.0));
    const std::vector<int> &route = route_;

    for (int idx = 0; idx < static_cast<int>(route_.size() - 1); ++idx) {
      additives[route[idx]][route[idx + 1]] += additive;
      additives[route[idx + 1]][route[idx]] += additive;
    }

    return additives;
  }

 private:
  int started_from_;
  int position_;
  Matrix distances_;
  Matrix pheromones_;
  std::set<int> to_visit_;
  std::vector<int> route_;
  double alpha_;
  double beta_;
};

class Aco {
 public:
  struct AntTourResult {
    std::vector<int> route;
    double route_length;
    std::vector<std::vector<double>> additive;
  };

  explicit Aco(const Matrix &input, double initial_pheromone = 0.2,
               double alpha = 1, double beta = 1,
               double pheromone_leftover = 0.6, int n_epochs_max = 1000,
               int n_epochs_stable = 50)
      : distances_(input),
        pheromones_(input.GetRows(), input.GetCols()),
        alpha_(alpha),
        beta_(beta),
        p_(pheromone_leftover),
        best_route_length_(0.0),
        n_epochs_max_(n_epochs_max),
        n_epochs_stable_(n_epochs_stable) {
    InitPheromonMatrix(distances_, initial_pheromone);
  };
  std::vector<int> GetRoute() { return best_route_; }

  double GetRouteLength() { return best_route_length_; }

  void Execute(bool in_parallel = false) {
    if (IsDisconnected(distances_)) {
      throw std::invalid_argument("Graph is disconnected");
    }
    int n_epochs_max = n_epochs_max_;
    int n_epochs_stable = n_epochs_stable_;
    int counter = 1;
    while (n_epochs_max > 0 && n_epochs_stable > 0) {
      bool is_it_a_new_route = ExecuteEpoch(in_parallel);
      if (!is_it_a_new_route) {
        --n_epochs_stable;
      } else {
        n_epochs_stable = n_epochs_stable_;
      }

      ++counter;
      --n_epochs_max;
    }
  }

 private:
  void DFS(Matrix &graph, int start, std::vector<bool> &visited) {
    int num_vertices = graph.GetRows();
    std::stack<int> stack;
    stack.push(start - 1);
    visited[start - 1] = true;
    while (!stack.empty()) {
      int current_vertex = stack.top();
      stack.pop();
      if (!visited[current_vertex]) {
        visited[current_vertex] = true;
      }
      for (int neighbor = num_vertices - 1; neighbor >= 0; --neighbor) {
        if (graph.GetMatrix()[current_vertex][neighbor] != 0 &&
            !visited[neighbor]) {
          stack.push(neighbor);
        }
      }
    }
  }

  bool IsDisconnected(Matrix &graph) {
    int numVertices = graph.GetRows();

    std::vector<bool> visited(numVertices, false);

    DFS(graph, 1, visited);

    for (bool v : visited) {
      if (!v) {
        return true;
      }
    }

    return false;
  }

  AntTourResult SendAntToTour(int vertexToStart) {
    Ant ant = CreateAnt(vertexToStart);
    std::vector<int> route = ant.FullTour();
    double route_length = ant.GetRouteLength();
    std::vector<std::vector<double>> additive = ant.GetPheromoneAdditive();

    AntTourResult result;
    result.route = route;
    result.route_length = route_length;
    result.additive = additive;

    return result;
  }

  std::vector<AntTourResult> SendAntsToTour() {
    int nAnts = distances_.GetRows();
    std::vector<AntTourResult> results;

    for (int idx = 0; idx < nAnts; ++idx) {
      AntTourResult result = SendAntToTour(idx);
      results.push_back(result);
    }

    return results;
  }

  std::vector<AntTourResult> SendAntsToTourParallel() {
    int nAnts = distances_.GetRows();
    std::vector<AntTourResult> results(nAnts);
    std::mutex mutex;

    unsigned int thread_number = std::min(
        (std::thread::hardware_concurrency() - 1), (unsigned int)nAnts);

    std::vector<std::thread> threads(thread_number);
    int idx = 0;
    while (idx < nAnts) {
      size_t t = 0;
      for (; t < thread_number && idx < nAnts; ++idx, ++t) {
        threads.at(t) = std::move(std::thread([&, idx]() {
          AntTourResult result = SendAntToTour(idx);
          std::lock_guard<std::mutex> lock(mutex);
          results[idx] = result;
        }));
      }
      for (auto &th : threads) {
        if (th.joinable()) th.join();
      }
      t = 0;
    }

    return results;
  }

  AntTourResult ReduceTourResults(const std::vector<AntTourResult> &results) {
    std::vector<std::vector<double>> additive = results[0].additive;
    std::vector<int> bestRoute = results[0].route;
    double bestRouteLength = results[0].route_length;

    for (size_t i = 1; i < results.size(); ++i) {
      if (results[i].route_length < bestRouteLength) {
        bestRoute = results[i].route;
        bestRouteLength = results[i].route_length;
      }

      for (size_t row = 0; row < additive.size(); ++row) {
        for (size_t col = 0; col < additive[row].size(); ++col) {
          additive[row][col] += results[i].additive[row][col];
        }
      }
    }

    AntTourResult result;
    result.route = bestRoute;
    result.route_length = bestRouteLength;
    result.additive = additive;

    return result;
  }

  bool UpdateState(const AntTourResult &result) {
    bool isItANewRoute = false;

    auto &pheromonesMatrix = pheromones_.GetMatrix();
    for (size_t row = 0; row < pheromonesMatrix.size(); ++row) {
      for (size_t col = 0; col < pheromonesMatrix[row].size(); ++col) {
        pheromonesMatrix[row][col] *= p_;
        pheromonesMatrix[row][col] += result.additive[row][col];
      }
    }

    if (best_route_length_ == 0.0 || result.route_length < best_route_length_) {
      best_route_ = result.route;
      best_route_length_ = result.route_length;
      isItANewRoute = true;
    }

    return isItANewRoute;
  }

  bool ExecuteEpoch(bool in_parallel = false) {
    std::vector<AntTourResult> results;
    if (in_parallel) {
      results = SendAntsToTourParallel();
    } else {
      results = SendAntsToTour();
    }
    AntTourResult reducedResult = ReduceTourResults(results);
    bool isItANewRoute = UpdateState(reducedResult);

    return isItANewRoute;
  }

  void InitPheromonMatrix(const Matrix &distances, double initial_pheromone) {
    VVDouble distances_matrix = distances_.GetMatrix();
    for (int i = 0; i < distances.GetRows(); ++i) {
      for (int j = 0; j < distances.GetCols(); ++j) {
        if (distances_matrix[i][j] != 0) {
          pheromones_.GetMatrix()[i][j] = initial_pheromone;
        }
      }
    }
  }

  Ant CreateAnt(int startPosition) {
    Matrix localDistances{distances_};
    Matrix localPheromones{pheromones_};

    return Ant(startPosition, localDistances, localPheromones, alpha_, beta_);
  }

  Matrix distances_;
  Matrix pheromones_;
  double alpha_;
  double beta_;
  double p_;
  std::vector<int> best_route_;
  double best_route_length_;
  double n_epochs_max_;
  double n_epochs_stable_;
};

class AcoExecutor {
 public:
  AcoExecutor() {}
  std::pair<std::vector<int>, double> RunSequential(const Matrix &distances,
                                                    int execute_iterations) {
    return Run(distances, execute_iterations);
  }
  std::pair<std::vector<int>, double> RunParallel(const Matrix &distances,
                                                  int execute_iterations) {
    return Run(distances, execute_iterations, true);
  }

 private:
  std::pair<std::vector<int>, double> Run(const Matrix &distances,
                                          int execute_iterations,
                                          bool in_parallel = false) {
    std::vector<int> best_route;
    double best_route_length = 0.0;
    for (int i = 0; i < execute_iterations; ++i) {
      Aco aco_instance(distances);
      aco_instance.Execute(in_parallel);
      if (i == execute_iterations - 1) {
        best_route = aco_instance.GetRoute();
        best_route_length = aco_instance.GetRouteLength();
      }
    }
    return std::make_pair(best_route, best_route_length);
  }
};
};  // namespace Parallels

#endif  // PARALLELS_ACO_H
