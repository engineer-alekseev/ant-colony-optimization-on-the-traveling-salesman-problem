#ifndef PARALLELS_CONSOLE_H
#define PARALLELS_CONSOLE_H

#include <chrono>
#include <functional>
#include <iostream>
#include <map>
#include <regex>

#include "../aco/aco.h"
#include "../matrix/matrix.h"

namespace Parallels {

class Console {
  const std::map<std::string, std::function<void()>> options_{
    {"ant", std::bind(&Console::OptionAnt, this)}};

 public:
  Console(const Console&) = delete;
  ~Console() = default;
  
  void ShowConsole(std::string& method_name) const {
    if (method_name.empty()) {
      method_name = "ant";
    }
    while (method_name != "q") {
      DoAlgorithms(method_name);
      std::cout << "\n=============================================\n";
      std::cout << "Enter 'q' to exit or press Enter to rerun 'ant': ";
      std::getline(std::cin >> std::ws, method_name);
      if (method_name.empty()) {
        method_name = "ant";
      }
    }
    std::cout << "Exit\n";
  }

  static Console& GetInstance() {
    static Console instance_;
    return instance_;
  }

 private:
  Console() = default;

  void DoAlgorithms(const std::string& method_name) const {
    std::system("clear");
    auto choice = options_.find(ToLower(method_name));
    if (choice != options_.end()) {
      choice->second();
    } else {
      std::cout << "Method not found\n";
    }
  }

  std::string ToLower(std::string s) const {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    return s;
  }

  void OptionAnt() {
    ShowMethodName("Ant Method");
    std::string answer{"yes"};
    while (std::regex_match(ToLower(answer), std::regex("yes|y"))) {
      try {
        std::cout << "Enter matrix size\n";
        int size{};
        std::cin >> size;
        ClearInput();
        matrix_ = EnterMatrix(size, true);
        matrix_.PrintToFile("AntMatrixGenerated.txt");
        std::cout << "\n=============================================\n";
        int iterations = GetIterations();
        time_ = CountTime([this, iterations](int) { CallAntUsual(iterations); },
                          iterations);
        ShowTime("Ant usual: ");
        time_ =
            CountTime([this, iterations](int) { CallAntParallel(iterations); },
                      iterations);
        ShowTime("Ant parallel: ");
      } catch (std::invalid_argument& e) {
        std::cout << "=============================================\n";
        std::cout << e.what() << std::endl;
        std::cout << "=============================================\n";
      }
      std::cout << "\n=============================================\n";
      std::cout << "Do you want to continue? Enter yes or no\n";
      std::cin >> answer;
    }
  }

  int GetIterations() const {
    std::cout << "Enter amount of iterations\n";
    int iterations{};
    std::cin >> iterations;
    if (iterations < 0)
      throw std::invalid_argument("Incorrect number of iterations");
    return iterations;
  }

  Matrix EnterMatrix(int size, bool graph = false) {
    if (size < 1) throw std::invalid_argument("Incorrect matrix size");
    Matrix result(size);
    if (std::regex_match(ToLower(CheckAnswer()), std::regex("yes|y"))) {
      if (graph) {
        result.FillRandomMatrixGraph();
      } else {
        result.FillRandomMatrix();
      }
    } else {
      std::cout << "Enter matrix value\n";
      double value{};
      for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
          if (std::cin >> value) {
            result.GetMatrix()[i][j] = value;
          } else {
            ClearInput();
            throw std::invalid_argument("Incorrect matrix value");
          }
        }
      }
    }
    ClearInput();
    return result;
  }

  std::string CheckAnswer() const {
    std::cout << "Do you want to create a random matrix? Enter yes or no\n";
    std::string answer{};
    std::cin >> answer;
    return answer;
  }

  void CallAntUsual(int iterations) {
    try {
      std::cout << "\t\tUsual Ant Result\n";
      auto result = aco_ex_.RunSequential(matrix_, iterations);
      std::cout << "Best path: ";
      PrintResultVector(result.first);
      std::cout << "Distance length: " << result.second << std::endl;
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
  }

  void CallAntParallel(int iterations) {
    try {
      std::cout << "\t\tParallel Ant Result\n";
      auto result = aco_ex_.RunParallel(matrix_, iterations);
      std::cout << "Best path: ";
      PrintResultVector(result.first);
      std::cout << "Distance length: " << result.second << std::endl;
    } catch (std::invalid_argument& e) {
      std::cout << e.what() << std::endl;
    }
  }

  double CountTime(std::function<void(int)> f, int iterations) {
    const auto start{std::chrono::system_clock::now()};
    f(iterations);
    const auto finish{std::chrono::system_clock::now()};
    return std::chrono::duration_cast<std::chrono::milliseconds>(finish - start)
        .count();
  }

  void ClearInput() const {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  void ShowMethodName(const std::string& method_name) const {
    std::cout << "=============================================\n";
    std::cout << "\t\t" << method_name << "\n";
    std::cout << "=============================================\n";
  }

  void ShowTime(const std::string& message) const {
    std::cout << "\n=============================================\n";
    std::cout << message << time_ << " ms \n";
    std::cout << "\n=============================================\n";
  }

  void PrintResultVector(const std::vector<int>& v) const {
    for (auto value : v) {
      std::cout << value << " ";
    }
    std::cout << std::endl;
  }

  AcoExecutor aco_ex_;
  Matrix matrix_;
  double time_{};
};

}  // namespace Parallels

#endif  // PARALLELS_CONSOLE_H
