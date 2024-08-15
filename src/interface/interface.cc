#include <iostream>

#include "console.h"

int main(int argc, char* argv[]) {
  std::string method{};
  if (argc == 2) method = std::string(argv[1]);
  Parallels::Console::GetInstance().ShowConsole(method);
  return 0;
}
