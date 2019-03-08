#define USE_NAMESPACE

#include <lat-io.h>

#include <string>
#include <vector>
#include <cassert>

int main(const int argc, const char* argv[])
{
  using namespace latio;
  if (argc != 2) {
    printf("usage: program-name filename\n");
  } else {
    const std::string fn(argv[1]);
    LatData ld;
    ld.load(fn);
    print(ld);
  }
  return 0;
};
