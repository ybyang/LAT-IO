#define USE_NAMESPACE

#include <lat-io.h>

#include <string>
#include <vector>
#include <cassert>

namespace latio {

inline std::string idx_name(const LatDim& dim, const long idx)
{
  using namespace qshow;
  if (idx < dim.indices.size()) {
    return dim.indices[idx];
  } else {
    return show(-idx-1);
  }
}

inline long lat_data_offset(const LatInfo& info, const std::vector<long>& idx)
{
  assert(idx.size() <= info.size());
  long ret = 0;
  for (int i = 0; i < idx.size(); ++i) {
    const long k = idx[i];
    ret = ret * info[i].size + k;
  }
  return ret;
}

inline void print(const LatData& ld)
{
  using namespace qshow;
  const LatInfo& info = ld.info;
  printf("%s", show(info).c_str());
  std::vector<long> idx(info.size(), 0);
  for (long k = 0; k < lat_data_size(info); ++k) {
    for (int a = 0; a < info.size(); ++a) {
      printf("%s[%8s] ", info[a].name.c_str(), idx_name(info[a], idx[a]).c_str());
    }
    printf("%24.17E\n", ld.res[lat_data_offset(info, idx)]);
    idx[info.size() - 1] += 1;
    for (int a = info.size() - 1; a > 0; --a) {
      if (idx[a] == info[a].size) {
        idx[a] = 0;
        idx[a-1] += 1;
      }
    }
  }
}

}

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
