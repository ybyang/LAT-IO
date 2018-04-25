#include <lat-io.h>

inline LatData mk_lat_data()
{
  LatData ld;
  LatInfo& info = ld.info;
  LatDim d1;
  d1.name = "d1";
  d1.size = 10;
  d1.indices.push_back("3");
  d1.indices.push_back("7");
  d1.indices.push_back("2");
  LatDim d2;
  d2.name = "d2";
  d2.size = 3;
  d2.indices.push_back("a");
  d2.indices.push_back("b");
  d2.indices.push_back("c");
  LatDim d3;
  d3.name = "d3";
  d3.size = 123;
  info.push_back(d1);
  info.push_back(d2);
  info.push_back(d3);
  lat_data_alloc(ld);
  for (long i = 0; i < ld.res.size(); ++i) {
    ld.res[i] = i;
  }
  return ld;
}

int main()
{
  LatData ld = mk_lat_data();
  ld.save("test.lat");
  ld.load("test.lat");
  ld.save("test.lat.2");
  return 0;
};
