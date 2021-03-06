#include <lat-io.h>
#include <iog_class.h>

inline LatData mk_lat_data()
{
  LatData ld;
  LatInfo& info = ld.info;
  LatDim d1;
  d1.name = "t";
  d1.size = 10;
  d1.indices.push_back("3");
  d1.indices.push_back("7");
  d1.indices.push_back("2");
  const LatDim d2 = lat_dim_string("conf", qutils::make_array<std::string>("a", "123", "c"));
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

  general_data_base iog("test.iog");
  convert_from_lat(ld,iog);
  iog.print();
  iog.save();
  convert_to_lat(iog,ld);
  ld.save("test.lat.3");
  return 0;
};
