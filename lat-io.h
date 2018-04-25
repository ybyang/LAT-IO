#pragma once

#include <string>
#include <vector>

struct LatDim
{
  std::string name;
  long size;
  std::vector<std::string> indices;
};

typedef std::vector<LatDim> LatInfo;

struct LatData
{
  LatInfo info;
  std::vector<double> res;
};

