#pragma once

#include <string>
#include <vector>
#include <cassert>

#include <endian.h>

#include <show.h>

const std::string lat_data_header = "#!/usr/bin/env lat-io-glimpse\n";

struct LatDim
{
  std::string name;
  long size = 0;
  std::vector<std::string> indices;
};

typedef std::vector<LatDim> LatInfo;

struct LatData
{
  LatInfo info;
  std::vector<double> res;
  //
  void load(const std::string& fn);
  void save(const std::string& fn) const;
};

inline long lat_data_size(const LatInfo& info)
{
  if (info.size() == 0) {
    return 0;
  }
  long total = 1;
  for (int i = 0; i < info.size(); ++i) {
    total *= info[i].size;
  }
  return total;
}

inline void lat_data_alloc(LatData& ld)
{
  ld.res.resize(lat_data_size(ld.info));
}

namespace qutils {

inline uint64_t flip_endian_64(uint64_t x)
{
  return
    ((x >> 56)) |
    ((x >> 40) & 0xFF00) |
    ((x >> 24) & 0xFF0000) |
    ((x >>  8) & 0xFF000000) |
    ((x <<  8) & 0xFF00000000) |
    ((x << 24) & 0xFF0000000000) |
    ((x << 40) & 0xFF000000000000) |
    ((x << 56));
}

inline void flip_endian_64(void* str, const size_t len)
{
  assert(0 == len % 8);
  uint64_t* p = (uint64_t*)str;
  for (size_t i = 0; i < len / 8; ++i) {
    p[i] = flip_endian_64(p[i]);
  }
}

inline bool is_big_endian()
{
#if defined(__BYTE_ORDER) && (__BYTE_ORDER != 0) && (__BYTE_ORDER == __BIG_ENDIAN)
  return true;
#else
  return false;
#endif
}

inline bool is_little_endian()
{
  return not is_big_endian();
}

inline void to_from_little_endian_64(void* str, const size_t len)
{
  assert(0 == len % 8);
  if (is_big_endian()) {
    flip_endian_64(str, len);
  }
}

inline void to_from_big_endian_64(void* str, const size_t len)
{
  assert(0 == len % 8);
  if (is_little_endian()) {
    flip_endian_64(str, len);
  }
}

inline std::string qgetline(FILE* fp)
{
  char* lineptr = NULL;
  size_t n = 0;
  if (getline(&lineptr, &n, fp) > 0) {
    std::string ret(lineptr);
    std::free(lineptr);
    return ret;
  } else {
    std::free(lineptr);
    return std::string();
  }
}

}

namespace qshow {

inline std::string show(const LatDim& dim)
{
  std::ostringstream out;
  out << ssprintf("\"%s\"[%ld]:", dim.name.c_str(), dim.size);
  for (long i = 0; i < dim.indices.size(); ++i) {
    out << ssprintf(" \"%s\"", dim.indices[i].c_str());
  }
  return out.str();
}

inline std::string show(const LatInfo& info)
{
  std::ostringstream out;
  out << ssprintf("ndim: %ld\n", info.size());
  for (int i = 0; i < info.size(); ++i) {
    out << show(info[i]) << "\n";
  }
  return out.str();
}

inline std::vector<std::string> split_into_lines(const std::string& str)
{
  const size_t len = str.length();
  std::vector<std::string> lines;
  size_t start = 0;
  size_t stop = 0;
  while (stop < len) {
    while (start < len && str[start] == '\n') {
      start += 1;
    }
    stop = start;
    while (stop < len && !(str[stop] == '\n')) {
      stop += 1;
    }
    if (stop > start) {
      lines.push_back(std::string(str, start, stop - start));
    }
    start = stop;
  }
  return lines;
}

inline bool parse_char(char& c, long& cur, const std::string& data)
{
  if (data.size() <= cur) {
    return false;
  } else {
    c = data[cur];
    cur += 1;
    return true;
  }
}

inline bool parse_string(std::string& str, long& cur, const std::string& data)
{
  char c;
  if (!parse_char(c, cur, data) or c != '"') {
    return false;
  } else {
    const long start = cur;
    char c;
    while (parse_char(c, cur, data) and c!= '"') {
    }
    str = std::string(data, start, cur - start - 1);
    return data[cur-1] == '"' && cur > start;
  }
}

inline bool parse_long(long& num, long& cur, const std::string& data)
{
  const long start = cur;
  char c;
  while (parse_char(c, cur, data)) {
    if ('0' > c or c > '9') {
      cur -= 1;
      break;
    }
  }
  if (cur <= start) {
    return false;
  } else {
    const std::string str = std::string(data, start, cur - start);
    num = read_long(str);
    return true;
  }
}

inline LatDim read_lat_dim(const std::string& str)
{
  LatDim dim;
  long cur = 0;
  char c;
  if (!parse_string(dim.name, cur, str)) {
    assert(false);
  } else if (!parse_char(c, cur, str) or c != '[') {
    assert(false);
  } else if (!parse_long(dim.size, cur, str)) {
    assert(false);
  } else if (!parse_char(c, cur, str) or c != ']') {
    assert(false);
  } else if (!parse_char(c, cur, str) or c != ':') {
    assert(false);
  } else {
    while (parse_char(c, cur, str)) {
      assert(c == ' ');
      std::string index;
      if (!parse_string(index, cur, str)) {
        assert(false);
      }
      dim.indices.push_back(index);
    }
  }
  return dim;
}

inline LatInfo read_lat_info(const std::string& str)
{
  LatInfo info;
  const std::vector<std::string> infos = split_into_lines(str);
  assert(infos.size() >= 1);
  const std::string ndim_prop = "ndim: ";
  assert(infos[0].compare(0, ndim_prop.size(), ndim_prop) == 0);
  assert(read_long(std::string(infos[0], ndim_prop.size())) == infos.size() - 1);
  for (int i = 1; i < infos.size(); ++i) {
    info.push_back(read_lat_dim(infos[i]));
  }
  return info;
}

}

inline void LatData::load(const std::string& fn)
{
  using namespace qshow;
  using namespace qutils;
  FILE* fp = fopen(fn.c_str(), "r");
  std::vector<char> check_line(lat_data_header.size(), 0);
  fread(check_line.data(), lat_data_header.size(), 1, fp);
  assert(std::string(check_line.data(), check_line.size()) == lat_data_header);
  std::vector<std::string> infos;
  infos.push_back(lat_data_header);
  while (infos.back() != "END_HEADER\n" && infos.back() != "") {
    infos.push_back(qgetline(fp));
  }
  std::ostringstream out;
  for (int i = 3; i < infos.size() - 1; ++i) {
    out << infos[i];
  }
  const std::string info_str = out.str();
  info = read_lat_info(info_str);
  lat_data_alloc(*this);
  assert(res.size() == lat_data_size(info));
  assert(res.size() * sizeof(double) == read_long(infos[2]));
  fread(res.data(), sizeof(double), res.size(), fp);
  to_from_little_endian_64(res.data(), res.size() * sizeof(double));
  fclose(fp);
}

inline void LatData::save(const std::string& fn) const
{
  using namespace qshow;
  using namespace qutils;
  FILE* fp = fopen(fn.c_str(), "w");
  assert(fp != NULL);
  const std::string data_size = ssprintf("data_size\n%ld\n", res.size() * sizeof(double));
  const std::string info_str = show(info);
  const std::string end_header = "END_HEADER\n";
  fwrite(lat_data_header.data(), lat_data_header.size(), 1, fp);
  fwrite(data_size.data(), data_size.size(), 1, fp);
  fwrite(info_str.data(), info_str.size(), 1, fp);
  fwrite(end_header.data(), end_header.size(), 1, fp);
  if (is_little_endian()) {
    fwrite(res.data(), sizeof(double), res.size(), fp);
  } else {
    std::vector<double> res_copy(res);
    assert(res_copy.size() == res.size());
    to_from_little_endian_64(res_copy.data(), res_copy.size() * sizeof(double));
    fwrite(res_copy.data(), sizeof(double), res_copy.size(), fp);
  }
  fclose(fp);
}

