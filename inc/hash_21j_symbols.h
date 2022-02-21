#pragma once

#include <string>
#include <string.h>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

#include "fastwigxj.h"
#include "wigxjpf.h"
#include "common.h"
#include "error.h"

#include "phmap.h"
#include "phmap_dump.h"
#include "progressbar.h"

// load fastwigx and wigxjpf tables found in given path
void init(std::string root_folder_string, const int verbosity);

// release fastwigx and wigxjpf tables
void release();

// TODO add dynamic path
void Hash_21j_symbols(std::string hash_tables_store_path, const int tj);

// test fastwig tables
int test_fastwig(std::string fastwing_tables_folder);


struct MyKey
{
#if 1
  friend size_t hash_value(const MyKey &k)
  {
    return phmap::HashState().combine(0, k.key[0], k.key[1], k.key[2], k.key[3], k.key[4], k.key[5], k.key[6], k.key[7], k.key[8]);
  }

#else
  friend size_t hash_value(const MyKey &k)
  {
    uint64_t x = *((uint64_t *)k.key.data());
    return phmap::HashState().combine(0, x, k.key[8]);
  }
#endif
  bool operator==(const MyKey o) const { return key == o.key; }

  uint8_t &operator[](int idx) { return key[idx]; }

  std::array<uint8_t, 9> key;
};

