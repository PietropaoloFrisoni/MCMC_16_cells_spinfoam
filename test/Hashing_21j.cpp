#include "hash_21j_symbols.h"

int main(int argc, char *argv[])
{

   std::cout << "\n\nTesting hashing of 21j symbols...\n\n" << std::endl;

   // in general I want these two strings to be different!
   std::string fastwig_tables_folder = argv[1];
   std::string store_path_assigned = argv[1];
   std::string hashed_tables_path_assigned = "./data_folder/hashed_21j";
   int dspin_assigned = 2;
   int verbosity = 2;

   init(fastwig_tables_folder, verbosity);

   Hash_21j_symbols(hashed_tables_path_assigned, dspin_assigned);

   release();
}