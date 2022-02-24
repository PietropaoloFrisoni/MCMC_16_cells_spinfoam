#include "hash_21j_symbols.h"

int main(int argc, char **argv)
{

   if (argc < 4)
   {
      fprintf(stderr, "Usage: %s [DSPIN] [HASH_TABLES_STORE_PATH] [FASTWIG_TABLES_PATH] [VERBOSITY]\n", argv[0]);
      exit(EXIT_FAILURE);
   }

   std::cout << "\n\nHashing of 21j symbols...\n\n"
             << std::endl;

   int dspin_assigned = std::stoi(argv[1]);
   std::string hashed_tables_path_assigned = argv[2];
   std::string fastwig_tables_folder = argv[3];
   int verbosity = std::stoi(argv[4]);

   init(fastwig_tables_folder, verbosity);

   Hash_21j_symbols(hashed_tables_path_assigned, dspin_assigned);

   release();
}