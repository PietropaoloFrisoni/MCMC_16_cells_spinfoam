#include "hash_21j_symbols.h"

int main(int argc, char **argv)
{

   if (argc < 3)
   {
      fprintf(stderr, "\n\nError: too few arguments.\nRun the script as: %s [DSPIN] [HASH_TABLES_STORE_PATH] [FASTWIG_TABLES_PATH]\n", argv[0]);
      exit(EXIT_FAILURE);
   }

   int dspin_assigned = std::stoi(argv[1]);
   std::string hashed_tables_path_assigned = argv[2];
   std::string fastwig_tables_folder = argv[3];

   init(fastwig_tables_folder);

   Hash_21j_symbols(hashed_tables_path_assigned, dspin_assigned);

   release();
}