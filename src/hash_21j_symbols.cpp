#include "hash_21j_symbols.h"

char *FASTWIG_3J_TABLE_PATH;
char *FASTWIG_6J_TABLE_PATH;
char *FASTWIG_9J_TABLE_PATH;

char *DATA_ROOT;

// searches for all files in root folder
// with extension ".3j" or ".6j".
// The largest tables found are set for loading.
static void find_fastwig_tables(char *folder)
{

  DIR *d;
  struct dirent *dir;
  d = opendir(folder);

  if (d == NULL)
    error("error opening root folder");

  char table_3j_path[strlen(folder) + 256];
  char table_3j_path_toload[strlen(folder) + 256];
  long table_3j_size = 0L;

  char table_6j_path[strlen(folder) + 256];
  char table_6j_path_toload[strlen(folder) + 256];
  long table_6j_size = 0L;

  char table_9j_path[strlen(folder) + 256];
  char table_9j_path_toload[strlen(folder) + 256];
  long table_9j_size = 0L;

  char *fname;
  while ((dir = readdir(d)) != NULL)
  {

    fname = dir->d_name;
    char *dot = strrchr(fname, '.');

    // 3j
    if (dot && !strcmp(dot, ".3j"))
    {

      strcpy(table_3j_path, folder);
      strcat(table_3j_path, "/");
      strcat(table_3j_path, fname);

      // check size
      struct stat st;
      stat(table_3j_path, &st);

      if (st.st_size > table_3j_size)
      {

        // (larger) found, set path and size
        strcpy(table_3j_path_toload, folder);
        strcat(table_3j_path_toload, "/");
        strcat(table_3j_path_toload, fname);

        table_3j_size = st.st_size;
      }
    }

    // 6j
    if (dot && !strcmp(dot, ".6j"))
    {

      strcpy(table_6j_path, folder);
      strcat(table_6j_path, "/");
      strcat(table_6j_path, fname);

      // check size
      struct stat st;
      stat(table_6j_path, &st);

      if (st.st_size > table_6j_size)
      {

        // (larger) found, set path and size
        strcpy(table_6j_path_toload, folder);
        strcat(table_6j_path_toload, "/");
        strcat(table_6j_path_toload, fname);

        table_6j_size = st.st_size;
      }
    }

    // 9j
    if (dot && !strcmp(dot, ".9j"))
    {

      strcpy(table_9j_path, folder);
      strcat(table_9j_path, "/");
      strcat(table_9j_path, fname);

      // check size
      struct stat st;
      stat(table_9j_path, &st);

      if (st.st_size > table_9j_size)
      {

        // (larger) found, set path and size
        strcpy(table_9j_path_toload, folder);
        strcat(table_9j_path_toload, "/");
        strcat(table_9j_path_toload, fname);

        table_9j_size = st.st_size;
      }
    }
  }

  closedir(d);

  if (table_3j_size == 0)
    error("no 3j table found in root folder");
  if (table_6j_size == 0)
    error("no 6j table found in root folder");
  if (table_9j_size == 0)
    error("no 9j table found in root folder");

  FASTWIG_3J_TABLE_PATH = strdup(table_3j_path_toload);
  FASTWIG_6J_TABLE_PATH = strdup(table_6j_path_toload);
  FASTWIG_9J_TABLE_PATH = strdup(table_9j_path_toload);
}

static inline void load_fastwig_tables()
{
  fastwigxj_load(FASTWIG_3J_TABLE_PATH, 3, NULL);
  fastwigxj_load(FASTWIG_6J_TABLE_PATH, 6, NULL);
  fastwigxj_load(FASTWIG_9J_TABLE_PATH, 9, NULL);
}

static inline void unload_fastwig_tables()
{
  fastwigxj_unload(3);
  fastwigxj_unload(6);
  fastwigxj_unload(9);
}

void init(std::string root_folder_string, int verbosity)
{

  char *root_folder = &root_folder_string[0];

  std::cout << root_folder << std::endl;

  // check root folder is accessible
  DIR *d = opendir(root_folder);
  if (d == NULL)
    error("error opening root folder");
  closedir(d);

  DATA_ROOT = strdup(root_folder);

  // initialize wigxjpf
  wig_table_init(2 * 100, 9);

  // load fastwig tables
  find_fastwig_tables(DATA_ROOT);
  load_fastwig_tables();
}

void release()
{

  // wigxjpf
  wig_table_free();

  // fastwigxj
  unload_fastwig_tables();

  // free paths
  free(DATA_ROOT);
  free(FASTWIG_6J_TABLE_PATH);
  free(FASTWIG_3J_TABLE_PATH);
  free(FASTWIG_9J_TABLE_PATH);
}

void Hash_21j_symbols(const int tj)
{

  // char* store_path_assigned = &store_path_assigned_string[0];

  // build path for output table
  //   char filename[128];

  //   sprintf(filename, "/Hashed_21j_symbols_dspin_%d", tj);
  //   strcat(store_path_assigned, filename);

  // initialize table
  using Hash = phmap::flat_hash_map<MyKey, double>;
  Hash h{};

  // loop over all possible boundary intw
  // TODO: consider parallelize the code

  uint8_t ti1, ti2, ti3, ti4; // boundary intw
  uint8_t tb1, tb2, tb3;      // blue spins
  uint8_t tp1, tp2;           // purple spins
  uint8_t tg1, tg2, tl;       // internal virtual spins

  int tb1_min, tb1_max;
  int tb2_min, tb2_max;
  int tb3_min, tb3_max;
  int tp1_min, tp1_max;
  int tp2_min, tp2_max;
  int tg1_min, tg1_max;
  int tg2_min, tg2_max;
  int tl_min, tl_max;

  double sj1, sj2, sj3, sj4, sj5, sj6, sj7, nj; // loaded symbols
  double phi, pho;                              // phases
  double sji, ci, yi, ti;                       // building steps for internal symbols
  double sym, c, y, t;
  double df;

  int ti_min = 0;
  int ti_max = 2 * tj;

  // rough estimate for the number of steps
  int steps = DIM(tj) * DIM(tj) * DIM(tj) * DIM(tj) *
              (DIM(tj) + 1) * (DIM(tj) + 1) * (DIM(tj) + 1) * (DIM(tj) + 1) *
              (DIM(tj) + 2);
  int pstep = (int)ceil((double)steps / 100.0);
  int done = 0;

  wig_temp_init(2 * 100);

  printf("Maximum amount is roughly %d symbols, starting... \n", steps);

  for (ti1 = ti_min; ti1 <= ti_max; ti1 += 2)
  {
    for (ti2 = ti_min; ti2 <= ti_max; ti2 += 2)
    {
      for (ti3 = ti_min; ti3 <= ti_max; ti3 += 2)
      {
        for (ti4 = ti_min; ti4 <= ti_max; ti4 += 2)
        {

          tb1_min = 0;
          tb1_max = 2 * tj;

          for (tb1 = tb1_min; tb1 <= tb1_max; tb1 += 2)
          {

            tb2_min = abs(tb1 - tj);
            tb2_max = tb1 + tj;

            for (tb2 = tb2_min; tb2 <= tb2_max; tb2 += 2)
            {

              tb3_min = abs(tb2 - tj);
              tb3_max = tb2 + tj;

              for (tb3 = tb3_min; tb3 <= tb3_max; tb3 += 2)
              {

                tp1_min = 0;
                tp1_max = 2 * tj;

                for (tp1 = tp1_min; tp1 <= tp1_max; tp1 += 2)
                {

                  tp2_min = abs(tp1 - tj);
                  tp2_max = tp1 + tj;

                  for (tp2 = tp2_min; tp2 <= tp2_max; tp2 += 2)
                  {

                    // now internal sums
                    sym = c = 0.0;

                    tg1_min = abs(ti2 - tj);
                    tg1_max = ti2 + tj;

                    for (tg1 = tg1_min; tg1 <= tg1_max; tg1 += 2)
                    {

                      tg2_min = max(abs(tg1 - tj), abs(tb2 - tj));
                      tg2_max = min(tg1 + tj, tb2 + tj);

                      for (tg2 = tg2_min; tg2 <= tg2_max; tg2 += 2)
                      {

                        tl_min = max(abs(ti4 - tj), abs(tg2 - tj));
                        tl_max = min(ti4 + tj, tg2 + tj);

                        sji = ci = 0.0;

                        for (tl = tl_min; tl <= tl_max; tl += 2)
                        {

                          // load innermost symbols
                          sj5 = fw6jja(tj, ti3, tj, tj, ti4, tl);
                          sj6 = fw6jja(ti4, tj, tl, tb1, tj, tj);
                          sj7 = fw6jja(tl, tj, tg2, tb2, tj, tb1);
                          nj = fw9jja(tl, tj, ti3, tj, ti2, tj, tg2, tg1, tj);

                          phi = real_negpow(
                              2 * tl + (tj + tj + ti4) + (ti4 + tj + tl) + // sj5
                              (ti4 + tj + tj) +                            // sj6
                              (tl + tj + tb1) + (tg2 + tj + tb2) +         // sj7
                              2 * tl + (ti3 + tj + tj)                     // nj
                          );

                          df = DIM(tl);

                          comp_sum(df * phi * sj5 * sj6 * sj7 * nj, sji, ci, yi, ti);

                        } // tl

                        // build 'upper' symbols
                        sj1 = fw6jja(ti2, tj, tj, ti1, tj, tg1);
                        sj2 = fw6jja(tj, tp1, tj, tj, ti1, tg1);
                        sj3 = fw6jja(tj, tb3, tp2, tj, tg2, tb2);
                        sj4 = fw6jja(tj, tp2, tp1, tj, tg1, tg2);

                        pho = real_negpow(
                            2 * tg1 + (ti1 + tj + tj) + (ti2 + tj + tg1) +                    // sj1
                            2 * tg1 + (ti1 + tj + tg2) + (tg2 + tp1 + tj) + (ti1 + tj + tj) + // sj2
                            2 * tg2 + 2 * tp2 + (tg2 + tj + tb2) + (tb2 + tb3 + tj) +         // sj3
                            2 * tp2 + (tg1 + tp1 + tj)                                        // sj4
                        );

                        df = DIM(tg1) * DIM(tg2);

                        comp_sum(df * pho * sj1 * sj2 * sj3 * sj4 * sji, sym, c, y, t);

                      } // tg2
                    }   // tg1

                    // symbol computed, global phase absorbed in the partials

                    // h[MyKey{5, 5, 5}] = 11.33;

                    h[MyKey{ti1, ti2, ti3, ti4, tb1, tb2, tb3, tp1, tp2}] = sym;

                    done++;

                    if (done % pstep == 0)
                    {
                      printf("Done %2d%%... \n", done / pstep);
                    }

                  } // tp2
                }   // tp1

              } // tb3
            }   // tb2
          }     // tb1

        } // ti4
      }   // ti3
    }     // ti2
  }       // ti1

  // dump hashtable to disk

  phmap::BinaryOutputArchive ar_out("/home/frisus95/Scrivania/Final_project/data_folder/Hashed_21j_symbols_dspin_2");
  h.phmap_dump(ar_out);

  wig_temp_free();
}

int test_fastwig(std::string fastwig_tables_folder)
{

  // TODO Questo fa ALTAMENTE schifo, va migliorato.
  //  A fastwing bisogna passare un char * perchè non sa che cazzo sia string
  //  Dato che sto allocando sullo stack, non devo liberare memoria
  std::string stringa_3j_tmp = fastwig_tables_folder;
  std::string stringa_6j_tmp = fastwig_tables_folder;
  std::string stringa_9j_tmp = fastwig_tables_folder;

  std::string tmp_string_3j = stringa_3j_tmp.append("/test_table_18.3j");
  std::string tmp_string_6j = stringa_6j_tmp.append("/test_table_8.6j");
  std::string tmp_string_9j = stringa_9j_tmp.append("/test_hashed_3_6_12.9j");

  const char *fastwing_tables_folder_fullpath_3j = tmp_string_3j.c_str();
  const char *fastwing_tables_folder_fullpath_6j = tmp_string_6j.c_str();
  const char *fastwing_tables_folder_fullpath_9j = tmp_string_9j.c_str();

  double val3j, val6j, val9j;
  int i;

  printf("FASTWIGXJ C test program\n");

  /* Load tables produced during build test. */

  fastwigxj_load(fastwing_tables_folder_fullpath_3j, 3, NULL);
  fastwigxj_load(fastwing_tables_folder_fullpath_6j, 6, NULL);
  fastwigxj_load(fastwing_tables_folder_fullpath_9j, 9, NULL);

  /* For fallback of 9j to sum-of-6j-products when
   * symbols too large for table.
   *
   * And a dirty hack to not test this when quadmath/float128
   * is not available.
   */

  /* Large values to allow large 3j and 6j testing. */
  /* More normal use may need instead e.g. 2*100. */
  wig_table_init(2 * 10000, 9);
  wig_temp_init(2 * 10000);

  /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
   * and 2*m.  To be able to handle half-integer arguments.
   *
   * Also note that for 3j symbols, the third m is not given,
   * as it is fully redundant (unless trivial-0 symbols are wanted).
   *
   * (A version which takes the third m exists as fw3jja6).
   */

  val3j = fw3jja(2 * 5, 2 * 7, 2 * 5,
                 2 * (-3), 2 * 5);

  printf("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja(2 * 10, 2 * 15, 2 * 10,
                 2 * (-3), 2 * 12);

  printf("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val3j = fw3jja6(2 * 5, 2 * 7, 2 * 5,
                  2 * (-3), 2 * 5, 2 * (-2));

  printf("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

  val3j = fw3jja6(2 * 10, 2 * 15, 2 * 10,
                  2 * (-3), 2 * 12, 2 * (-9));

  printf("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

  val6j = fw6jja(2 * 3, 2 * 4, 2 * 2,
                 2 * 2, 2 * 2, 2 * 3);

  printf("6J{ 3   4   2;  2   2   3}:      %#25.15f\n", val6j);

  val6j = fw6jja(2 * 10, 2 * 15, 2 * 10,
                 2 * 7, 2 * 7, 2 * 9);

  printf("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

  for (i = 0; i < 3; i++)
  {
    /* Test with dynamic tables on second and third round. */

    /* Normally, one would set up the tables before evaluating
     * any symbols.
     */

    if (i == 1)
    {
      fastwigxj_dyn_init(3, 128);
      fastwigxj_dyn_init(6, 128);
      fastwigxj_dyn_init(9, 128);
    }

    val9j = fw9jja(1, 2, 3,
                   2, 3, 5,
                   3, 3, 6);

    printf("9J{0.5 1 1.5; 1 1.5 2.5; 1.5 1.5 3}:%#22.15f\n", val9j);

    /* The following symbol is too large for the C14N routine. */

    val9j = fw9jja(140, 140, 140,
                   140, 140, 140,
                   140, 140, 140);

    printf("9J{70 70 70; 70 70 70; 70 70 70}:%#25.15f\n", val9j);

    /* The following symbol is not in the 9j table, and will
     * automatically go via 6j fallback (through float128).
     */

    val9j = fw9jja(1, 2, 3,
                   4, 6, 8,
                   3, 6, 9);

    printf("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);

    /* This has arguments to become too large for the intermediate
     * index variables in the Rasch-Yu routine, so wont insert into
     * the dynamic table.
     */

    val3j = fw3jja6(2 * 8000, 2 * 8000, 2 * 8000,
                    2 * 1, 2 * 2, 2 * (-3));

    printf("3J(8000 8000 8000; 1   2  -3):   %#25.15f\n", val3j);

    /* Not exceeding 64-bit Rasch-Yu index. */

    val3j = fw3jja6(2 * 100, 2 * 100, 2 * 100,
                    2 * 1, 2 * 2, 2 * (-3));

    printf("3J(100 100 100; 1   2  -3):      %#25.15f\n", val3j);

    /* Not exceeding 64-bit Rasch-Yu index. */

    val6j = fw6jja(2 * 300, 2 * 400, 2 * 200,
                   2 * 200, 2 * 200, 2 * 300);

    printf("6J{300 400 200;  200 200 300}:   %#25.15f\n", val6j);

    /* Exceeding 64-bit Rasch-Yu index. */

    val6j = fw6jja(2 * 3000, 2 * 4000, 2 * 2000,
                   2 * 2000, 2 * 2000, 2 * 3000);

    printf("6J{3000 4000 2000;  2000 2000 3000}:%#22.15f\n", val6j);

    /* On last round, release values. */

    if (i == 2)
    {
      fastwigxj_dyn_free(3);
      fastwigxj_dyn_free(6);
      fastwigxj_dyn_free(9);
    }
  }

  {
    /* Call the lookup functions with argument lists instead.
     * The lists are not destroyed, so can be partially reused.
     */
    int arg3j[5], arg6j[6], arg9j[9];

    arg3j[0] = 2 * 5;
    arg3j[1] = 2 * 7;
    arg3j[2] = 2 * 5;
    arg3j[3] = 2 * (-3);
    arg3j[4] = 2 * 5;

    val3j = fw3jjl(arg3j);

    printf("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

    arg3j[0] = 2 * 10;
    arg3j[1] = 2 * 15;
    arg3j[2] = 2 * 10;
    arg3j[3] = 2 * (-3);
    arg3j[4] = 2 * 12;

    val3j = fw3jjl(arg3j);

    printf("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

    arg3j[0] = 2 * 11; /* Partial reuse. */

    val3j = fw3jjl(arg3j);

    printf("3J(11  15  10; -3  12  -9):      %#25.15f\n", val3j);

    arg6j[0] = 2 * 3;
    arg6j[1] = 2 * 4;
    arg6j[2] = 2 * 2;
    arg6j[3] = 2 * 2;
    arg6j[4] = 2 * 2;
    arg6j[5] = 2 * 3;

    val6j = fw6jjl(arg6j);

    printf("6J{ 3   4   2;  2   2   3}:      %#25.15f\n", val6j);

    arg6j[0] = 2 * 10;
    arg6j[1] = 2 * 15;
    arg6j[2] = 2 * 10;
    arg6j[3] = 2 * 7;
    arg6j[4] = 2 * 7;
    arg6j[5] = 2 * 9;

    val6j = fw6jjl(arg6j);

    printf("6J{10  15  10;  7   7   9}:      %#25.15f\n", val6j);

    arg9j[0] = 1;
    arg9j[1] = 2;
    arg9j[2] = 3;
    arg9j[3] = 4;
    arg9j[4] = 6;
    arg9j[5] = 8;
    arg9j[6] = 3;
    arg9j[7] = 6;
    arg9j[8] = 9;

    val9j = fw9jjl(arg9j);

    printf("9J{0.5 1 1.5; 2 3 4; 1.5 3 4.5}: %#25.15f\n", val9j);
  }

  {
    /* By using separate canonicalisation (and prefetch) and value
     * retrieval, several lookups can be interleaved.
     */

    int arg3j_1[5];
    int arg3j_2[5];
    int arg3j_3[5];
    uint64_t x_1;
    uint64_t x_2;
    uint64_t x_3;

    arg3j_1[0] = 2 * 5;
    arg3j_1[1] = 2 * 7;
    arg3j_1[2] = 2 * 5;
    arg3j_1[3] = 2 * (-3);
    arg3j_1[4] = 2 * 5;

    fw3jj_canon(arg3j_1, &x_1);
    fw3jj_prefetch(x_1);

    arg3j_2[0] = 2 * 10;
    arg3j_2[1] = 2 * 15;
    arg3j_2[2] = 2 * 10;
    arg3j_2[3] = 2 * (-3);
    arg3j_2[4] = 2 * 12;

    fw3jj_canon(arg3j_2, &x_2);
    fw3jj_prefetch(x_2);

    arg3j_3[0] = 2 * 11;
    arg3j_3[1] = 2 * 15;
    arg3j_3[2] = 2 * 10;
    arg3j_3[3] = 2 * (-3);
    arg3j_3[4] = 2 * 12;

    fw3jj_canon(arg3j_3, &x_3);
    fw3jj_prefetch(x_3);

    /* The list of original arguments must be provided at lookup
     * also, as they are used in case the value is not in the table
     * and a fallback calculation is required.
     */

    val3j = fw3jj_get(arg3j_1, x_1);

    printf("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);

    val3j = fw3jj_get(arg3j_2, x_2);

    printf("3J(10  15  10; -3  12  -9):      %#25.15f\n", val3j);

    val3j = fw3jj_get(arg3j_3, x_3);

    printf("3J(11  15  10; -3  12  -9):      %#25.15f\n", val3j);
  }

  /* Print statistics. */

  printf("\n");
  fastwigxj_print_stats();

  /* Remove tables from memory. */

  fastwigxj_unload(3);
  fastwigxj_unload(6);
  fastwigxj_unload(9);

  /* Release WIGXJPF memory. */

  wig_temp_free();
  wig_table_free();

  return 0;
}
