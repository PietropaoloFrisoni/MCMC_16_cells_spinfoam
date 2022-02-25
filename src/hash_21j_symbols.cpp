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

void init(std::string root_folder_string)
{

  char *root_folder = &root_folder_string[0];

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

void Hash_21j_symbols(std::string hash_tables_store_path, const int tj)
{

  char tmp[1024];

  sprintf(tmp, "/hashed_21j_symbols_j_%.8g", ((double)(tj) / 2.0));

  hash_tables_store_path += std::string(tmp);

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

                    h[MyKey{ti1, ti2, ti3, ti4, tb1, tb2, tb3, tp1, tp2}] = sym;

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

  phmap::BinaryOutputArchive ar_out(&hash_tables_store_path[0]);
  h.phmap_dump(ar_out);

  wig_temp_free();
}