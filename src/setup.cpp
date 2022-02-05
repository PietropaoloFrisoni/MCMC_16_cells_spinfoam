#include <iostream>
#include <string.h>
#include <dirent.h> 
#include <sys/stat.h>

#include "wigxjpf.h"
#include "fastwigxj.h"
#include "utilities.h"
#include "error.h"
#include "common.h"


// path of found fastwig tables
char* FASTWIG_3J_TABLE_PATH;
char* FASTWIG_6J_TABLE_PATH;


// searches for all files in root folder 
// with extension ".3j" or ".6j". 
// The largest tables found are set for loading.
static void find_fastwig_tables(char* folder) {
  
    DIR *d;
    struct dirent *dir;
    d = opendir(folder);

    if (d == NULL) error("error opening root folder");

    char table_3j_path[strlen(folder) + 256];
    char table_3j_path_toload[strlen(folder) + 256];
    long table_3j_size = 0L;
    
    char table_6j_path[strlen(folder) + 256];
    char table_6j_path_toload[strlen(folder) + 256];
    long table_6j_size = 0L;

    char* fname;
    while ((dir = readdir(d)) != NULL) {

        fname = dir->d_name;
        char *dot = strrchr(fname, '.');

        // 3j
        if (dot && !strcmp(dot, ".3j")) {

            strcpy(table_3j_path, folder);
            strcat(table_3j_path, "/");
            strcat(table_3j_path, fname);

            // check size
            struct stat st;
            stat(table_3j_path, &st);

            if (st.st_size > table_3j_size) {

                // (larger) found, set path and size
                strcpy(table_3j_path_toload, folder);
                strcat(table_3j_path_toload, "/");
                strcat(table_3j_path_toload, fname);

                table_3j_size = st.st_size;

            }

        }

        // 6j
        if (dot && !strcmp(dot, ".6j")) {

            strcpy(table_6j_path, folder);
            strcat(table_6j_path, "/");
            strcat(table_6j_path, fname);

            // check size
            struct stat st;
            stat(table_6j_path, &st);

            if (st.st_size > table_6j_size) {

                // (larger) found, set path and size
                strcpy(table_6j_path_toload, folder);
                strcat(table_6j_path_toload, "/");
                strcat(table_6j_path_toload, fname);

                table_6j_size = st.st_size;

            }

        }

    }
    closedir(d);

    if (table_3j_size == 0) error("no 3j table found in root folder");
    if (table_6j_size == 0) error("no 6j table found in root folder");

    FASTWIG_3J_TABLE_PATH = strdup(table_3j_path_toload);
    FASTWIG_6J_TABLE_PATH = strdup(table_6j_path_toload);

}

static inline void load_fastwig_tables() {
    fastwigxj_load(FASTWIG_3J_TABLE_PATH, 3, NULL);
    fastwigxj_load(FASTWIG_6J_TABLE_PATH, 6, NULL);
}

static inline void unload_fastwig_tables() {
    fastwigxj_unload(3);
    fastwigxj_unload(6);
}
