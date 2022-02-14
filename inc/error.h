#pragma once

/**********************************************************************/

////////////////////////////////////////////////////////////////////////
// Functions to handle exceptions.
////////////////////////////////////////////////////////////////////////

// Prints the file and line of the error location, an optional message
// and then interrupts program.
#define error(format, ...) {\
        fprintf(stderr, "16 cell spinfoam: error at file %s, line %d.\n", __FILE__, __LINE__);\
        fprintf(stderr, "           : ");\
        fprintf(stderr, format, ##__VA_ARGS__);\
        fprintf(stderr, "\n");\
        exit(EXIT_FAILURE);\
        }

// Prints the file and line of the error location and an optional message.
#define warning(format, ...) {\
        fprintf(stderr, "16 cell spinfoam: warning at file %s, line %d.\n", __FILE__, __LINE__);\
        fprintf(stderr, "           : ");\
        fprintf(stderr, format, ##__VA_ARGS__);\
        fprintf(stderr, "\n");\
        }


/**********************************************************************/