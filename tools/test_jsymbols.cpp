#include "hash_21j_symbols.h"

using namespace std;

int main(int argc, char *argv[])
{

std::cout << "\n\nTesting fastwig tables...\n\n" << std::endl;

std::string fastwig_tables_folder = argv[1];
 
test_fastwig(fastwig_tables_folder);  

}
