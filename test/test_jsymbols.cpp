#include <string>
#include <iostream>

#include "fastwigxj.h"
#include "wigxjpf.h"

#include "jsymbols.h"

using namespace std;

int main(int argc, char *argv[])
{

std::cout << "\n\nTESTING FASTWIG TABLES...\n\n" << std::endl;

std::string fastwig_tables_folder = argv[1];
 
test_fastwig(fastwig_tables_folder);  

  
}
