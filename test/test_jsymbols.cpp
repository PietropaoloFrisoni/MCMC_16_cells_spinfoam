#include <string>
#include "fastwigxj.h"
#include "wigxjpf.h"
#include <iostream>

#include "jsymbols.h"

using namespace std;

int main(int argc, char *argv[])
{

std::cout << "Testing non threaded version" << std::endl;

std::string fastwing_tables_folder = argv[1] ;
 
test_fastwig(fastwing_tables_folder);  
  
}
