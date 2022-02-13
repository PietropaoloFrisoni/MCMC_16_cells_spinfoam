
#include <string>
#include <iostream>

#include "fastwigxj.h"
#include "wigxjpf.h"

#include "phmap.h"
#include "phmap_dump.h"


using phmap::flat_hash_map;

void lezzo(){
 // Create an unordered_map of three strings (that map to strings)
  flat_hash_map<std::string, std::string> email =
      {
          {"tom", "tom@gmail.com"},
          {"jeff", "jk@gmail.com"},
          {"jim", "jimg@microsoft.com"}};

  // Iterate and print keys and values
  for (const auto &n : email)
    std::cout << n.first << "'s email is: " << n.second << "\n";

  // Add a new entry
  email["bill"] = "bg@whatever.com";

  // and print it
  std::cout << "bill's email is: " << email["bill"] << "\n";


   wig_temp_init(2 * 10000);

  /* Note that the arguments to wig3jj, wig6jj and wig9jj are 2*j
   * and 2*m.  To be able to handle half-integer arguments.
   *
   * Also note that for 3j symbols, the third m is not given,
   * as it is fully redundant (unless trivial-0 symbols are wanted).
   *
   * (A version which takes the third m exists as fw3jja6).
   */

  double val3j = fw3jja(2 * 5, 2 * 7, 2 * 5,
                 2 * (-3), 2 * 5);

  printf("3J( 5   7   5; -3   5  -2):      %#25.15f\n", val3j);


}