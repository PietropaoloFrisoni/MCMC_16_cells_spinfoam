#include <iostream>


class Chain {      

  public:       

    int dspin;
    int length;
    int * chain;
    int * indices;
    int burnin;
    double * amplitudes;
    double sigma;
    std::string store_path;

    Chain(int dspin_assigned, int length_assigned, double sigma_assigned, int burnin_assigned, std::string store_path_assigned) {
          dspin = dspin_assigned;
          indices = new int [16];
          length = length_assigned;  
          chain = new int [length];
          sigma = sigma_assigned;
          burnin = burnin_assigned;
          store_path = store_path_assigned;
          std::cout << "chain built with length " <<  length  << std::endl;
         };        

    ~Chain(){      
           delete [] chain;
           delete [] indices;
           std::cout << "chain destroyed " << std::endl; 
          };

};


int main () {

Chain catenazza(10, 10, 10, 10, "intucul");

std::cout << catenazza.dspin << std::endl;  

catenazza.dspin = 4;

std::cout << catenazza.dspin << std::endl;   

}