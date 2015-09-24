#include "parameters.h"
#include <sstream>
#include <string>

char *parameter_mpi_device(int i) {
	    std::stringstream s;
	    s << "mpi_device_" << i;
	    char *a;
	    a = new char[s.str().size()+1];
	    strcpy(a,s.str().c_str());
	    return a;
}

