#include "fileFormat.h"
#include <string>
using namespace KHMNDU005;

int main(int argc, char** argv)
{
    //parameters
    std::string outfile;
    std::string infile;
    
    infile  =std::string( argv[1]);
    outfile = std::string(argv[2]);
    
    fileFormat ff(infile, outfile);
    
    ff.read_write();
    
    return 0;


}