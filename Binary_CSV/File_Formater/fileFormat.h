#ifndef FILEFORMAT_H_
#define FILEFORMAT_H_
#include <string>
#include <vector>



namespace KHMNDU005
{
    class fileFormat
    {
         public:
             
             
             //input file name
             std::string inputfile;
             //output file name
             std::string outputfile;
             
             fileFormat(std::string inputFile, std::string outputFile);
             ~fileFormat();
             //to hold a row of data
            
             
             //method to read from input and write to output csv file
             void  read_write();
             
             //method to compare time instances
             bool sameTime(std::string timeA, std::string timeB);
             
             
    };
    


}
#endif /* FILEFORMAT_H_*/