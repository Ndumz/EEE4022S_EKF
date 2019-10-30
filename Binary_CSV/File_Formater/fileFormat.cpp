#include "fileFormat.h"
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>

//start

namespace KHMNDU005
{
    fileFormat::fileFormat(std::string infile, std::string outfile):inputfile(infile),outputfile(outfile){
            
    }
    
    fileFormat::~fileFormat(){}
    
    void fileFormat::read_write()
    {
        std::ifstream in(inputfile);
        std::string line;
        std::string::size_type sz;
        
        //writing parameters
        std::ofstream out(outputfile);
       // out.open(outputfile);
        out<<"Time,XGyro,YGyro,ZGyro,XAccel,YAccel,ZAccel,XMag,YMag,ZMag,lat,long,height"<<"\n";
        
        
        //tells us how many lines we have gone through , resets at 3 lines
        int numLines = 0;
        while(std::getline(in,line))
        {
            std::stringstream lineStream(line);
            std::string cell;
            //current time
            std::string currTime = "0";
            
            std::getline(lineStream,cell,',');
           // std::cout<<"First input: "<<cell<<"\n";
            if(cell == "$PCHRS" || cell =="PCHRS")
            {
                std::vector<std::string> gyro;
                std::vector<std::string> accel;
                std::vector<std::string> mag;
                   
                std::getline(lineStream,cell,',');
                //std::cout<<"Sensor: "<<cell<<"\n";
                if(cell== "0")
                {
			
                    while(std::getline(lineStream,cell,','))
                    {
                        gyro.push_back(cell);
                    }
                }
                else if(cell== "1")
                {
		 std::getline(lineStream,cell,',');
                 while(std::getline(lineStream,cell,','))
                    {
                        accel.push_back(cell);
                    }
                }
                else if(cell== "2")
                {
		 std::getline(lineStream,cell,',');
                 while(std::getline(lineStream,cell,','))
                    {
                        mag.push_back(cell);
                    }    
                }
                //write line to file
                for(int i =0;i<gyro.size();i++)
                {
                    if((i<gyro.size()-1))
                    {
                    out<<gyro[i]<<",";
                    }
                    //std::cout<<"GYRO: "<<gyro.size()<<"\n";
                }
                
                for(int i =0 ; i<accel.size();i++)
                {
                    if(i<(accel.size()-1))
                    {
                    out<<accel[i]<<",";
                    }
                }
                
                for(int i =0 ; i< mag.size();i++)
                {
                     
                     if(i < mag.size()-2)
                     {
                         out<<mag[i]<<",";
                     }
                     if(i == (mag.size()-2))
                     {
                         out<<mag[i]<<",";
                     }
                }
                //out<<"\n";
                //gyro.erase();
                //accel.erase();
                //mag.erase();
                
                
            }
	if(cell=="$PCHRG")
	{
		std::string lat;
		std::string lng;
		std::string alt;

		std::getline(lineStream,cell,',');
		std::getline(lineStream,cell,',');
		lat = cell;
		out<<lat<<",";
		std::getline(lineStream,cell,',');
		lng = cell;
		out<<lng<<",";
		std::getline(lineStream,cell,',');
		alt = cell;
		out<<alt<<"\n";

			
	}
            //out<<"\n";
        }
        out.close();
    }
    


} 
