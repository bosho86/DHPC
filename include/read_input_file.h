/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   read_input_file.h
 * Author: hamilton
 *
 * Created on November 3, 2019, 7:30 AM
 */

#ifndef READ_INPUT_FILE_H
#define READ_INPUT_FILE_H

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

template <typename T>
class readInputFile{

public:
    
    readInputFile(){
    };
    
    std::vector< std::vector<T> > readFileToVector(const char *filename)
    {
        bool Val_bool = true;
        std::ifstream source;
        source.open(filename);
        std::string line;
        std::vector< std::vector<T> > voutter;

        while (std::getline(source, line))
        {
            std::vector<T> vinner;
            std::stringstream stream(line);
            
            while(Val_bool) {
                T n;
                stream >> n;
                
                if(!stream)
                    break;
               
                if(n > -1.0){
                    
                    vinner.push_back(n);

                }
            }            
            
            voutter.push_back(vinner);
        }
        
        return voutter;
    }

    void displayVector(const std::vector< std::vector <T> > v, const int index)
    {
        for (int i=0; i < v[index].size(); i++){
            std::cout<<v[index][i]<< std::endl;
        }
    }

};



#endif /* READ_INPUT_FILE_H */

