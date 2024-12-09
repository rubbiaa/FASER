/* 
  FASERCAL File Mask class
  A. Rubbia September 2024 
*/

#ifndef _TFILEMASK_
#define _TFILEMASK_ 1

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class TFileMask : public TObject {

public:
    const char *filename;
    int run_number;
    std::vector<int> event_numbers;


    TFileMask() : filename(nullptr), run_number(0) {};

    TFileMask(int run, const char *filename) : filename(filename), run_number(run) {};

    void addEvent(int event) { event_numbers.push_back(event);};

    void Dump() const {
        std::string fname("FileMask_run" + std::to_string(run_number) + "_" + std::string(filename) + ".fmask");
        std::ofstream file_out(fname);
        if (!file_out.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }
        // Write the header
        file_out << "Run Number,Event Number\n";
        for (const auto& event : event_numbers) {
            file_out << run_number << "," << event << "\n";
        }
        file_out.close();
    };

};

#endif