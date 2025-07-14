#include <iostream>
#include <getopt.h> // For getopt
#include <TFile.h>

void print_keys() {
    std::string filePath = "/home/sawan/check_k892/data/kstar/LHC22o_pass7/447406.root"; // Default file path

    // // Parse command line arguments
    // int option;
    // while ((option = getopt(argc, argv, "d:")) != -1) {
    //     switch (option) {
    //         case 'd':
    //             filePath = optarg; // Set file path to the provided argument
    //             break;
    //         default:
    //             std::cerr << "Usage: " << argv[0] << " [-d file_path]\n";
    //             return 1;
    //     }
    // }

    // Open the ROOT file
    TFile *file = TFile::Open(filePath.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filePath << "\n";
        return;
    }
    std::cout << "Successfully opened file: " << filePath << "\n";

    TIter nextkey(file->GetListOfKeys());
    TKey *key;
    // print all keys in the file
    while ((key = (TKey *)nextkey())) {
        std::cout << key->GetName() << "\n";
    }

    // Perform operations on the file...

    file->Close();
    delete file;

}
