#include <iostream>
#include <getopt.h> // For getopt
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TCollection.h>
#include <TH1.h>

namespace
{
    std::string Indent(int level)
    {
        return std::string(level * 2, ' ');
    }

    void PrintHistogramInfo(TH1 *hist, int level)
    {
        if (!hist)
            return;

        std::cout << Indent(level) << "[HIST] " << hist->GetName()
                  << " | class=" << hist->ClassName()
                  << " | title=\"" << hist->GetTitle() << "\""
                  << " | nbinsX=" << hist->GetNbinsX()
                  << " | entries=" << hist->GetEntries()
                  << "\n";
    }

    void PrintObjectInfo(TObject *obj, int level)
    {
        if (!obj)
            return;

        if (obj->InheritsFrom(TH1::Class()))
        {
            PrintHistogramInfo(static_cast<TH1 *>(obj), level);
            return;
        }

        std::cout << Indent(level) << "[OBJ ] " << obj->GetName()
                  << " | class=" << obj->ClassName() << "\n";

        if (obj->InheritsFrom(TDirectory::Class()))
        {
            TDirectory *dir = static_cast<TDirectory *>(obj);
            TIter nextKey(dir->GetListOfKeys());
            TKey *subKey = nullptr;
            while ((subKey = static_cast<TKey *>(nextKey())))
            {
                TObject *subObj = subKey->ReadObj();
                std::cout << Indent(level + 1) << "[KEY ] " << subKey->GetName()
                          << " | class=" << subKey->GetClassName() << "\n";
                PrintObjectInfo(subObj, level + 2);
                delete subObj;
            }
            return;
        }

        if (obj->InheritsFrom(TCollection::Class()))
        {
            TCollection *coll = static_cast<TCollection *>(obj);
            TIter next(coll);
            TObject *subObj = nullptr;
            while ((subObj = next()))
            {
                PrintObjectInfo(subObj, level + 1);
            }
        }
    }
}

void print_keys() {
    // std::string filePath = "/home/sawan/check_k892/data/kstar/LHC22o_pass7/447406.root"; // Default file path
    std::string filePath = "/home/sawan/check_k892/macro/AnalysisMCpass1pass2.root"; // Default file path

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
    // Print all keys and recursively inspect objects under them
    while ((key = (TKey *)nextkey())) {
        std::cout << "[KEY ] " << key->GetName()
                  << " | class=" << key->GetClassName() << "\n";

        TObject *obj = key->ReadObj();
        PrintObjectInfo(obj, 1);
        delete obj;
    }

    file->Close();
    delete file;

}
