#include "../include/file_tools.h"

#include <filesystem>

namespace fs = std::filesystem;

std::vector<std::string> get_files(const char *folder_path, const char *ext)
{
    std::vector<std::string> files;

    for (const auto& file : fs::directory_iterator(folder_path)) 
    {
        if (fs::is_regular_file(file.status()) && file.path().extension() == ext)
        {
            files.push_back(file.path().filename().string());
        }
    }

    return files;
}

void check_for_directories(const char *folder_path, std::vector<std::string> &directories)
{
    std::string path("");
    for(const auto& dir : directories)
    {
        path = folder_path + dir;
        if(!fs::exists(path))
        {
            if(fs::create_directory(path))
            {
                std::cout << "Directory " << dir << " has been created in " << path << std::endl;
            }
            else
            {
                std::cerr << "Failed to create directory " << dir << " in " << path << std::endl;
            }
        }
    }
}

void check_for_file(const char *path)
{
    if(!fs::exists(path))
    {
        std::ofstream file(path);
        if(file)
        {
            std::cout << "File has been created in " << path << std::endl;
        }
        else
        {
            std::cerr << "Failed to create file in " << path << std::endl;
        }
    }
}
