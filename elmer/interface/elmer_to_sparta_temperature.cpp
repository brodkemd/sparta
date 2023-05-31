#include "elmer_to_sparta_temperature.h"


ElmerToSpartaTemperature::ElmerToSpartaTemperature(std::string _data_file, std::string _file_stem) : ElmerGridToSpartaSurf(_file_stem) {
    this->file_stem = _file_stem;
    this->data_file = _data_file;

    this->load_data();
}


void ElmerToSpartaTemperature::load_data() {
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream data_file;

    // opening the node file
    data_file.open(this->data_file);

    // temporary vector to store the split strings
    std::vector<std::string> v;
    // std::pair<std::string, std::vector<std::string>> temp;

    if (data_file.is_open()) {
        // reading the file
        while (data_file) {
            // getting the latest line
            std::getline(data_file, line);

            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));

            // making sure the size is right
            if (v.size() == 1) {
                // adding the data to the class list 
                this->data.push_back(v);
            }
        }
    } else {
        error(FLERR, "node file did not open");
    }

    // Close the file
    data_file.close();

}


void ElmerToSpartaTemperature::make_sparta_surf_data() {

    std::ofstream surf_data_file;
    surf_data_file.open(this->file_stem + ".surf.data");

    surf_data_file << "ITEM: TIMESTEP\n";
    surf_data_file << "NAN\n";
    surf_data_file << "ITEM: NUMBER OF SURFS\n";
    surf_data_file << "NAN\n";// std::to_string(tris.size())//(len(tris)),
    surf_data_file << "ITEM: BOX BOUNDS NAN NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "ITEM: SURFS id c\n";

}