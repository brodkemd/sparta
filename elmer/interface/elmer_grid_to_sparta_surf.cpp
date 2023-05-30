#include "elmer_grid_to_sparta_surf.h"


ElmerGridToSpartaSurf::ElmerGridToSpartaSurf(std::string _file_stem) {
    this->file_stem = _file_stem;

    this->load_boundary();
    this->load_nodes();
}
        

void ElmerGridToSpartaSurf::load_boundary() {
    std::cout << "Loading Boundary\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream boundary_file;

    // opening the boundary file
    boundary_file.open(this->file_stem + ".boundary");

    // temporary vector to store the split strings
    std::vector<std::string> v;

    if (boundary_file.is_open()) {
        while (boundary_file) {
            std::getline(boundary_file, line);
            boost::split(v, line, boost::is_any_of(" "));
            this->boundary_data.push_back(v);
        }
    }

    // Close the file
    boundary_file.close();
}


void ElmerGridToSpartaSurf::load_nodes() {
    std::cout << "Loading Nodes\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream node_file;

    // opening the node file
    node_file.open(this->file_stem + ".boundary");

    // temporary vector to store the split strings
    std::vector<std::string> v;
    std::pair<std::string, std::vector<std::string>> temp;

    if (node_file.is_open()) {
        while (node_file) {
            std::getline(node_file, line);
            boost::split(v, line, boost::is_any_of(" "));
            temp.first = v[0];
            v.erase(v.begin());
            temp.second = v;
            this->node_data.push_back(temp);
        }
    }

    // Close the file
    node_file.close();

}


void ElmerGridToSpartaSurf::make_sparta_surf() {
    std::vector<std::pair<std::string, std::vector<std::string>>> tris;
    std::vector<std::string> boundary_element_id_list;
    std::pair<std::string, std::vector<std::string>> temp_pair;
    std::vector<std::string> temp_vec;

    for (int i = 0; i < this->boundary_data.size(); i++) {
        temp_vec = this->boundary_data[i];
        temp_pair.first = temp_vec[0];
        temp_vec.erase(temp_vec.begin(), temp_vec.begin() + 5);
        temp_pair.second = temp_vec;

        tris.push_back(temp_pair);
        for (std::string id : temp_vec) {
            boundary_element_id_list.push_back(id);
        }

    }
    
    for (std::pair<std::string, std::vector<std::string>> val : tris) {
        std::cout << val
    }
}