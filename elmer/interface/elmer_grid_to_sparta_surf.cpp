#include "elmer_grid_to_sparta_surf.h"

/**
 * Constructor for the class
*/
ElmerGridToSpartaSurf::ElmerGridToSpartaSurf(std::string _file_stem) {
    // setting the file
    this->file_stem = _file_stem;

    // loading the boundary
    this->load_boundary();

    // loading the node data
    this->load_nodes();
}


/**
 * Loads boundary element ids and nodes from boundary file
*/
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
        unsigned int count = 1;
        // reading the file
        while (boundary_file) {
            // getting the latest line
            std::getline(boundary_file, line);

            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));

            if (v[4] != (std::string)"303") {
                error(FLERR, ((std::string)("element is not a triangle in boundary file at line " + std::to_string(count))).c_str());
            }

            // if the line contains a certain number of values, erase the first 5
            if (v.size() > 5) {
                v.erase(v.begin(), v.begin() + 5);
            }

            // adding the data to the class vector
            this->boundary_data.push_back(v);

            // incrementing line counter
            count++;
        }
    } else error(FLERR, "boundary file did not open");
    

    // Close the file
    boundary_file.close();
}


/**
 * Loads node locations from the node file
*/
void ElmerGridToSpartaSurf::load_nodes() {
    std::cout << "Loading Nodes\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream node_file;

    // opening the node file
    node_file.open(this->file_stem + ".nodes");

    // temporary vector to store the split strings
    std::vector<std::string> v;
    std::pair<std::string, std::vector<std::string>> temp;

    if (node_file.is_open()) {
        // reading the file
        while (node_file) {
            // getting the latest line
            std::getline(node_file, line);

            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));

            // setting the first element in the pair to the id
            temp.first = v[0];
            
            // erasing unneeded values
            v.erase(v.begin(), v.begin()+2);

            // setting the second to be the data
            temp.second = v;

            // adding the data to the class list
            this->node_data.push_back(temp);
        }
    } else {
        error(FLERR, "node file did not open");
    }

    // Close the file
    node_file.close();

}


/**
 * the main function in this class
 */
void ElmerGridToSpartaSurf::make_sparta_surf() {
    // stores triangle data
    std::vector<std::pair<std::string, std::vector<std::string>>> tris;

    // vector of boundary node ids, used to make sure they ascend with no gaps
    std::vector<unsigned int> boundary_element_id_list;

    // temporary variables
    std::pair<std::string, std::vector<std::string>> temp_pair;
    std::vector<std::string> temp_vec;

    // getting the boundary triangles and boundary node ids
    for (int i = 0; i < this->boundary_data.size(); i++) {
        // getting the vector from the class vector
        temp_vec = this->boundary_data[i];

        // setting parameters of the pair
        temp_pair.first = temp_vec[0];
        temp_pair.second = temp_vec;

        // adding to tris list
        tris.push_back(temp_pair);

        // adding each node id to the list to check for missing ids
        for (std::string id : temp_vec) {
            boundary_element_id_list.push_back(atoi(id.c_str()));
        }
    }

    // converting to a set sorts and filter duplicates from the vector
    std::set<int> s; // making a set
    for( unsigned i = 0; i < boundary_element_id_list.size(); ++i ) {
        s.insert(boundary_element_id_list[i]); // inserting id into list
    }

    // putting the data back into the vector
    boundary_element_id_list.assign(s.begin(), s.end());

    // making sure the ids increment by one in the boundary node ids
    for (int i = 0; i < boundary_element_id_list.size() - 1; i++) {
        if (boundary_element_id_list[i+1] - boundary_element_id_list[i] != 1) {
            std::cout << boundary_element_id_list[i+1] << " " << boundary_element_id_list[i] << "\n";
            error(FLERR, "Ids do not ascend, will need to modify this code");
        }
    }

    std::cout << "writing to: " << this->file_stem + ".surf" << "\n";

    // writing to file
    std::ofstream surf_file;
    surf_file.open(this->file_stem + ".surf");

    if (!(surf_file.is_open())) error(FLERR, "Surf file could not be opened");

    // the header
    surf_file << "# Surface element file written by SPARTA\n\n";
    surf_file << boundary_element_id_list.size() << " points\n";
    surf_file << tris.size() << " triangles\n\nPoints\n\n";

    // the point data
    for (std::pair<std::string, std::vector<std::string>> it : this->node_data) {
        surf_file << it.first << " ";
        for (int i = 0; i < it.second.size() - 1; i++) {
            surf_file << it.second[i] << " ";
        }
        surf_file << it.second.back() << "\n";
    }

    surf_file << "\nTriangles\n\n";

    // the triangles
    for (std::pair<std::string, std::vector<std::string>> val : tris) {
        surf_file << val.first << " ";
        for (int i = 0; i < val.second.size() - 1; i++) {
            surf_file << val.second[i] << " ";
        }
        surf_file << val.second.back() << "\n";
    }

    // closing the file
    surf_file.close();
}