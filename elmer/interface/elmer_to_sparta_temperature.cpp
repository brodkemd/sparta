#include "elmer_to_sparta_temperature.h"


/**
 * The constructor for this class, also makes the parent class
*/
ElmerToSpartaTemperature::ElmerToSpartaTemperature(std::string _data_file, std::string _file_stem) : ElmerGridToSpartaSurf(_file_stem) {
    // setting the attributes
    this->file_stem = _file_stem;
    this->data_file = _data_file;
}

/**
 * loads data from the elmer output file
*/
void ElmerToSpartaTemperature::load_data() {
    // std::cout << "loading data\n";
    // Create a text string, which is used to output the text file
    std::string line;

    // Read from the text file
    std::ifstream data_file;

    // opening the node file
    data_file.open(this->data_file);

    // temporary vector to store the split strings
    std::vector<std::string> v;


    if (data_file.is_open()) {
        bool go = false;

        // reading the file
        while (data_file) {
            // getting the latest line
            std::getline(data_file, line);

            // trimming whitespaces off of the line
            boost::algorithm::trim(line);
            
            // filters empty lines
            if (!(line.length())) continue;

            // splitting the line at spaces
            boost::split(v, line, boost::is_any_of(" "));
            
            // if good to record data
            if (go) {
                if (v.size() == 1) {
                    // adding the data to the class list 
                    this->data.push_back(v[0]);
                }
            }
            // sets a bool if based on if the 
            if (!go) go = (v[0] == (std::string)"Perm:");
        }
    } else {
        // catching if the file did not open
        error(FLERR, ((std::string)"data file did not open, " + this->data_file).c_str());
    }

    // Close the file
    data_file.close();

}

/**
 * the main function for this class
*/
void ElmerToSpartaTemperature::make_sparta_surf_data() {
    // loading the boundary data
    this->load_boundary();

    // loading the data from the elmer case
    this->load_data();

    // opening the output file
    std::ofstream surf_data_file;
    surf_data_file.open(this->file_stem + ".surf.data");

    // catching if the file opened
    if (!(surf_data_file.is_open())) error(FLERR, "surf data file did not open");

    // setting the number of digits of a double to write to a file
    surf_data_file << std::fixed << std::setprecision(max_precision);

    // writing the header to the file
    surf_data_file << "ITEM: TIMESTEP\n";
    surf_data_file << "NAN\n";
    surf_data_file << "ITEM: NUMBER OF SURFS\n";
    surf_data_file << "NAN\n";// std::to_string(tris.size())//(len(tris)),
    surf_data_file << "ITEM: BOX BOUNDS NAN NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "-NAN NAN\n";
    surf_data_file << "ITEM: SURFS id c\n";

    // used later
    double avg;
    // std::cout << "averaging\n";
    for (int i = 0; i < this->boundary_data.size(); i++) {
        // computes the average temperature of the nodes that make up the surface element
        // this value is used to set the surface element temperature
        avg = 0;
        for (int j = 1; j < this->boundary_data[i].size(); j++) {
            // std::cout << this->data[std::stoi(boundary_data[i][j])-1] << " ";
            // gets the data point corresponding to node id and adds it to the rolling sum
            avg += std::stod(this->data[std::stoi(boundary_data[i][j])-1]);
        }
        // std::cout << "\n";
        // computing the average by dividing the sum by the number of points
        avg = avg/(this->boundary_data[i].size() - 1);

        // writing the id as well as the avg to the sparta data file
        surf_data_file << this->boundary_data[i][0] << " " << avg << "\n";
    }
    // closing the file
    surf_data_file.close();
}