#ifndef MATERIALS_H
#define MATERIALS_H

#include <string>
#include <vector>

class Material {
    public:
        std::string Name = "";
        std::string Heat_Conductivity = "";
        std::string Heat_Capacity = "";
        std::string Poisson_ratio = "";
        std::string Heat_expansion_Coefficient = "";
        std::string Youngs_modulus = "";
        std::string Sound_speed = "";
        std::string Density = "";
        
        // Material(std::string _name) {
        //     this->Name = _name;
        // }

        void to_string(std::vector<std::string>& out) {
            std::vector<std::string> data = {
                this->Name,
                this->Heat_Conductivity,
                this->Heat_Capacity,
                this->Poisson_ratio,
                this->Heat_expansion_Coefficient,
                this->Youngs_modulus,
                this->Sound_speed,
                this->Density
            };

            std::vector<std::string> labels = {
                "Name",
                "Heat Conductivity",
                "Heat Capacity",
                "Poisson ratio",
                "Heat expansion Coefficient",
                "Youngs modulus",
                "Sound speed",
                "Density"
            };

            for (int i = 0; i < data.size() - 1; i++) {
                out.push_back(labels[i]+" = "+data[i]);
            }
        }
};

class Data {
    public:
        Material Aluminum_Generic;
        void _set();

        Data() {
            this->_set();
        }
};

void Data::_set() {
    Aluminum_Generic.Name = ("Aluminum (generic)");
    Aluminum_Generic.Heat_Capacity = "897.0";
    Aluminum_Generic.Heat_expansion_Coefficient = "23.1e-6";
    Aluminum_Generic.Density = "2700.0";
    Aluminum_Generic.Heat_Conductivity = "237.0";
    Aluminum_Generic.Sound_speed = "5000.0";
    Aluminum_Generic.Youngs_modulus = "70.0e9";
    Aluminum_Generic.Poisson_ratio = "0.35";
}

#endif