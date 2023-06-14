#ifndef UTIL_H
#define UTIL_H

#include "toml.h"
#include "classes.h"

void format_name(std::string_view _s, std::string& _out);

void table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<std::string>& _var, std::string _sep = "=");

bool is_number(const std::string& s);

template<typename T>
void id_table_value_parser(std::string _caller, toml::node_t __tbl, std::vector<T>& _var, std::string _sep = "=");

#endif