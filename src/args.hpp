#pragma once

#include <unordered_map>
#include <string>
#include <filesystem>
#include <memory>
#include <functional>

class Args;
#include "Solver.hpp"

using std::literals::operator""s;

class Args{
  protected:
    const char* command;
    std::unordered_map<std::string, std::function<char**(char**, char**)>> options{
      {"-h",               [this](char** begin, char** end)->char** {return this->print_help(begin, end);}},
      {"--help",           [this](char** begin, char** end)->char** {return this->print_help(begin, end);}},
      {"-m",               [this](char** begin, char** end)->char** {return this->set_model(begin, end);}},
      {"--model",          [this](char** begin, char** end)->char** {return this->set_model(begin, end);}},
      {"-o",               [this](char** begin, char** end)->char** {return this->set_output_file_name(begin, end);}},
      {"--output",         [this](char** begin, char** end)->char** {return this->set_output_file_name(begin, end);}},
      {"-s",               [this](char** begin, char** end)->char** {return this->set_solver_type(begin, end);}},
      {"--solver",         [this](char** begin, char** end)->char** {return this->set_solver_type(begin, end);}},
      {"--delta_t",        [this](char** begin, char** end)->char** {return this->set_delta_t(begin, end);}},
      {"--max_time",       [this](char** begin, char** end)->char** {return this->set_max_time(begin, end);}},
      {"--mesh_size",      [this](char** begin, char** end)->char** {return this->set_mesh_size(begin, end);}},
      {"--theta",          [this](char** begin, char** end)->char** {return this->set_theta(begin, end);}},
      {"--implicit_euler", [this](char** begin, char** end)->char** {return this->set_implicit_euler(begin, end);}},
      {"--explicit_euler", [this](char** begin, char** end)->char** {return this->set_explicit_euler(begin, end);}},
      {"--crank_nicolson", [this](char** begin, char** end)->char** {return this->set_crank_nicolson(begin, end);}}
    };

    std::string model = DEFAULT_MODEL;
    std::string output_file_name = DEFAULT_OUTPUT_FILE;
    std::unique_ptr<SolverAdapter> solver;
    
    double delta_t = DEFAULT_DELTA_T; //ms
    double max_time = DEFAULT_MAX_TIME; //ms;
    double mesh_size = DEFAULT_MESH_SIZE; //mm
    double theta = DEFAULT_THETA; 

    void print_help_and_exit(int exit_status);

    char** print_help(char** begin, char** end);
    char** set_model(char** begin, char** end);
    char** set_output_file_name(char** begin, char** end);
    char** set_solver_type(char** begin, char** end);
    char** set_delta_t(char** begin, char** end);
    char** set_max_time(char** begin, char** end);
    char** set_mesh_size(char** begin, char** end);
    char** set_theta(char** begin, char** end);
    char** set_implicit_euler(char** begin, char** end);
    char** set_explicit_euler(char** begin, char** end);
    char** set_crank_nicolson(char** begin, char** end);


  public:
    Args(int argc, char** argv);
    const std::string& get_model() const;
    std::string get_output_file_name() const;
    SolverAdapter* get_solver();
    double get_delta_t() const;
    double get_max_time() const;
    double get_mesh_size() const;
    std::string get_mesh_filename() const;
    double get_theta() const;

    friend std::ostream& operator<<(std::ostream& os, const Args& args);
};
