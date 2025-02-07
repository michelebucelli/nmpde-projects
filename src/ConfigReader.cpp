#include "../include/ConfigReader.hpp"
#include <iostream>
#include <fstream>
#include <string>
ConfigReader::ConfigReader()
{
    requiredVariables = {
        "mesh_2d_path",
        "mesh_3d_path",
        "degree_velocity",
        "degree_pressure",
        "T",
        "deltat",
        "Re"};

    pattern = std::regex(R"(^([\w_]+)\s*=\s*([^\s]+)\s*$)");

    if (!readConfigFile())
    {
        std::cerr << "Error: Failed to read config file." << '\n';
    }

    // Ensure all variables are set
    ensureVariablesSet();
}

auto ConfigReader::readConfigFile() -> bool
{
    std::ifstream file("../parameters.config");
    if (!file.is_open())
    {
        std::cerr << "Error: Config file not found." << '\n';
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue;
        }

        std::smatch match;
        if (std::regex_match(line, match, pattern))
        {
            std::string variableName = match[1];
            std::string variableValue = match[2];

            if (variableName == "mesh_2d_path")
            {
                mesh2DPath = variableValue;
            }
            else if (variableName == "mesh_3d_path")
            {
                mesh3DPath = variableValue;
            }
            else if (variableName == "degree_velocity")
            {
                degreeVelocity = std::stoi(variableValue);
            }
            else if (variableName == "degree_pressure")
            {
                degreePressure = std::stoi(variableValue);
            }
            else if (variableName == "T")
            {
                simulationPeriod = std::stod(variableValue);
            }
            else if (variableName == "deltat")
            {
                timeStep = std::stod(variableValue);
            }
            else if (variableName == "Re")
            {
                Re = std::stod(variableValue);
            }
            else
            {
                std::cerr << "Warning: Unknown variable '" << variableName << "' in config file." << '\n';
            }
        }
        else
        {
            std::cerr << "Error: Invalid format in line: " << line << '\n';
            return false;
        }
    }

    file.close();

    return true;
}

auto ConfigReader::ensureVariablesSet() -> void
{
    for (const auto &var : requiredVariables)
    {
        if (!checkVariableExists(var))
        {
            promptUserForVariable(var);
        }
    }
}

auto ConfigReader::promptUserForVariable(const std::string &variableName) -> void
{
    std::string input;
    while (true)
    {
        try
        {
            std::cout << "Please provide a value for " << variableName << ": ";
            std::getline(std::cin, input);

            if (variableName == "mesh_2d_path")
            {
                mesh2DPath = input;
                if (mesh2DPath.empty())
                    throw std::invalid_argument("Path cannot be empty.");
            }
            else if (variableName == "mesh_3d_path")
            {
                mesh3DPath = input;
                if (mesh3DPath.empty())
                    throw std::invalid_argument("Path cannot be empty.");
            }
            else if (variableName == "degree_velocity")
            {
                degreeVelocity = std::stoi(input);
                if (degreeVelocity <= 0)
                    throw std::invalid_argument("Degree must be a positive integer.");
            }
            else if (variableName == "degree_pressure")
            {
                degreePressure = std::stoi(input);
                if (degreePressure <= 0)
                    throw std::invalid_argument("Degree must be a positive integer.");
            }
            else if (variableName == "T")
            {
                simulationPeriod = std::stod(input);
                if (simulationPeriod <= 0)
                    throw std::invalid_argument("Simulation period must be positive.");
            }
            else if (variableName == "deltat")
            {
                timeStep = std::stod(input);
                if (timeStep <= 0)
                    throw std::invalid_argument("Time step must be positive.");
            }
            else if (variableName == "Re")
            {
                Re = std::stod(input);
                if (Re <= 0)
                    throw std::invalid_argument("Reynolds number must be positive.");
            }
            else
            {
                throw std::invalid_argument("Unknown variable: " + variableName);
            }

            break; // Exit loop if input is valid
        }
        catch (const std::exception &e)
        {
            std::cerr << "Invalid input: " << e.what() << " Please try again." << '\n';
        }
    }
}

auto ConfigReader::checkVariableExists(const std::string &variableName) const -> bool
{
    if (variableName == "mesh_2d_path")
    {
        return !mesh2DPath.empty();
    }
    else if (variableName == "mesh_3d_path")
    {
        return !mesh3DPath.empty();
    }
    else if (variableName == "degree_velocity")
    {
        return degreeVelocity > 0;
    }
    else if (variableName == "degree_pressure")
    {
        return degreePressure > 0;
    }
    else if (variableName == "T")
    {
        return simulationPeriod > 0.0;
    }
    else if (variableName == "deltat")
    {
        return timeStep > 0.0;
    }
    else if (variableName == "Re")
    {
        return Re > 0.0;
    }
    return false;
}

auto ConfigReader::getMesh2DPath() const -> std::filesystem::path
{
    return mesh2DPath;
}

auto ConfigReader::getMesh3DPath() const -> std::filesystem::path
{
    return mesh3DPath;
}

auto ConfigReader::getDegreeVelocity() const -> int
{
    return degreeVelocity;
}

auto ConfigReader::getDegreePressure() const -> int
{
    return degreePressure;
}

auto ConfigReader::getSimulationPeriod() const -> double
{
    return simulationPeriod;
}

auto ConfigReader::getTimeStep() const -> double
{
    return timeStep;
}

auto ConfigReader::getRe() const -> double
{
    return Re;
}
