#ifndef CONFIG_HANDLER_HPP
#define CONFIG_HANDLER_HPP

#include "./Preconditioners.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cctype>
#include <stdexcept>

/**
 * @brief Trims leading and trailing whitespace from a string.
 *
 * @param str The input string to trim.
 * @return The trimmed string.
 */
std::string trim(const std::string &str)
{
    size_t first = str.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t");
    return str.substr(first, (last - first + 1));
}

/**
 * @brief Converts a string to an InflowRegime enum value.
 *
 * @param str The string to convert ("Steady" or "Unsteady").
 * @return The corresponding InflowRegime enum value.
 * @throws std::runtime_error if the string is invalid.
 */
InflowRegime stringToInflowRegime(const std::string &str)
{
    if (str == "Steady")
        return InflowRegime::Steady;
    if (str == "Unsteady")
        return InflowRegime::Unsteady;
    throw std::runtime_error("Invalid InflowRegime value: " + str);
}

/**
 * @brief Converts a string to a Preconditioner enum value.
 *
 * @param str The string to convert (e.g., "SIMPLE", "YOSIDA").
 * @return The corresponding Preconditioner enum value.
 * @throws std::runtime_error if the string is invalid.
 */
Preconditioner stringToPreconditioner(const std::string &str)
{
    if (str == "SIMPLE")
        return Preconditioner::SIMPLE;
    if (str == "YOSIDA")
        return Preconditioner::YOSIDA;
    if (str == "APPROX_SIMPLE")
        return Preconditioner::APPROX_SIMPLE;
    if (str == "APPROX_YOSIDA")
        return Preconditioner::APPROX_YOSIDA;
    if (str == "IDENTITY")
        return Preconditioner::IDENTITY;
    if (str == "BLOCK_IDENTITY")
        return Preconditioner::BLOCK_IDENTITY;
    if (str == "PCD")
        return Preconditioner::PCD;
    if (str == "APPROX_PCD")
        return Preconditioner::APPROX_PCD;
    throw std::runtime_error("Invalid Preconditioner value: " + str);
}

/**
 * @brief Parses a configuration file and extracts simulation parameters.
 *
 * The configuration file should contain key-value pairs in the format:
 * @code
 * key = value
 * @endcode
 * Example:
 * @code
 * mesh_file_name = ../mesh/Navier_Stokes_2D_fine.msh
 * degree_velocity = 2
 * degree_pressure = 1
 * regime = Steady
 * peak_velocity = 1.5
 * preconditioner = SIMPLE
 * T = 8.0
 * dt = 0.005
 * @endcode
 *
 * @param filename The path to the configuration file.
 * @param mesh_file_name Output parameter for the mesh file path.
 * @param degree_velocity Output parameter for the velocity polynomial degree.
 * @param degree_pressure Output parameter for the pressure polynomial degree.
 * @param regime Output parameter for the inflow regime.
 * @param peak_velocity Output parameter for the peak inlet velocity.
 * @param preconditioner Output parameter for the preconditioner type.
 * @param T Output parameter for the final simulation time.
 * @param dt Output parameter for the time step size.
 * @throws std::runtime_error if the file cannot be opened or parsed.
 */
void parseConfigFile(const std::string &filename,
                     std::string &mesh_file_name,
                     unsigned int &degree_velocity,
                     unsigned int &degree_pressure,
                     InflowRegime &regime,
                     double &peak_velocity,
                     Preconditioner &preconditioner,
                     double &T,
                     double &dt)
{

    std::ifstream configFile(filename);
    if (!configFile.is_open())
    {
        throw std::runtime_error("Error opening config file: " + filename);
    }

    std::map<std::string, std::string> configMap;

    std::string line;
    while (std::getline(configFile, line))
    {
        line = trim(line);
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#')
            continue;

        // Find the delimiter
        size_t delimiterPos = line.find('=');
        if (delimiterPos == std::string::npos)
            continue;

        std::string key = trim(line.substr(0, delimiterPos));
        std::string value = trim(line.substr(delimiterPos + 1));

        configMap[key] = value;
    }

    // Extract values from the map with error checking
    try
    {
        mesh_file_name = configMap.at("mesh_file_name");
        degree_velocity = std::stoi(configMap.at("degree_velocity"));
        degree_pressure = std::stoi(configMap.at("degree_pressure"));
        regime = stringToInflowRegime(configMap.at("regime"));
        peak_velocity = std::stod(configMap.at("peak_velocity"));
        preconditioner = stringToPreconditioner(configMap.at("preconditioner"));
        T = std::stod(configMap.at("T"));
        dt = std::stod(configMap.at("dt"));
    }
    catch (const std::out_of_range &e)
    {
        throw std::runtime_error("Missing required parameter in config file: " + std::string(e.what()));
    }
}

#endif