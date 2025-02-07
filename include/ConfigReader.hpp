#ifndef CONFIG_READER_HPP
#define CONFIG_READER_HPP

#include <string>
#include <regex>
#include <unordered_set>
#include <filesystem>

/*
    ConfigReader Class
    ------------------
    * Loads simulation config (mesh paths, polynomial degrees, simulation period, time step, Reynolds number) using regex and filesystem.
    * Validates all required parameters and prompts user if any are missing.
    * Provides public getters for easy access to these configuration values.
    ------------------
*/


class ConfigReader
{
private:
    std::filesystem::path mesh2DPath;                           ///< Path to the 2D mesh
    std::filesystem::path mesh3DPath;                           ///< Path to the 3D mesh
    int degreeVelocity = 0;                                     ///< Degree of polynomial for velocity
    int degreePressure = 0;                                     ///< Degree of polynomial for pressure
    double simulationPeriod = 0.0;                              ///< Simulation period (T)
    double timeStep = 0.0;                                      ///< Time step (deltat)
    double Re = 0.0;                                            ///< Reynolds number
    std::regex pattern{};                                       ///< Regular expression for parsing config
    std::unordered_set<std::string> requiredVariables = {};     ///< Required variables

    auto checkVariableExists(const std::string &variableName) const     -> bool;
    auto ensureVariablesSet()                                           -> void;
    auto promptUserForVariable(const std::string& variableName)         -> void;

public:
    ConfigReader();

    auto getMesh2DPath() const        -> std::filesystem::path;
    auto getMesh3DPath() const        -> std::filesystem::path;
    auto getDegreeVelocity() const    -> int;
    auto getDegreePressure() const    -> int;
    auto getSimulationPeriod() const  -> double;
    auto getTimeStep() const          -> double;
    auto readConfigFile()             -> bool;
    auto getRe() const                -> double;
    };

#endif 