# Set the minimum version of CMake that's required
cmake_minimum_required(VERSION 3.12)

# Set the project name
project(TrigScint VERSION 2.1.0
                  DESCRIPTION "Module for all things trigger scintillator."
                  LANGUAGES CXX
)

setup_library(name TrigScint
              dependencies DARK::Framework DARK::Event DARK::DetDescr DARK::Tools
              python_install_path ${PYTHON_INSTALL_PREFIX}
)
