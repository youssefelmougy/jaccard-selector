#!/bin/bash

# Setup Selector version and its libraries+modules
echo "-----------------------------------------------------"
echo "Setting up Selector version and its libraries+modules"
echo "-----------------------------------------------------"

source $PWD/scripts/Selector_setup.sh

echo "---------------------------------------------------------"
echo "Selector version and its libraries+modules setup complete"
echo "---------------------------------------------------------"

# Setup Selector version and its libraries+modules
echo "----------------------------------------------------------"
echo "Setting up GenomeAtScale version and its libraries+modules"
echo "----------------------------------------------------------"

source $PWD/scripts/GenomeAtScale_setup.sh

echo "--------------------------------------------------------------"
echo "GenomeAtScale version and its libraries+modules setup complete"
echo "--------------------------------------------------------------"

# All setups are complete
echo "-----------------------------------------------------------------------------------------"
echo "Selector and GenomeAtScale setup complete! Please see /scripts/run.sh to run experiments."
echo "-----------------------------------------------------------------------------------------"

