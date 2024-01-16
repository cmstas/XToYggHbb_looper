#!/bin/bash

# Start from XToYggHbb_looper/.
rm tmp_create_package/ -rf
mkdir -p tmp_create_package
cd tmp_create_package

mkdir -p XToYggHbb_looper
cp ../NanoCORE ../data/ ../utils/ XToYggHbb_looper/. -r # Copy relevant folders
tar -cf - --exclude=temp_data* ../cpp | tar -xf - -C XToYggHbb_looper/. # Copy cpp folder without the plot folders

cd XToYggHbb_looper/NanoCORE
make -j8 # If there is a need to make clean, do it locally
cd ../.. # Go back to tmp_create_package/
cd XToYggHbb_looper/cpp
make clean; make -j4
cd ../.. # Go back to tmp_create_package/
tar -chJf package.tar.gz XToYggHbb_looper

mv package.tar.gz ../.
#rm -rf tmp_create_package # Leave the package be to able to easily check what is done
