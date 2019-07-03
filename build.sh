#!/bin/bash

# if it does not exist, create a build folder
mkdir -p build
cd build

#initialize variable
unset CMAKE_OPTIONS;
CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=Debug"

for i in "$@"
do
case $i in
    -release)
    CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=Release"
    shift # past argument=value
    ;;
    *)
            # unknown option
    ;;
esac
done

# run cmake
echo $CMAKE_OPTIONS
echo "cmake $CMAKE_OPTIONS .."
cmake $CMAKE_OPTIONS ..

# run make and return to root directory
make
cd ..

