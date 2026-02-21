# 1. Activate your conda environment
conda activate hep

# 2. Navigate to the project directory
cd ~/Desktop/cmsModel

# 3. Clean any old build and create a fresh build directory
rm -rf build && mkdir build && cd build

# 4. Run CMake with the correct paths to Geant4 and ROOT (via conda prefix)
cmake -DGeant4_DIR=$HOME/software/geant4/geant4-v11.4.0/lib/cmake/Geant4 \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/hep \
      ..

# 5. Build
make -j$(sysctl -n hw.ncpu)

# 6. Run (interactive mode)
./cmsmodel

# Or run in batch mode with a macro:
./cmsmodel run.mac

