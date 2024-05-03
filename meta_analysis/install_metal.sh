git clone git@github.com:statgen/METAL.git

# Ensure homebrew is installed
brew install cmake

mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
