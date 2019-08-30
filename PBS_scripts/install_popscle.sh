module load gcc/4.9.2
module load cmake
module load zlib

mkdir build
cd build

cmake --disable-lzma -DHTS_INCLUDE_DIRS=/projects/stitzel-lab/lawlon/Software/htslib/include -DHTS_LIBRARIES=/projects/stitzel-lab/lawlon/Software/htslib/lib/libhts.a ..
make
