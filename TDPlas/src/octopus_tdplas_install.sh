# Change tdplas and octopus paths below
TDPLASDIR=/home/gabriel/my-wavet/WaveT/TDPlas/
OCTOPUSDIR=/home/gabriel/my-forked-octopus/octopus/
cd $OCTOPUSDIR
make clean
rm libtdplas.a
rm tdplas.o
rm tdplas.mod
cd ./src
make clean
rm libtdplas.a
rm tdplas.o
rm tdplas.mod
cd ./hamiltonian
rm libtdplas.a
rm tdplas.o
rm tdplas.mod
cd $TDPLASDIR
make clean
rm libtdplas.a
rm tdplas.o
rm tdplas.mod
make
cp $TDPLASDIR/src/tdplas.mod $OCTOPUSDIR/src/hamiltonian/.
cd $OCTOPUSDIR/src
make
cd $OCTOPUSDIR
./install
