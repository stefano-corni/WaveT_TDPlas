# Change tdplas and octopus paths below
TDPLASDIR=/home/gabriel/my-wavet/WaveT/TDPlas/
OCTOPUSDIR=/home/gabriel/my-forked-octopus/octopus/
cd $OCTOPUSDIR
make clean
cd ./src
make clean
cd ./hamiltonian
rm tdplas.mod
cd $TDPLASDIR
make clean
make
cp $TDPLASDIR/src/tdplas.mod $OCTOPUSDIR/src/hamiltonian/.
cd $OCTOPUSDIR
make -j
make install
