rm api/NeighborKMC.*
rm api/modules.rst
sphinx-apidoc -o ./api ../NeighborKMC
make clean
make html
make html
