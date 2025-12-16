OUTPUT = 2d_diffusive_convection 2d_finger_convection 2d_finger_convection \
		3d_finger_convection 2d_finger data yang2021jfm_case3_2d binary_fluid_convection \
		*flags.txt *args processinfo *.asc *.nc snapshots search* initial*\
		bisections energy* failures profiles tBisections yang2021jfm*


.PHONY: clean update build
update:
	sudo apt-get update
	sudo apt-get install libopenmpi-dev
	sudo apt install libeigen3-dev
	sudo apt-get install libfftw3-dev
	sudo apt-get install libfftw3-mpi-dev
	sudo apt install netcdf-bin libnetcdff-dev

clone:
	rm -rf channelflow
	git clone --depth=1 https://github.com/epfl-ecps/channelflow.git

cloneilc:
	rm -rf channelflow
	git clone --single-branch --branch module/ilc --depth=1 https://github.com/epfl-ecps/channelflow.git

build:
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	rm -rf build
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/user/local/ -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable " -DWITH_SHARED=ON -DWITH_HDF5CXX=OFF -DWITH_PYTHON=OFF;\
	make -j16
buildpython:
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	rm -rf build
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/user/local/ -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable " -DWITH_SHARED=ON -DWITH_HDF5CXX=OFF -DWITH_PYTHON=ON -DUSE_MPI=OFF;\
	make -j16
builduconn:
	git clone --depth=1 https://github.com/epfl-ecps/channelflow.git
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/gpfs/sharedfs1/admin/hpc2.0/apps/openmpi/5.0.5/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16
buildpsc:
	git clone --depth=1 https://github.com/epfl-ecps/channelflow.git
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/opt/packages/openmpi/gnu/5.0.3-gcc13.2.1-cpu/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16

buildnrel:
	git clone --depth=1 https://github.com/epfl-ecps/channelflow.git
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/opt/cray/pe/mpich/8.1.28/ofi/gnu/10.3/bin/mpicxx -DSET_RPATH=OFF -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16

clean:
	rm -rf build data $(OUTPUT)

allclean:
	rm -rf channelflow build $(OUTPUT)

compute_nu:
	mpiexec -n 16 ./build/modules/ddc/tools/ddc_nusselt -Ua -0.5 -Ub 0.5 -Ta 0.5 -Tb -0.5 -Sa 0.5 -Sb -0.5 t20000
	mpiexec -n 16 ./build/modules/ddc/tools/ddc_nusselt -scalar "salt" -Ua -0.5 -Ub 0.5 -Ta 0.5 -Tb -0.5 -Sa 0.5 -Sb -0.5 s20000

rotatingrbc_test:
	mpiexec -n 16 ./build/modules/ddc/tools/ddc_randomfield -Lx 1 -Lz 1 -ymin -0.5 -ymax 0.5 -Nx 64 -Ny 65 -Nz 64 iniU iniT
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "out_rrbc_test/" -Ra 2e5 -Pr 1 -Ek 1e-2 -Ua 0 -Ub 0 -Ta 0.5 -Tb -0.5 -dt 0.01 -dT 1 -s 10 -T 101 iniU iniT iniT

rotating_saltfinger_test:
	mpiexec -n 16 ./build/modules/ddc/tools/ddc_randomfield -Lx 1 -Lz 1 -ymin -0.5 -ymax 0.5 -Nx 64 -Ny 65 -Nz 64 iniU iniT iniS
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "out_rrbc_test/" -Ra 2e5 -Pr 1 -Ek 1e-2 -Le 100 -Rr  -Ua 0 -Ub 0 -Ta -0.5 -Tb 0.5 -Sa -0.5 -Sb 0.5 -dt 0.01 -dT 1 -s 10 -T 101 iniU iniT iniS