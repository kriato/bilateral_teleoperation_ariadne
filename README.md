# REQUIREMENTS

1) Install [Ariadne](https://github.com/ariadne-cps/ariadne)
2) Install [yaml-cpp](https://github.com/jbeder/yaml-cpp)
	```
	sudo apt install libyaml-cpp-dev
	```

# BUILD
1) Clone the repository
	```
	git clone https://github.com/kriato/bilateral_teleoperation_ariadne
	```

2) Build
	```
	mkdir build
	cd build
	cmake ..
	make -j
	```
3) The executable will be generated in build/bin, along with plots folder