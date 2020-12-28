# INSTALL

1) Install [Ariadne](https://github.com/ariadne-cps/ariadne)
2) Install [yaml-cpp](https://github.com/jbeder/yaml-cpp)
	```
	sudo apt install libyaml-cpp-dev
	```
3) Clone the repository
	```
	git clone https://github.com/kriato/bilateral_teleoperation_ariadne
	```

4) Build
	```
	mkdir build
	cd build
	cmake ..
	make -j
	```

# USAGE
To read the parameters from config.yml:
```
#define PROTOTYPE 1 \\ in main.cpp
```
To read the parameters from defines.hpp:
```
#define PROTOTYPE 0 \\ in main.cpp
```
**Remember to recompile after doing it!**