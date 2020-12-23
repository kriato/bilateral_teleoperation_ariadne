# INSTALL

1) Install [Ariadne](https://github.com/ariadne-cps/ariadne)

2) Clone the repository
	```
	git clone https://github.com/kriato/bilateral_teleoperation_ariadne
	```

3) Clone yaml-cpp repository (used for faster prototyping)
	```
	mkdir deps
	cd deps
	git clone https://github.com/jbeder/yaml-cpp
	```
4) Build everything
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