# Magnetic Line Plotting

### Configuration

This project is configured using `Visual Studio`, using external libraries such as `Eigen`, `Yaml`, etc., these external libraries are managed by `vcpkg`

In my environment, the project should be built under `Release` and `x86` mode, and `vcpkg` should have the packages as follow:

> eigen3:x64-windows-static 
>
> fmt:x64-windows-static
>
> yaml-cpp:x64-windows-static

Every time, in the initialization process, calculate the magnetization process and calculate the field around the surface of each micro-robot, then based on the field generate the magnetic line’s start points. In the first step, move the 

### Run the project

To run the project, you should enter into `./MagneticLinePlotting/Binary/Release/` and run `Demo.exe` with the command line parameters.

- `-t 0` two balls and each ball contains three chains
- `-t 1` two balls and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 2` two uniform balls
- `-t 3` two peanuts and each peanut is formed by two balls, each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 4` seven balls and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 5` seven uniform balls
- `-t 6` same as `-t 1` but with different orientation
- `-t 7` three balls arranged like an equilateral triangle and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 8` four balls arranged like a square and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 9` five balls arranged like three equilateral triangles and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 10` six balls arranged like an big equilateral triangle and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 11` three balls arranged in a special shape and each ball contains three disks with the ratio of short axis and long axis `ba`
- `-t 12` seven balls arranged in a hexagon and each ball contains three chains
- `-t 13` seven uniform balls arranged in a hexagon
- `-t 14` one ball contains three disks
- `-t 15` like `-t 11` but another special shape
- `-t 16` like `-t 7` but uniform balls
- `-t 17` like `-t 8` but uniform balls
- `-t 18` like `-t 17` but another shape
- `-t 19` like `-t 9` but uniform balls
- `-t 20` like `-t 10` but uniform balls
- `-t 21` two balls and each ball contains single bundle



### Process

- If you just want the magnetic line drawn in the plane not all the space, please turn on `PLANARMOTION`
- If you want to draw the total field not only the field by the micro-robots, please turn on `DRAWALLFIELD`
- If you want to change the allowed space to draw the magnetic line, feel free to change the range in `bool collide2particle(const VectorDd &) const`
- In the `void printTorque()`, there is some options for you to choose which torque to print
- If you want to rotate the external field, please turn on `ROTATEEXTERNALFIELD`
- If you want to test the interaction just for one ball, please turn on `TESTONEBALL`
- If you want to make the lines denser, please add the limitation of the `cnt` in `void generateSampleParticle()`
- If the distribution has any symmetric, please turn on the corresponding manipulation after generate the test particles in `void generateSampleParticle()`
