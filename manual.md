# multiRobotSimulation
This is a multi-Robot simulation project based on Player/Stage platform
# Preparation
## Environment Requirements
**OS:** Ubuntu 14.04(other version not testes, may crush due to lack of maintenance for Player/Stage)
**Player Version:** [3.1.0](https://github.com/playerproject/player/archive/release-3-1-0.tar.gz)
**Stage Version:** [4.1.1](https://github.com/rtv/Stage/archive/v4.1.1.tar.gz)
**Build Tools& Libs:**
* git
* cmake
* g++
* fltk1.1-dev
* libjpeg8-dev
* libpng12-dev
* libglu1-mesa-dev
* libltdl-dev
* libgtk2.0-dev
* build-essential
* libgnomecanvas2-dev
* gsl-bin
* libgsl0-dev
* libopencv-dev

**environment installation command**
```
sudo apt-get install git cmake g++ fltk1.1-dev libjpeg8-dev libpng12-dev libglu1-mesa-dev libltdl-dev libgtk2.0-dev build-essential libgnomecanvas2-dev gsl-bin libgsl0-dev libopencv-dev
```
## Installation
**Note:** Only way to install it is building it from source. As mentioned earlier, the project is not under maintenance. One must be careful with the environment and version
1. Setup Installation Directories
```
mkdir ~/PlayerStage
```
1. Setup Environmen Variables
```
export PLY=$HOME/StagePlayer/player 
export PATH=$PATH:$PLY/bin 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PLY/lib:$PLY/lib/x86_64-linux-gnu
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PLY/lib/x86_64-linux-gnu/pkgconfig

export STG=$HOME/StagePlayer/stage
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$STG/lib:$STG/lib64 
export PATH=$PATH:$STG/bin
source ~/.bashrc
```
1. Build Player
    1. Download Player source code to `~/PlayerStage`
    1. Change directory to Source file: `cd ~/PlayerStage/player-release-3-1-0`
    1. Make building work directory: `mkdir build`
    1. Goto building work directory: `cd build`
    1. Build and install
        * `cmake -DCMAKE_INSTALL_PREFIX=$PLY ..`
        * `make`
        * `make install`
1. Build Stage
    1. Download Stage source code to `~/PlayerStage`
    1. Change directory to Source file: `cd ~/PlayerStage/Stage-4.1.1`
    1. Make building work directory: `mkdir build`
    1. Goto building work directory: `cd build`
    1. Build and install
        * `cmake -DCMAKE_INSTALL_PREFIX=$STG ..`
        * `make`
        * `make install`
1. Test installation
    * `cd ~/PlayerStage/Stage-4.1.1/worlds`
    * `player worlds/simple.cfg`
    * `stage worlds/simple.world`
# Usage
## Player/Stage Platform
    The Following Picture illustrates the relation between control code, Player and Stage.
![simulation relation](./pics/ServerClient_sim.png)
    The Stage doing the acutal simulation computation, e.g. locations, speeds, collisions, etc. In player's perspective thru drivers, Stage is the real physical world. While Player is acting as proxy for control code and provides with useful interfaces. Stage and Player communicate thru sockets, so it can run on multiple devices with different physical and IP locations. 
    Developer can control the simulation thru configuring the simulation model and control code which is explained in next two section.
## Setup robot world model
* `.world` file
        this file contains the information about the robot world
* `.inc` file
        same with `.world` file except can be included by other world files
*  `.cfg` file
        this file will tell player which(*Dirver_name*) and where(*IP:Socket_Num*) driver to use, i.e. talks to which kind of devices
    * Setup a robot model
        1. Configure geometry shape
        1. Define physical properties:
            * visibility to sensors
            * collision mode
            * color
        1. Define sensors processed by robot
    * Setup world file
        1. Import bitmap of the map
        1. Setup physical properties
        1. Setup Robots and Obstacles positions
    * Setup cfg file
        1. Setup driver for Stage(or real robot) and sensors
        1. include `.world file`
    *Use scripts to create large amount of robots*
    *See details in [Player-Stage-Manual](http://player-stage-manual.readthedocs.io/en/latest/) chapter 2, 3, 4 and [official model guide](http://rtv.github.io/Stage/group__model.html)*
## Interfaces for control code
Control code talk to Player's proxies and thus talk to Stage, the simulation. There are two types of proxies, Simulation proxy(Top-down Knowledge) and Robot-perspective proxy(Distributed Knowledge).
* Simulation Proxy:
    **name:** `SimulationProxy`
    **Methdods:** 
        ```
        GetPose2d(char *item_name, double &x, double &y, double &yaw)  
        setPose2d(char *item_name, double x, double y, double yaw)  
        GetProperty(char *item_name, char *property, void *value, size_t value_len)  
        SetProperty(char *item_name, char *property, void *value, size_t value_len)
        ```
* Sensor proxies, methods omitted here
    * `Position2dProxy`
    * `RangerProxy`
    * `BlobfinderProxy`
    * `GripperProxy`
    * `FiducialProxy`
    *See details in [Player-Stage-Manual](http://player-stage-manual.readthedocs.io/en/latest/) Chapter 3&5 and [official proxy guide](http://playerstage.sourceforge.net/doc/Player-3.0.2/player/classPlayerCc_1_1ClientProxy.html)*
## Control Code
1. Setup player client
1. Setup proxies
1. Use info from proxies to make path & motion planning
## Build your code
Use makefile:  
    ```
    CXXFLAGS = `pkg-config --cflags playerc++` -g -std=gnu++11
    CFLAGS = `pkg-config --cflags playerc` -g
    LDLIBS = `pkg-config --libs playerc++`
    CC = g++
    ```
## Play with it!  
  1. First run cfg file: `player swan_robot.cfg`  
  1. Then run the control code `./main`
# MultiRobot Model Explanation
## Motion Model
Regarding robot as a vehicle, it has speed `V` meter/seconds moving forward and turning speed `W` rads/seconds. In simulation, the Player will set `V` and `W` every `dt` seconds. Therefore, we predict trajectory in this time window to avoid obstcales
## Robot Model
There are a number of robots. Each of them has two type of sensors, fiducial and laser
* Fiducial: Each robot carries with a beacon and a neighbour-finder that can locate neighbour with `10m`. Thus, we can achieve locating neighbours from robot perspective rather than get global/supervision information.
* Laser: Each robot carries a laser with 180 degree FOV at front of them.
## Obstcale Avoidance: Dynamic Window Approach
As we stated in Motion model, we have a time window `dt` to make decision about `V` and `W`.
1. Determine obstcale position relative to the robot
1. Sample out `V`s and `W`s in within sepcific range
1. Predict trajectory for each conbination of `V` and `W`
1. Evaluate each trajectory based on the end status of robot in that trajectory:
    * Dist: Summation of all Euclidean distances to obstcales
    * Heading: the orientation toward goal
    * Velocity: prefer higher velocity
    * Swarm behavior: diagonal behavior follow Finsler manifold pattern  
    Evaluation = ALPHA* Dist + BETA* Heading + ZETA* Velocity + GAMMA* Swarm Behavior
1. Set Each robot into selected `V` and `W`
## Trajectory Prediction
In `dt` time window, we calculate the final status of robot based on discrete time sampling.  
At each discrete time point, we compute the new position relative to original robot coordinating system((the one at the start of predicting trajectory)).   
```
current_X = current_X + V*time_interval*cos(current_Yaw)
current_Y = current_Y + V*time_interval*sin(current_Yaw)
current_YAW = current_YAW + W*time_interval
```
## Obstcales
Currently there are eight unmoving obstcales. Their positions can be changed at .world file or during the code running
The DWA need properly tuned parameters, and may not able to pass thru dense obstcale but take a way detour the obstcale area
## Result
![result_pic](./pics/result.gif)  
# Trouble Shooting
## Build problems
* Wrong operating system, versions or dependencies
    Solv: Use a clean OS setup again  
    Solv: Carefully try different version of source code, mayber newest one fixed bug.
* `Boost.Signals is no longer being maintained and is now deprecated. Please switch to Boost.Signals2`  
    Solv: [substitute "boost::TIME_UTC_" for "boost::TIME_UTC"](https://sourceforge.net/p/playerstage/mailman/message/32992365/)
* `CMakeFiles/playerprop.dir/playerprop.o: undefined reference to symbol`  
    Solv: go to the folder where you extracted the .tar file and search for "playerprop". Click on the playerprop.dir and look for a "link.txt" file. Add the line -lboost_system next to the first -lm tag you find.Repeat until the command "make" no longer encounters an error that makes it stop.[Src](https://forums.gentoo.org/viewtopic-t-1050322-start-0.html) 
## Robot world configuring problems
* `Segmentation fault while running playerstage`  
    Solv:[Array problems](https://sourceforge.net/p/playerstage/mailman/message/27504842/)
* Sensor problems  
    Solv: check the definition of sensor
## Control code problems
## Others
* `VMware: vmw_ioctl_command error Invalid argument`crush in Virtual Machine  
    Solv: `echo "export SVGA_VGPU10=0" >> ~/.bashrc`
# Future Work
* Add swarm learning behavior distributedly for each robot with respect to evaluation weight
* Scale up to thousands of robots and random obstacles with dynamical socket
* Migrate to ROS, since Player/Stage is not well-maintened and documented
# Resources
[Installation Blog](http://www.cnblogs.com/joseph-linux/p/5103403.html)  
[Online Player-Stage-Manual](http://player-stage-manual.readthedocs.io/en/latest/)  
[Github eaxmple codes](https://github.com/jennyhasahat/Player-Stage-Manual)  
[Player docs](http://playerstage.sourceforge.net/doc/Player-3.0.2/player/classPlayerCc_1_1ClientProxy.html)  
[Stage docs](http://rtv.github.io/Stage/group__model.html)  
[Mailing list in sourceforge](https://sourceforge.net/p/playerstage/mailman/)  
[The Dynamic Window Approach To Collision Avoidance](https://www.ri.cmu.edu/pub_files/pub1/fox_dieter_1997_1/fox_dieter_1997_1.pdf)  
[dynamic window based approach tomobile robot motion control in the presence of moving obstacles](http://ieeexplore.ieee.org/document/4209377/)  
[DWA Blog](http://blog.csdn.net/heyijia0327/article/details/44983551)