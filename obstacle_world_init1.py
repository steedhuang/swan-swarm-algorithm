import sys

def write_obstacle_world_cfg(col, row):
    with open('obstacle_world.world', 'w') as f:
        header = \
"""
# Authors: Dabu ZHANG, Jun Steed HUANG
# This is a world containing multiple robots and a set of random obstacle like head wind
# the goal is for majority robots move against wind by team formation
include "map.inc"
include "robot.inc"

# windowGUI adjustment
window
( 
   size [700.000 700.000] 
   scale 35
)

floorplan
(   
    bitmap "bitmaps/goal_map.png"
    size [15 15 1.5]
)

# define robot one by one, should be automated by a loop or a shell script
"""
        robot_info = \
"""    ranger_return -1
    obstacle_return 0
)
"""
        f.write(header)
        
        for col_j in range(col):
            for row_i in range(row):
                x = -3 - col_j/1.4 
                y =  - row / 2.8 + row_i/1.4 +0.3
                i = row_i + col_j * row 
                f.write('robot\n')
                f.write('(\n')
                f.write('    name "robot' + str(i) + '"\n')
                f.write('    pose [' + str(x) + ' ' + str(y) + ' 0 0]\n')
                f.write('    color "red"\n')
                fiducial_return = i + 1 # TODO fiducial代表什么？
                f.write('    fiducial_return ' + str(fiducial_return) + '\n')
                f.write(robot_info)
        end = \
"""

# this is the head wind drag 
define obstcale_b model
(   
    bitmap "bitmaps/circle.png"
    gui_nose 0
    gui_grid 0
    gui_move 0
    gui_outline 0
    gripper_return 0
    fiducial_return 0
    ranger_return 2
	size [0.25 0.25 1.5]
	color "orange"
)

# distributed with 0.618 ratio
obstcale_b(pose [-2.1 -6.21 0 0]) 
obstcale_b(pose [-1.7 -5.37 0 0]) 
obstcale_b(pose [-1 -4 0 0])  
obstcale_b(pose [0.1 -1.79 0 0])   
obstcale_b(pose [0.1 1.79 0 0])
obstcale_b(pose [-1 4 0 0])
obstcale_b(pose [-1.7 5.37 0 0])
obstcale_b(pose [-2.1 6.21 0 0])

# The end.

"""
        f.write(end)

if __name__ == '__main__':
    col = int(sys.argv[1]) # 列数
    row = int(sys.argv[2]) # 行数
    write_obstacle_world_cfg(col, row)
