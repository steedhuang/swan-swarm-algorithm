# This is a robot definition file
# Authors: Ce MA, Jun Steed Huang, 2018
# The swan shape is different on front/back/left/right to see orientation

define laser_component sensor
(
    # define the size of each transducer [xsize ysize zsize] in meters
    size [0.025 0.025 0.025]
    # define the range bounds [min max]
    range [0 5]
    # define the angular field of view (by default in degrees)
    fov 180
    # define the color that ranges are drawn in the gui
    color_rgba [ 1 0 0 0.5] 
    samples 180
)

define robot_laser ranger
(   
    laser_component(pose [0.75 0 0 0])
)

define neighbour_finder fiducial
(   
    size [0.1 0.1 0.1]
    range_min 0.0
    range_max 10.0
    range_max_id 10.0
    fov 360
    ignore_zloc 1
)

# draw the robot shape
define robot position
(   size [0.6 0.5 0.4]
    origin [0.125 0 0 0]
    block
    (
        points 9
        point[0] [0.25 0]
        point[1] [0.625 0.25]
        point[2] [1 0]
        point[3] [1.25 0.5]
        point[4] [1 1]
        point[5] [0.625 0.75]
        point[6] [0.25 1]
        point[7] [0 0.75]
        point[8] [0 0.25]
        z [0 1]
    )
    drive "diff"
    robot_laser(pose [0.75 0.5 0 0] )
    neighbour_finder( )
)

# The end.
