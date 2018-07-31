import sys

def write_robots(num):
    with open('swan_robot.cfg', 'w') as f:
        init_driver = """
driver
(		
    name "stage"
    plugin "stageplugin"

    provides ["simulation:0"]

    # load the named file into the simulator
    worldfile "obstacle_world.world"	
)
        """
        f.write(init_driver)    
        f.write('\n')
        for i in range(num):
            f.write("# robot" + str(i) + '\n')
            f.write('driver\n')
            f.write("(\n")
            f.write('    name "stage"\n')
            f.write('	provides ["' + str(7000+i) + ':position2d:0" "' + str(7000+i) + ':ranger:0" "' + str(7000+i) + ':fiducial:0"]\n')
            f.write('	model "robot' + str(i) + '"\n')
            f.write(")\n\n")
        f.write('# The end.\n')


if __name__ == '__main__':
    num = int(sys.argv[1])
    write_robots(num)