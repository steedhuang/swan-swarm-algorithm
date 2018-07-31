import sys

def write(num):
    with open('main.cpp', 'w') as f:
        header = \
"""
//Authors: Ce MA, Jun Steed HUANG
//Date: 2018-5-25
//Many robot group a few obstacle avoidance based on Swarm Dynamical Window Approach

#include <stdio.h>
#include <unistd.h>
#include <time.h>
#include <libplayerc++/playerc++.h>
#include <vector>
#include <chrono>
#include <ctime>
#include <math.h>
#include <float.h>
//#include <values.h>
#define W_MIN -45 // Minimum of W
#define W_MAX 45 // Maximum of W
#define X_MIN 0.1 // Minimum of X
#define X_MAX 1   // Maximum of X
#define X_SAMPLE_GRAIN 0.1 // Sample size of X
#define W_SAMPLE_GRAIN 1 // Sample size of W
#define ROBOT_R 0.75 // the logic robot Radius 
#define OBSTACLE_R 0.25 // the logic obstcale Radius
#define ALPHA 0.2 // the weight of dist in evaluation
#define BETA 0.05 // the weight of heading in evaluation
#define ZETA 0.1  // the weight of velocity in evaluation
#define GAMMA 0.1 // the weight of swarm behavior in evaluation
#define DIST_PARA 100 // the predefined dist value for no obstcale situation
#define TIME_GRAIN 0.1  // time interval size
#define LASER_RANGE 4 // probe reach
"""
        f.write(header)
        f.write("#define SWARM_SIZE "+ str(num) + "// number of robots\n")
        temp = \
"""
struct location
{
     double X;
     double Y;
     double YAW;
}typedef robot_loc_t;
using namespace PlayerCc;
using namespace std;
int flag = 0; // flag for shut down
double time_window;
double degree2radians(double degree){
      return degree*M_PI/180.0;
}
double radians2degree(double radian){
      return radian*180.0/M_PI;
}
double heading_evaluation(double Yaw){
      return 180-fabs(Yaw);
}
double dist_evaluation(double traj_end_x, double traj_end_y, RangerProxy &laser, int loc)
{     double dist = 100;
      double angle = loc - 90;
      double ob_x = laser[loc]*cos(angle);
      double ob_y = laser[loc]*sin(angle);
      double temp_dist = sqrt((ob_x-traj_end_x)*(ob_x-traj_end_x)+(ob_y-traj_end_y)*(ob_y-traj_end_y)) - ROBOT_R - OBSTACLE_R;
      return temp_dist;
}
double velocity_evaluation(double v){
      return fabs(v);
}
vector<int> obstacle_detection(RangerProxy &laser, double *min_dist)
{
      vector<int> obstacle_location;
      double min_distance = LASER_RANGE;
      for(int i=0;i<laser.GetRangeCount();i++)
      {
            if(laser[i]<LASER_RANGE)
            {     if(min_distance>laser[i])
                  {
                        min_distance = laser[i];
                  }
                  obstacle_location.push_back(i);
            }
      }
      *min_dist = min_distance;
      return obstacle_location;
}
vector<robot_loc_t> trajectory_calculation(double X_speed, double W_speed, double *x_end, double *y_end, double *Yaw_status, double time_window)
{     double current_time = 0.0;
      double current_xx = 0.0;
      double current_yy = 0.0;
      double current_yaw = 0.0;
      vector<robot_loc_t> trajs;
      while(current_time<time_window)
      {     current_time = current_time + TIME_GRAIN;
            current_xx = current_xx + X_speed*TIME_GRAIN*cos(degree2radians(fabs(current_yaw)));
            current_yy = current_yy + X_speed*TIME_GRAIN*sin(degree2radians(fabs(current_yaw)));
            current_yaw = current_yaw + TIME_GRAIN*W_speed;
            robot_loc_t lo;
            lo.X = current_xx;
            lo.Y = current_yy;
            lo.YAW = current_yaw;
            trajs.push_back(lo);
      }
      *x_end = current_xx;
      *y_end = current_yy;
      *Yaw_status = *Yaw_status + current_yaw;
      return trajs;
}
/**
 * determine the distance away from all neighbors
 **/
double follow_neighbour(double x_end, double y_end,\
                        FiducialProxy &neighborFinder)
{     double total_dist_to_neighbor = 0.0;
      for(int i =0; i< neighborFinder.GetCount(); i++){
            double x_dist = neighborFinder.GetFiducialItem(i).pose.px - x_end;
            double y_dist = neighborFinder.GetFiducialItem(i).pose.py - y_end;
            total_dist_to_neighbor =total_dist_to_neighbor+ sqrt(x_dist*x_dist+y_dist*y_dist)/(x_dist*y_dist); // Swan distance
      }
      return total_dist_to_neighbor;

}
void DynamicWindowApproach(double *X_speed, double *W_speed, double time_window,\
                              RangerProxy &laser,vector<int> obstcale_locs,\
                              Position2dProxy &p2d,\
                              FiducialProxy &neighborFinder)
{// the sampling window size and settings
      int W_windowSize = (int)((W_MAX-W_MIN)/W_SAMPLE_GRAIN);
      int X_windowSize = (int)((X_MAX-X_MIN)/X_SAMPLE_GRAIN);
      vector<double> Ws(W_windowSize);
      vector<double> Xs(X_windowSize);
      vector<vector<double> > dists(W_windowSize);
      vector<vector<double> > headings(W_windowSize);
      vector<vector<double> > velocitys(W_windowSize);
      vector<vector<double> > swarm_behavior_eva(W_windowSize);
      vector<vector<double> > Evalu(W_windowSize);
      double optimal_eva = -DBL_MAX;
      double total_dist = 0.0;
      double total_heading = 0.0;
      double total_velocity = 0.0;
      double total_swarm_behavior_eva =0.0;
      double optimal_x_speed = *X_speed;
      double optimal_w_speed = *W_speed;

      for(int i=0;i<W_windowSize; i++){
            Ws[i]= W_MIN + i*W_SAMPLE_GRAIN;
      } 
      for(int j=0;j<X_windowSize; j++){
            Xs[j]= X_MIN + j*X_SAMPLE_GRAIN;
      }
      // evaluating all possible v and w
      for(int i=0;i<W_windowSize;i++){
            dists[i]= vector<double>(X_windowSize);
            headings[i] = vector<double>(X_windowSize);
            velocitys[i] = vector<double>(X_windowSize);
            swarm_behavior_eva[i] = vector<double>(X_windowSize);
            for(int j=0;j<X_windowSize;j++){
                  //first generate the trajectory
                  double xx_end = 0.0;
                  double yy_end = 0.0;
                  double yaw_end = radians2degree(p2d.GetYaw());
                  vector<robot_loc_t> trajectory;
                  trajectory = trajectory_calculation(Xs[j], Ws[i], &xx_end, &yy_end, &yaw_end, time_window);
                  double dist_ = DIST_PARA;
                  //distance calculation
                  for(int k=0;k<obstcale_locs.size();k++){
                        double dist_temp = dist_evaluation(xx_end, yy_end, laser, obstcale_locs[k]);
                        if(dist_>dist_temp){
                            dist_ =  dist_temp;
                        }
                  }
                  if(dist_>2*OBSTACLE_R){
                        dist_ = 2*OBSTACLE_R;
                  }
                  if(dist_>1*OBSTACLE_R){
                        dists[i][j] = dist_;
                        total_dist = total_dist + dist_;
                  }else{
                        dists[i][j] = -1;
                  }
                  //heading calculation
                  headings[i][j] = heading_evaluation(yaw_end);

                  total_heading = total_heading + headings[i][j];
                  //velocity calculation
                  velocitys[i][j] = velocity_evaluation(Xs[j]);
                  total_velocity = total_velocity + velocitys[i][j];
                  swarm_behavior_eva[i][j] = follow_neighbour(xx_end, yy_end, neighborFinder);
                  total_swarm_behavior_eva = total_swarm_behavior_eva + swarm_behavior_eva[i][j];
                  printf("swarm: %f\\n",swarm_behavior_eva[i][j] );
            }
      }

      for(int i=0;i<W_windowSize;i++){
            Evalu[i] = vector<double>(X_windowSize);
            for(int j=0;j<X_windowSize;j++){
                  // do normalization to reduce outlier effects
                  Evalu[i][j] = ALPHA*(dists[i][j]/total_dist)+\
                                    BETA*(headings[i][j]/total_heading)+\
                                    ZETA*(velocitys[i][j]/total_velocity)+\
                                    GAMMA*(swarm_behavior_eva[i][j]/total_swarm_behavior_eva);
                  if(Evalu[i][j]>optimal_eva&&dists[i][j]!=-1){
                        optimal_eva = Evalu[i][j];
                        optimal_x_speed = Xs[j];
                        optimal_w_speed = Ws[i];
                  }
            }
      }
      *X_speed = optimal_x_speed;
      *W_speed = optimal_w_speed;


}
int main(int argc, char *argv[])
{	
      /*need to do this line in c++ only*/
      using namespace PlayerCc;
"""
        f.write(temp)
        for i in range(num):
            f.write('      PlayerClient    robot' + str(i) + '("localhost", ' + str(7000 + i) + ');\n')
            f.write('      Position2dProxy p2dProxy_robot' + str(i) + '(&robot' + str(i) + ',0);\n')
            f.write('      RangerProxy      laserProxy_robot' + str(i) + '(&robot' + str(i) + ',0);\n')
            f.write('      FiducialProxy neighbor_finderProxy_robot' + str(i) + '(&robot' + str(i) + ',0);\n')
        f.write('PlayerClient plyclnts[SWARM_SIZE] ={')
        for i in range(num-1):
            f.write('robot' + str(i) + ', ')
        f.write('robot' + str(num - 1))
        f.write('};\n')
        f.write('Position2dProxy p2dProxys[SWARM_SIZE] = {')
        for i in range(num-1):
            f.write('p2dProxy_robot' + str(i) + ', ')
        f.write('p2dProxy_robot' + str(num - 1))        
        f.write('};\n')
        f.write('RangerProxy rngProxys[SWARM_SIZE] = {')
        for i in range(num-1):
            f.write('laserProxy_robot' + str(i) + ', ')
        f.write('laserProxy_robot' + str(num - 1))
        f.write('};\n')
        f.write('FiducialProxy fidProxys[SWARM_SIZE]={')
        for i in range(num-1):
            f.write('neighbor_finderProxy_robot' + str(i) + ', ')
        f.write('neighbor_finderProxy_robot' + str(num - 1))
        f.write('};\n')
        temp = \
"""
     // SimulationProxy SimProxys[SWARM_SIZE]={simProxy_robot0, simProxy_robot1, simProxy_robot2, simProxy_robot3, simProxy_robot4, simProxy_robot5, simProxy_robot6, simProxy_robot7, simProxy_robot8, simProxy_robot9, simProxy_robot10, simProxy_robot11};
      vector<double> X_speeds(SWARM_SIZE);
      vector<double> W_speeds(SWARM_SIZE);
      for(int i=0;i<SWARM_SIZE;i++){
            X_speeds[i] = 0.25;
            W_speeds[i] = (SWARM_SIZE/2-i)*0.0017; // roughly 3.1415926/180
            p2dProxys[i].SetMotorEnable(1);
            p2dProxys[i].RequestGeom();
            rngProxys[i].RequestGeom();
            fidProxys[i].RequestGeometry();
            rngProxys[i].RequestConfigure();

      }
      srand(time(NULL));
	double Avoid_dist = 5-1.25*2; //this is twice size of length of robot
      double current_min_dist;
      time_window = 3.25;// every time_window to update once
//======================================LOOP=====================================================================
      while(true)
      {	
	    for(int i = 0;i<SWARM_SIZE;i++){
            // read from the proxies
            // robot.Read();
            plyclnts[i].Read();
            vector<int> locs;
            // Determine obstcales for each robot
            locs = obstacle_detection(rngProxys[i], &current_min_dist);
            DynamicWindowApproach(&X_speeds[i], &W_speeds[i], time_window, rngProxys[i], locs, p2dProxys[i], fidProxys[i]);
            p2dProxys[i].SetSpeed(X_speeds[i], dtor(W_speeds[i]));
            // Stop when some at destination
            if ( p2dProxys[i].GetXPos()>6.18)   
		flag = flag+1;
            }
            // More than 61.8% 
            if (flag > SWARM_SIZE*0.618)
            break;
            sleep(time_window);
      }
}

"""
        f.write(temp)
if __name__ == '__main__':
    num = int(sys.argv[1])
    write(num)
