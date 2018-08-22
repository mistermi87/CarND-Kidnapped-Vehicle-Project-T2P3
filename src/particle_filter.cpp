/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles=500;

	//creating random distributions for x,y,theta
	random_device rd;
    mt19937 gen(rd());
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// pick random values from above distributions and create vector
	for(int i=0; i<num_particles;++i){
        Particle p;
        p.x=dist_x(gen);
        p.y=dist_y(gen);
        p.theta=dist_theta(gen);
        //keep orientation within +/-PI
        if(p.theta>M_PI) p.theta-=2*M_PI;
        if(p.theta<-M_PI) p.theta+=2*M_PI;
        p.weight=1.0;
        particles.push_back(p);
	}
	//cout<<"Initialized"<<endl;
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    //creating random noise distributions for x,y,z after forward step (with zero-mean)
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist_pos_x(0, std_pos[0]);
    normal_distribution<double> dist_pos_y(0, std_pos[1]);
    normal_distribution<double> dist_pos_theta(0, std_pos[2]);

    //cout<<"Prediction starting"<<endl;

    for(int i=0; i<num_particles; ++i){
        Particle p = particles[i];
        if (fabs(yaw_rate) < 0.00001) {
            p.x = p.x + velocity * delta_t * cos(p.theta) + dist_pos_x(gen);
            p.y = p.y + velocity * delta_t * sin(p.theta) + dist_pos_y(gen);
        }
        else{
            p.x = p.x + velocity/yaw_rate*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta))+dist_pos_x(gen);
            p.y = p.y + velocity/yaw_rate*(-cos(p.theta+yaw_rate*delta_t)+cos(p.theta))+dist_pos_y(gen);
            p.theta=p.theta+yaw_rate*delta_t+dist_pos_theta(gen);
        }
    particles[i]=p;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                std::vector<LandmarkObs> observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	//for each particle:
    //cout<<"updating weights"<<endl;
	double sig_x=std_landmark[0];
	double sig_y=std_landmark[1];

	for(int i=0; i<num_particles;++i){
        Particle p = particles[i];
        //reinitialize the weight
        p.weight=1.0;

        // predicted measurement:
        // for each landmark calculate the predicted measurement taking into account the position and angle of the particle
        std::vector<LandmarkObs> predicted;

        for(auto landmark : map_landmarks.landmark_list){
            //distances between particle and landmark in map coordinates --> already covered translation part here as well
            double distm_x_pl=landmark.x_f-p.x;
            double distm_y_pl=landmark.y_f-p.y;
            //only create predicted observations of the landmarks within the particles sensor range (incl. measurement uncertainty) (save computational time)
            //if(sqrt(pow(distm_x_pl+2*sig_x,2)+pow(distm_y_pl+2*sig_y,2))<=sensor_range){
            if((fabs(distm_x_pl)+sig_x)<=sensor_range && (fabs(distm_y_pl)+sig_y)<=sensor_range){
                LandmarkObs pre;
                pre.id=landmark.id_i;
                //transformation from map coordinates to vehicle coordinates (only rotation)
                pre.x=cos(p.theta)*distm_x_pl+sin(p.theta)*distm_y_pl;
                pre.y=-sin(p.theta)*distm_x_pl+cos(p.theta)*distm_y_pl;

                predicted.push_back(pre);
            }
        }


        //compare the observed measurement and assign an ID to it



        //a temporary prediction that is matching the observation of the particle --> for final transfer to id of observation and for weight calculation
        int found_p_id=99;
        double found_p_x=0.0;
        double found_p_y=0.0;


        for(int j=0; j<observations.size(); ++j){
            LandmarkObs obs= observations[j];

            double minimum_dist=numeric_limits<double>::max();
            double distance;

            for(int k=0; k<predicted.size(); ++k){
                LandmarkObs pre= predicted[k];
                distance=dist(obs.x,obs.y,pre.x,pre.y);
                if(distance<minimum_dist){
                    minimum_dist=distance;
                    found_p_id=pre.id;
                    found_p_x=pre.x;
                    found_p_y=pre.y;
                }
            }
            observations[j].id=found_p_id;

            //for each match the probability is directly computed and contributes to the final weight
            double multi_variate_prob=1/(2*M_PI*sig_x*sig_y) * exp(-( pow(obs.x-found_p_x,2)/(2*pow(sig_x,2)) + pow(obs.y-found_p_y,2)/(2*pow(sig_y,2)) ));
            p.weight=p.weight*multi_variate_prob;
        }

        //add weight of current particle to weights vector for easier resampling later
        weights.push_back(p.weight);

        particles[i]=p;

        //Delete "predicted" to not run into problems when "pushing back" or reinitializing (varying vector length depending on observable landmarks in sensor range of each particle)
        vector<LandmarkObs>().swap(predicted);
	}


}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    //new vector of particles
    vector<Particle> resampled_particles;

    //discrete distribution according to weights of particles
    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(weights.begin(), weights.end());

    // resampling to new vector of particles using the distribution "d" created above
    for (int i=0; i<num_particles; ++i){
        resampled_particles.push_back(particles[d(gen)]);
     }

    //replace old particles by new, resampled particles
    particles = resampled_particles;

    //vector weights is now emptied, in order to be filled again in the next step.
    vector<double>().swap(weights);
    //...same with the "resampled_particles"
    vector<Particle>().swap(resampled_particles);
}


string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
