/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

double multiv_prob_2(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 50;  // TODO: Set the number of particles
  is_initialized = true;
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);
  

  	for (int i = 0; i < num_particles; i++) {
	  Particle particle;
	  particle.id = i;
	  particle.x = dist_x(gen);
	  particle.y = dist_y(gen);
	  particle.theta = dist_theta(gen);
	  particle.weight = 1.0;	  
	  particles.push_back(particle);
	  weights.push_back(particle.weight);
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/std::normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  std::default_random_engine gen;
  
   for (int i = 0; i < num_particles; i++) {
     
		if (fabs(yaw_rate) < 0.0001) {
     	particles[i].x = particles[i].x + (velocity * delta_t )*sin(particles[i].theta);
     	particles[i].y = particles[i].y + (velocity* delta_t )* cos(particles[i].theta);
          } else 
        {
        particles[i].x = particles[i].x + (velocity/yaw_rate)* (sin(particles[i].theta + yaw_rate * delta_t ) - sin(particles[i].theta));
     	particles[i].y = particles[i].y + (velocity/yaw_rate)* (- cos(particles[i].theta + yaw_rate * delta_t ) + cos(particles[i].theta));
     	particles[i].theta = particles[i].theta + yaw_rate * delta_t   ;       
        }
 	 std::normal_distribution<double> dist_x(0, std_pos[0]);
 	 std::normal_distribution<double> dist_y(0, std_pos[1]);
 	 std::normal_distribution<double> dist_theta(0, std_pos[2]);
     
     particles[i].x += dist_x(gen);
     particles[i].y += dist_y(gen);
     particles[i].theta += dist_theta(gen);  
     
     
	}

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  	
	for (int i = 0; i < observations.size(); i++) {
		double closest_dist = 999999;
		int closest_landmark_id = -1;

		for (int j = 0; j < predicted.size(); j++) {
		  double current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
		  if (current_dist < closest_dist) {
		    closest_dist = current_dist;
		    closest_landmark_id = predicted[j].id;;
		  }
		}
		observations[i].id = closest_landmark_id;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_std::normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
   for (int i = 0; i < num_particles; i++) {
        double p_x = particles[i].x;
    	double p_y = particles[i].y;
    	double p_theta = particles[i].theta;
     	particles[i].weight = 1.0;
     	vector<LandmarkObs> observations_m; 
     	  for (int j = 0; j < observations.size(); j++) {
             double obs_x = observations[j].x ;
             double obs_y = observations[j].y ;
   			 LandmarkObs obj ;
             obj.id = observations[j].id;
   			 obj.x = p_x + obs_x * cos(p_theta) - obs_y * sin(p_theta);
             obj.y = p_y + obs_x * sin(p_theta) + obs_y * cos(p_theta);
             observations_m.push_back(obj);
 			 }
     	
    vector<LandmarkObs> landmarks;
    for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
      LandmarkObs ld_obj;
      
      if (dist (map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f, p_x,  p_y ) <= sensor_range) {
          ld_obj.id = map_landmarks.landmark_list[k].id_i;
          ld_obj.x = map_landmarks.landmark_list[k].x_f;
          ld_obj.y = map_landmarks.landmark_list[k].y_f;
         landmarks.push_back(ld_obj);
      }
    }
     
    dataAssociation(landmarks, observations_m);
     for (int l=0; l < observations_m.size(); l++) {
       
      	 for  (int m=0; m< landmarks.size(); m++) {
           
           
           double prob;
           if (observations_m[l].id == landmarks[m].id){
           prob = multiv_prob_2(std_landmark[0], std_landmark[1], observations_m[l].x, observations_m[l].y, landmarks[m].x, landmarks[m].y);
           particles[i].weight = particles[i].weight * prob;
           //std::cout << "found prob = " << prob << "\n";
           }
         }
       
     }
     
     weights[i] = particles[i].weight;
     
   }

  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  	vector<Particle> resamples;

	
	std::default_random_engine gen;
	
	
	std::uniform_int_distribution<int> dist_particle(0, num_particles - 1);
	
	int idx = dist_particle(gen);
	
	double beta = 0.0;
	
	 
	
	for (int i = 0; i < particles.size(); i++) {
		std::uniform_real_distribution<double> random_weight(0.0, 2.0 * *max_element(weights.begin(), weights.end()));
		beta = beta + random_weight(gen);

	  while (beta > weights[idx]) {
	    beta = beta -  weights[idx];
	    idx = (idx + 1) % num_particles;
	  }
	  resamples.push_back(particles[idx]);
	}
	particles = resamples;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

