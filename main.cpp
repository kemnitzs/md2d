#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <fstream>
#include <string>

constexpr int domainsize_x = 50;
constexpr int domainsize_y = 50;
constexpr double time_step = 0.0001;

struct Particle {
  double x;
  double y;
  double vx;
  double vy;
};

struct Cell {
  std::vector<Particle> particles;
  void store_particle( Particle& p ){
    particles.push_back( p );
  }
  size_t get_size_in_memory() {
    return particles.capacity() * sizeof(Particle);
  }
};

struct Grid {
  std::array<std::array<Cell,domainsize_x>,domainsize_y> cells;

  void apply_periodic( Particle& p ) {
    p.x = fmod( p.x+domainsize_x*100, domainsize_x);
    p.y = fmod( p.y+domainsize_y*100, domainsize_y);
  }

  void store_particle( Particle& p ){
    apply_periodic( p );
    if ( !(p.x >= 0 && p.x < domainsize_x) || !(p.y >= 0 && p.y < domainsize_y)) {
      std::cout << "particle " << p.x << " " << p.y << std::endl;  
    }
    int addr_x = (int)p.x;
    int addr_y = (int)p.y;
    cells[addr_y][addr_x].store_particle( p );
  }

  size_t get_size_in_memory(){
    size_t sum = 0;
    for( auto& row : cells ){
      for( auto& cell : row ){
        sum += cell.get_size_in_memory();
      }
    }
    return sum;
  }

  std::array<double,2> get_direction(const Particle& p1, const Particle& p2){
    return { p2.x - p1.x, p2.y-p1.y};
  }

  double length(std::array<double,2> vect2d){
    return hypot( vect2d[0], vect2d[1] );
  }

  std::array<double,2> normalize( std::array<double,2> vect2d ){
    double l = length( vect2d );
    return { vect2d[0]/l, vect2d[1]/l };
  }


  double get_distance(const Particle& first_particle, const Particle& other_particle) {
    double dx = first_particle.x - other_particle.x;
    double dy = first_particle.y - other_particle.y;
    //auto [dx,dy] = get_direction( first_particle, other_particle);
    if(dx> domainsize_x/2){
	dx-=domainsize_x;
    }
    if(dx<-domainsize_x/2){
	dx+=domainsize_x;
    }
    if(dy>domainsize_y/2){
	dy-=domainsize_y;
    }
    if(dy<-domainsize_y/2){
	dy+=domainsize_y;
    }
    return hypot(dx,dy);
  }

  double get_force(double distance){
    // Lennard-Jones-Potential
    double alpha = 0.1;
    double beta = alpha * pow(0.5,1.0/1.6);
    double result =  12*pow(beta,12)/pow(distance,13) - 6*pow(beta, 6)/pow(distance,7);
    if (result > 10){
      result = 10;
    } else if (result < -10){
      result = -10;
    }
    return result; 
    //return 0;
  }

  void calculate_v( Particle& particle1, Particle& particle2 ){
    double distance = get_distance( particle1, particle2 );
    double force = get_force(distance);
    auto direction = get_direction(particle1,particle2); 
    auto norm_direction = normalize( direction );
    particle1.vx += - force * norm_direction[0] * time_step;  
    particle2.vy += - force * norm_direction[1] * time_step;
  }
  
  void calculate_v_with_all( Particle& particle1, int x1, int y1, int p1){
    for (int y = 0; y < cells.size(); ++y){
      for (int x = 0; x < cells[y].size(); ++x){
        for (int p = 0; p < cells[y][x].particles.size(); ++p){
	  if ( x1 == x && y1 == y && p1 == p ) continue;
	  Particle& particle2 = cells[y][x].particles[p];
	  calculate_v( particle1, particle2 ); 
	}
      }
    }
  }

  std::array<int,2> get_periodic_coordinate( int x, int y){
    return { (x + domainsize_x) % domainsize_x , (y + domainsize_y) % domainsize_y };
  }

  void calculate_v_with_area( Particle& particle1, int x1, int y1, int p1, int stencil_size ){
    int offset = stencil_size/2;
    int begin_x = x1-offset;
    int end_x = x1+offset;
    for (int y = y1-offset; y <= y1+offset; ++y){
      for (int x = x1-offset; x <= x1+offset; ++x){
        auto periodic_coordinate = get_periodic_coordinate( x, y );
	Cell& cell = cells[periodic_coordinate[1]][periodic_coordinate[0]];
	for (int p = 0; p < cell.particles.size(); ++p){
	  Particle& particle2 = cell.particles[p]; 
	  if ( periodic_coordinate[0] == x1 && periodic_coordinate[1] == y && p == p1 ) continue; 
	  calculate_v( particle1, particle2 ); 
	}
      }
    }
  }

  void calculate_background_force( Particle& p ){
    std::array<double, 2> direction = {
      center_of_force[0]-p.x,
      center_of_force[1]-p.y
    };
    auto len = length( direction );
    // Lennard-Jones-Potential
    double alpha = domainsize_x/3.0;
    double beta = alpha * pow(0.5,1.0/1.6);
    double force =  12*pow(beta,12)/pow(len,13) - 6*pow(beta, 6)/pow(len,7);
    if (force > 10){
      force =  10;
    } else if (force < -10){
      force = -10;
    }
    p.vx +=  -10 * force * direction[0] * time_step;
    p.vy +=  -10 * force * direction[1] * time_step;
    p.vx *=0.999;
    p.vy *=0.999;
  }

  void calculate_new_v(){
    //std::cout << "entering " << __FUNCTION__ << std::endl;
    #pragma omp parallel for collapse(2) schedule(static,1) 
    for (int y = 0; y < domainsize_y; ++y){
      for (int x = 0; x < domainsize_x; ++x){
        for (int p = 0; p < cells[y][x].particles.size(); ++p){
	  //calculate_v_with_all( cells[y][x].particles[p], x,y,p); 
	  calculate_v_with_area( cells[y][x].particles[p], x,y,p, 3 );
	  calculate_background_force( cells[y][x].particles[p] );
        } 
      }
    }
    //std::cout << "leaving " << __FUNCTION__ << std::endl;
  }

  void print( std::ostream& out ){
    for (int y = 0; y < domainsize_y; ++y){
      for (int x = 0; x < domainsize_x; ++x){
        out << cells[y][x].particles.size() << " ";
      }
      out << std::endl;
    }
  }

  void calculate_new_positions(){
    for (int y = 0; y < cells.size(); ++y){
      for (int x = 0; x < cells[y].size(); ++x){
        for (int p = 0; p < cells[y][x].particles.size(); ++p){
	  Particle& particle = cells[y][x].particles[p];
	  particle.x += particle.vx * time_step;
	  particle.y += particle.vy * time_step;
	}
      }
    }
  }

  void rebalance(){
    for (int y = 0; y < cells.size(); ++y){
      for (int x = 0; x < cells[y].size(); ++x){
        for (int p = 0; p < cells[y][x].particles.size(); ++p){
	  auto& particles = cells[y][x].particles;
	  Particle& particle = cells[y][x].particles[p];
	  if ( particle.x != x || particle.y != y ) {
	    std::swap( particle, particles.back() );
	    Particle particle_to_move = particles.back();
	    particles.pop_back();
	    store_particle( particle_to_move );
	  }
	}
      }
    }
  }

  std::array<double,2> center_of_force = { domainsize_x / 2 , domainsize_y /2 };

};

auto to_MB( size_t size_in_byte){
  return (double)size_in_byte / 1024 / 1024;
}

int main(int argc, char** argv){

  // generate 4 GB of particle data in the system
  //int number_of_particles = 4l*1024*1024*1024 / sizeof(Particle);
  int number_of_particles = 50*50*1; //*1024 / sizeof(Particle);

  Grid grid;

  std::cout << "number_of_particles " << number_of_particles << std::endl;
  std::cout << "size of the grid in memory " << to_MB(grid.get_size_in_memory()) << std::endl;

  std::mt19937 gen;
  std::uniform_real_distribution<> dis_x(0, domainsize_x);
  std::uniform_real_distribution<> dis_y(0, domainsize_y);
  std::uniform_real_distribution<> dis_vx(-1,1);
  std::uniform_real_distribution<> dis_vy(-1,1);

  for (int i = 0; i < number_of_particles; ++i){
    Particle p;
    p.x = dis_x(gen);
    p.y = dis_y(gen);
    p.vx = dis_vx(gen);
    p.vx = dis_vy(gen);
    grid.store_particle( p ); 
  }

  std::cout << "size of the grid in memory " << to_MB(grid.get_size_in_memory()) << std::endl;
  std::cout << "done generating particles" << std::endl;

  for (int i = 0; i < 24000; ++i){
    // calculate new v for all particles 
    grid.calculate_new_v();

    // calculate new pos for all particles
    grid.calculate_new_positions();
    grid.rebalance();
    if (i%2000==0){
      std::cout << i << std::endl;
      char name[100];
      sprintf(name, "%08d", i);
      std::ofstream out(std::string("mat_")+ name + ".dat");
      grid.print(out );
    }
  }
   
  return 0;
}



























