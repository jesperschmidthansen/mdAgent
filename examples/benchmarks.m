clear all

addpath("../inst/");
addpath("../src/");

system("make -C ../src/");

################################################################

max_num_agents = 1e4; 
interaction_cutoff = 3.0; 
interaction_param = 2.0; 

spawn_probability = 1e-2; 
convert_probability = 1e-3;

av_life_time = 1.0;

dt = 5e-3; 
ekin0 = 0.05; 

nthreads = 8;

################################################################

a = agents(nthreads);

# Initial configuration of one agent
a.r(1,:) = [0,0,0]; a.v(1,:) = [0,0,0]; 

tic();
profile on 
while a.nagents<max_num_agents

	# Calculate the force acting on agents 
	a.pforce(interaction_cutoff, interaction_param);

	# Integrate forward in time - here only proliferating agents 
	ekin = a.integrate('p', dt, ekin0); 

	# Spawn new proliferating agents
	a.spawn('p', spawn_probability);

	# Convert proliferating agents to necrotic agents
	a.convert("p->n", av_life_time, convert_probability);
end
profile off
toc()
profshow

