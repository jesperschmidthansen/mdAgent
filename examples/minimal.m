clear all

addpath("../src/mfiles");
addpath("../src/cppfiles");

system("make -C ../src/cppfiles/");

################################################################

max_nagents = 1e3; cfa = 3.0; A = 2.0; 
spawn_prob0 = 1e-2; 
dt = 5e-3; ekin0 = 0.05; nthreads = 8;

################################################################

a = agents(nthreads);
a.r(1,:) = [0,0,0]; 
a.v(1,:) = [0,0,0]; 

while this.nagents<max_nagents
	# Agents 
	a.pforce(cfa, A);
	ekin = a.integrate('p', dt, ekin0); 

	# Spawn 
	a.spawn('p', spawn_prob0);
end




