clear all

addpath("../src/mfiles");
addpath("../src/cppfiles");

system("make -C ../src/cppfiles/");

################################################################

max_nagents = 2e5; cfa = 3.0; A = 2.0; 
spawn_prob0 = 1e-3; 

z_wall = -2; cfs = 5; sw = 4; 

dt = 5e-3; ekin0 = 0.1;

################################################################

a = agents(nthreads);
a.r(1,:) = [0,0,0]; 
a.v(1,:) = [0,0,0]; 

n=1; 
while a.nagents<max_nagents

	# Agent dynamics
	a.pforce(cfa, A);
	a.sforce(z_wall, cfs, sw);

	ekin = a.integrate('p', dt, ekin0); 
	# Spawn 
	a.spawn(spawn_prob0);

	# Visualization
	if rem(n, 250)==0
		a.plot("p", "b"); 
	
		printf("\r n=%d, %.1f  num. agents=%d Ekin=%.3f ", n, dt*n, a.nagents, ekin); 
		fflush(stdout);		
	end

	n++;
end

printf("\n");



