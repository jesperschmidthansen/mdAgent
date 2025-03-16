clear all

addpath("../src/mfiles");
addpath("../src/cppfiles");

system("make -C ../src/cppfiles/");

################################################################

max_nagents = 1e4; cfa = 3.0; A = 2.0; 
spawn_prob0 = 1e-2; convert_prob0 = 1e-3; life_time = 1.0;
dt = 5e-3; ekin0 = 0.05; nthreads = 8;

################################################################

a = agents(nthreads);
a.r(1,:) = [0,0,0]; 
a.v(1,:) = [0,0,0]; 

tic();
profile on 
while a.nagents<max_nagents
	# Agents 
	a.pforce(cfa, A);
	ekin = a.integrate('p', dt, ekin0); 

	# Spawn 
	a.spawn('p', spawn_prob0);

	# Convert
	a.convert("p->n", life_time, convert_prob0);
end
profile off
toc()
profshow

