clear all

addpath("../src/mfiles");
addpath("../src/cppfiles");

system("make -C ../src/cppfiles/");

################################################################

max_nagents = 1e5; cfa = 3.0; A = 2.0; 

spawn_prob0 = 1e-3; 
crit_lifetime = 30; convert_prob0 = 1e-3;

dt = 5e-3; ekin0 = 0.05; nthreads = 8;

# Field 
ngrid = 40; lbox = 20; L0 = 1;
Db = 0.01; k0 = 0.005;

################################################################

a = agents(nthreads);
a.r(1,:) = [0,0,0]; 
a.v(1,:) = [0,0,0]; 

spawn_prob = spawn_prob0;
convert_prob = convert_prob0;

a.agentfield(ngrid, lbox);
rho_nutrient = field(ngrid, lbox, L0);

slice_id = int64(ngrid/2);

m = 1; n=1; 
while a.nagents<max_nagents

	# Agents 
	a.pforce(cfa, A);
	ekin = a.integrate('p', dt, ekin0); 

	# Fields
	a.getfield('p');	
	rhs = -k0*a.field.value.*rho_nutrient.value + Db*rho_nutrient.laplace(); 
	rho_nutrient.value = rho_nutrient.value + dt*rhs;

	# Spawn 
	avalues = rho_nutrient.posvalues(a);
	spawn_prob = spawn_prob0*avalues; 
	a.spawn(spawn_prob);

	# Convert
	avalues = rho_nutrient.posvalues(a);
	convert_prob = convert_prob0*(1-avalues); 
	a.convert("p->n", crit_lifetime, convert_prob);

	# Visualization
	if rem(n, 250)==0
		np(m) = length(find(a.types == 'p'));
		nn(m) = length(find(a.types == 'n' ));
		rsq = a.r(1:a.nagents,1).^2 + a.r(1:a.nagents,2).^2 + a.r(1:a.nagents,3).^2;
		radius(m) = sqrt(max(rsq));	
		t = linspace(0, dt*n, m);

		subplot(2,2,1); a.plot("pn", "br", 10, false, 3); 	
		subplot(2,2,3); a.plot("pn", "br", 10, true, 3); 
		subplot(2,2,2);	semilogy(t, np, '-o', t, nn+1, '-o'); 
		subplot(2,2,4); surf(rho_nutrient.value(:,:,slice_id));	axis([0 ngrid, 0, ngrid, 0, 1.1]);	
		pause(0.001);

		m++;
	
		printf("\r n=%d, %.1f  num. agents=%d Ekin=%.3f ", n, dt*n, a.nagents, ekin); 
	end

	n++;
end

printf("\n");



