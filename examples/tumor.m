clear all

addpath("../src/mfiles");
addpath("../src/cppfiles");

system("make -C ../src/cppfiles/");

################################################################

max_nagents = 1e5; cfa = 3.0; A = 2.0; 

spawn_prob0 = 1e-3; 
crit_lifetime = 30; convert_prob0 = 1e-3;

dt = 5e-3; ekin0 = 0.05; nthreads = 8;

# Fields 
ngrid = 40; lbox = 20; L0 = 1;
Db = 0.02; k0 = 0.001;

################################################################

a = agents(nthreads);
a.r(1,:) = [0,0,0]; 
a.v(1,:) = [0,0,0]; 

spawn_prob = spawn_prob0;
convert_prob = convert_prob0;

b = L0.*ones(ngrid, ngrid, ngrid);
h = lbox/ngrid;
slice_id = int64(ngrid/2);

m = 1; n=1; 
while a.nagents<max_nagents

	# Agents 
	a.pforce(cfa, A);
	ekin = a.integrate('p', dt, ekin0); 

	# Fields
	p_field = a.field('p', ngrid, lbox);	
	rhs = -k0*p_field.*b + Db*laplace(b, h); 
	b = b + dt*rhs;

	# Spawn 
	agent_field = fieldtoa(a.r(1:a.nagents, :), b, lbox); 
	spawn_prob = spawn_prob0*agent_field; 
	a.spawn(spawn_prob);

	# Convert
	agent_field = fieldtoa(a.r(1:a.nagents, :), b, lbox); 
	convert_prob = convert_prob0*(1-agent_field); 
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
		subplot(2,2,4); surf(b(:,:,slice_id));	axis([0 ngrid, 0, ngrid, 0, 1.1]);	
		pause(0.001);

		m++;
	
		printf("\r n=%d, %.1f  num. agents=%d Ekin=%.3f ", n, dt*n, a.nagents, ekin); 
	end

	n++;
end

printf("\n");



