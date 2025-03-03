classdef agents < handle
	
	properties
		nagents; max_nagents;
		r, v, F; # position, vel., force
		types; # types init to 'p'
		life_times; # life times	
		field; # Agent field
	endproperties

	methods
		
		function this = agents(nthreads=4, max_nagents=3e5)
			this.nagents = 1;
			this.max_nagents = max_nagents;

			this.r = zeros(max_nagents, 3);
			this.v = zeros(max_nagents, 3);
			this.life_times = zeros(max_nagents, 1);
			this.types = char(zeros(max_nagents, 1)); this.types(1) = 'p';

			setnthreads(nthreads);
		end

		function agentfield(this, ngrid, lbox);

			this.field = field(ngrid, lbox);

		end

		function pforce(this, cf, A)
			this.F = zeros(this.max_nagents, 3);
			this.F = pair_force(this.F, this.r, cf, A, this.nagents); 
		end

		function sforce(this, z0, cf, A)
			this.F = surface_force(this.F, this.r, z0, cf, A, this.nagents);
		end

		function ekin = integrate(this, atype, dt, ekin0)
			
			idx = find( this.types == atype );			
			ntype = length(idx);

			if ntype == 1
				ekin = 0.0;
			else 
				for k=1:3
					this.v(idx,k) += this.F(idx,k).*dt; 
					this.r(idx,k) += this.v(idx,k).*dt;
				end
				ekin = 0.5*sum(this.v(idx,1).^2 + this.v(idx,2).^2 + this.v(idx,3).^2);
				ekin = ekin/ntype;
				fact = sqrt(1.0 + 0.1*(ekin0/ekin - 1.0));
				for k=1:3
					this.v(idx,k) = fact.*this.v(idx,k); 
				end
			end

			this.life_times(1:this.nagents) += dt;
		end

		function spawn(this, critr)
			
			if length(critr)==1
				critr = ones(1, this.nagents)*critr;
			end
		
			idx = find( this.types=='p' ); 
			np = length(idx);
			random = rand(np, 1);
		
			for n=1:np 
				i =  idx(n);	
				if random(n) < critr(i) 
					na = this.nagents + 1;
					if na > this.max_nagents 
						error("Agent limit exceeded\n");
					end

					rnew = this.r(i,:);
					this.r(na, :) = rnew;
					
					vnew = rand(1,3) - 0.5; 
					this.v(na,:) = vnew;
					for k=1:3
						this.v(1:na,k) = this.v(1:na,k) - sum(this.v(1:na,k))./na;
					end	

					this.types(na) = 'p';
					this.nagents = na;
				end
			end

		end


		function convert(this, types, ltime, critval)
			
			if length(critval)==1 
				critval = ones(1, this.nagents)*critval;
			end
			if length(ltime)==1
				ltime = ones(1, this.nagents)*ltime;
			end

			idx = find( this.types==types(1) );
			np = length(idx);
			random = rand(np, 1);
	
			for n=1:np
				i = idx(n);
				if random(n) < critval(i) && this.life_times(i) > ltime(i)
					this.types(i) = types(end);
				end
			end

		end

		function getfield(this, type)
			idx = find( this.types==type );
			r = this.r(idx,:);   
			this.field.value = atofield(r, this.field.ngrids, this.field.lbox);
		end

		function plot(this, types, colors, lim=10, slice=false, ssize=5)
				
			if length(colors)!=length(types)
				error("Color and types are not of same length");
			end

			ntypes = length(types);
			for n=1:ntypes

				idx = find(this.types==types(n)); 	

				if slice
					m=1; idxt=[]; 
					for k=1:length(idx) 
						if this.r(idx(k),2)> 0 
					 		idxt(m)=idx(k); 
							m++;
						end 
					end
					idx = idxt; 
				end

				plot3(this.r(idx,1), this.r(idx,2), this.r(idx,3), ... 
						'o', 'markersize', ssize, 'markerfacecolor', colors(n));
				hold on
			end

			if slice
				plot3([-lim lim],[0 0], [-lim, -lim], 'k', 'linewidth', 3);
				plot3([-lim -lim],[0 0], [-lim, lim], 'k', 'linewidth', 3);
				plot3([lim lim],[0 0], [-lim, lim], 'k', 'linewidth', 3)
 				plot3([lim -lim],[0 0], [lim, lim], 'k', 'linewidth', 3)
			end

			hold off
			
			axis([-lim, lim, -lim, lim, -lim, lim]);				
			pause(0.001);
		end

		function save(this, filename, write_opt="w", max_agents=1e5)
	
			format = filename(end-2:end);
			if strcmp(format, "xyz")
				fptr = fopen(filename, write_opt);
			
				if this.nagents > max_agents 
					warning("Increasing variable max_agents- VMD will not support this when appending file");
					max_agents = this.nagents;
				end
	
				fprintf(fptr, "%d\n", max_agents);  
				fprintf(fptr, "Written by mdAgent\n");	
				for n=1:this.nagents
					fprintf(fptr, "%c %f %f %f \n", this.types(n), this.r(n,1), this.r(n,2), this.r(n,3));
				end 	
				for n=this.nagents+1:max_agents
					fprintf(fptr, "D 0.0 0.0 0.0\n");
				end	

				fclose(fptr);
			elseif strcmp(format, "mat")		
				r = this.r; v = this.v; types = this.types; 
				save(filename, "r", "v", "types");	
			else
				error("Format not supported");
			end
				
		end


	endmethods

endclassdef


