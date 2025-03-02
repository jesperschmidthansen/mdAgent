classdef field < handle
	
	properties
		value;
		ngrids;
		lbox, dx;
	endproperties

	methods
		
		function this = field(ngrid, lbox, value=0)
			this.ngrids = ngrid;
			this.lbox = lbox; this.dx = lbox/ngrid;
			this.value = value*ones(ngrid, ngrid, ngrid);
		end

		function pvalue = posvalues(this, agent)
			pvalue = fieldtoa(agent.r(1:agent.nagents, :), this.value, this.lbox); 
		end
		
		function d2dx2 = laplace(this)
			d2dx2 = laplace(this.value, this.dx);
		end
	
	endmethods

endclassdef


