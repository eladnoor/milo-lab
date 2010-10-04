import pylab, cplex, sys, common

class YeastStoichiometricLP():
    
    def __init__(self, name='Stoichiometric_LP', milp=False, T=300, log_file=sys.stderr):
        self.cpl = cplex.Cplex()
        self.cpl.set_problem_name(name)
        self.cpl.set_log_stream(log_file)
        self.cpl.set_results_stream(log_file)
        self.cpl.set_warning_stream(log_file)
        
        self.compounds = []
        self.unbounded_compounds_with_dG0_f = set()
        self.reactions = []
        self.solution_index = 0
        self.flux_upper_bound = 100
        self.milp = milp
        self.T = T

    def export(self, fname):
        self.cpl.write(fname, filetype='lp')

    def add_flux_constraint_list(self, reactions):
        """
            adds a list of reactions to the problem.
            reactions is a list of pairs of (reaction ID, sparse reaction).
            
            *   both compound IDs and reaction IDs must be unique string representations.
            **  the 'sparse reaction' is a dictionary with keys that are compound IDs and values which are any numeric primitive.
            *** this program assumes all reactions are irreversible. If you want a reversible one, add the opposite reaction as well.
        """
        # reaction fluxes are the continuous variables
        for (rid, sparse) in reactions:
            self.add_flux_constraint(rid, sparse, 0, self.flux_upper_bound)

    def add_flux_constraint(self, rid, sparse, lb=1, ub=1):
        """
            adds a variable with a specified rid that has the value of the flux through a reaction
            and adds the effect that it has on the mass-balance formula for each of the participating
            compounds.
        """
        self.cpl.variables.add(names=[rid], lb=[lb], ub=[ub])
        for (cid, coeff) in sparse.iteritems():
            if (cid not in self.compounds):
                self.compounds.append(cid)
                self.cpl.linear_constraints.add(names=[cid + "_massbalance"], senses='E', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(cid + "_massbalance", rid, coeff)
        self.reactions.append((rid, sparse))
        
        if (self.milp):
            self.add_indicator_variable(rid)

    def add_indicator_variable(self, rid):
        """
            add boolean indicator variables for each reaction, usually used for thermodynamic
            constraints and/or minimizing the number of reactions rather than their total flux

            the new integer variable is constrained, so that they it will be 1 if the corresponding flux is > 0
            using v_i - M*gamma_i <= 0
        """
        self.cpl.variables.add(names=[rid + "_gamma"], types='B')
        self.cpl.linear_constraints.add(names=[rid + "_bound"], senses='L', rhs=[0])
        self.cpl.linear_constraints.set_coefficients(rid + "_bound", rid, 1)
        self.cpl.linear_constraints.set_coefficients(rid + "_bound", rid + "_gamma", -self.flux_upper_bound)

    def add_reaction_num_constraint(self, max_num_reactions):
        if (self.milp == False):
            raise Exception("Cannot add reaction no. constraint without the MILP variables")
        self.cpl.linear_constraints.add(names=['num_reactions'], senses='L', rhs=[max_num_reactions])
        for (rid, sparse) in self.reactions:
            self.cpl.linear_constraints.set_coefficients('num_reactions', rid + "_gamma", 1)
    
    def add_dGr_constraints(self, cid2dG0_f):
        if (self.milp == False):
            raise Exception("Cannot add thermodynamic constraints without the MILP variables")
        
        for cid in self.compounds:
            self.cpl.variables.add(names=[cid + "_concentration"], lb=[-1e6], ub=[1e6])

        for (rid, sparse) in self.reactions:
            dG0_r = 0
            self.cpl.linear_constraints.add(names=[rid + "_thermo"], senses='L')
            for (cid, coeff) in sparse.iteritems():
                self.cpl.linear_constraints.set_coefficients(rid + "_thermo", cid + "_concentration", coeff)
                dG0_r += coeff * cid2dG0_f.get(cid, 0)
            self.cpl.linear_constraints.set_coefficients(rid + "_thermo", rid + "_gamma", 1e6)
            self.cpl.linear_constraints.set_rhs(rid + "_thermo", 1e6 - dG0_r/(common.R*self.T))

        # if the dG0_f is unknown, the concentration must remain unbound
        # therefore we leave only the known compounds in the encountered_cids set
        self.unbounded_compounds_with_dG0_f = set(self.compounds).intersection(cid2dG0_f.keys())
        
    def add_specific_dGf_constraints(self, cid2bounds):
        """
            add constraints on dG_f for compounds that have a dG0_f and have a specific bound on the concentrations
            should always follow add_dGr_constraints() 
        """
        if (self.unbounded_compounds_with_dG0_f == None):
            raise Exception("Cannot add concentration constraints before calling add_dGr_constraints()")
        bounded_compounds = set()
        for cid in sorted(self.unbounded_compounds_with_dG0_f.intersection(cid2bounds.keys())):
            (curr_c_min, curr_c_max) = cid2bounds[cid]
            if (curr_c_min != None):
                self.cpl.variables.set_lower_bounds(cid + "_concentration", pylab.log(curr_c_min))
                bounded_compounds.add(cid)
            if (curr_c_max != None):
                self.cpl.variables.set_upper_bounds(cid + "_concentration", pylab.log(curr_c_max))
                bounded_compounds.add(cid)
        
        self.unbounded_compounds_with_dG0_f = self.unbounded_compounds_with_dG0_f.difference(bounded_compounds)

    def add_global_dGf_constraints(self, global_c_range=(1e-6, 1e-2)):
        """
            add global constraints on dG_f for compounds that have a dG0_f (and don't have a specific bound already)
            should always follow add_dGr_constraints() 
        """
        if (self.unbounded_compounds_with_dG0_f == None):
            raise Exception("Cannot add concentration constraints before calling add_dGr_constraints()")
        
        # add constraints on dG_f for compounds that have a dG0_f
        for cid in sorted(self.unbounded_compounds_with_dG0_f):
            self.cpl.variables.set_lower_bounds(cid + "_concentration", pylab.log(global_c_range[0]))
            self.cpl.variables.set_upper_bounds(cid + "_concentration", pylab.log(global_c_range[1]))
    
    def add_margin_dGf_constraints(self, c_mid=1e-4):
        """
            Sets the thermodynamic constraints on each of the reactions.
            Note that when using margin optimization, there is no incentive to minimize the number of reactions or the flux,
            and this can cause the emergence of futile cycles in the solutions.

            should always follow add_dGr_constraints() 
        """
        if (self.unbounded_compounds_with_dG0_f == None):
            raise Exception("Cannot add concentration constraints before calling add_dGr_constraints()")
        
        self.margin = True        
        self.cpl.variables.add(names=["pCr"], lb=[0], ub=[1e6])

        for cid in sorted(self.unbounded_compounds_with_dG0_f):
            self.cpl.linear_constraints.add(names=[cid + "_conc_minimum"], senses='G', rhs=[pylab.log(c_mid)])
            self.cpl.linear_constraints.set_coefficients(cid + "_conc_minimum", cid + "_concentration", 1)
            self.cpl.linear_constraints.set_coefficients(cid + "_conc_minimum", "pCr", 1)

            self.cpl.linear_constraints.add(names=[cid + "_conc_maximum"], senses='L', rhs=[pylab.log(c_mid)])
            self.cpl.linear_constraints.set_coefficients(cid + "_conc_maximum", cid + "_concentration", 1)
            self.cpl.linear_constraints.set_coefficients(cid + "_conc_maximum", "pCr", -1)

    def add_localized_dGf_constraints(self, cid2dG0_f, cid2bounds, c_range, T=300):
        self.T = T
        for (rid, sparse) in self.reactions:
            dG0_r = 0
            for (cid, coeff) in sparse.iteritems():
                if (cid in cid2dG0_f):
                    dG0_r += coeff * cid2dG0_f[cid]
                else:
                    dG0_r = None
                    break
                
                (curr_c_min, curr_c_max) = cid2bounds.get(cid, (None, None))
                if (curr_c_min == None):
                    curr_c_min = c_range[0]
                if (curr_c_max == None):
                    curr_c_max = c_range[1]

                if (coeff < 0):
                    dG0_r += coeff * common.R*T*pylab.log(curr_c_max)
                else:
                    dG0_r += coeff * common.R*T*pylab.log(curr_c_min)
            
            if (dG0_r != None and dG0_r > 0):
                # this reaction is a localized bottleneck, add a constraint that its flux = 0
                constraint_name = rid + "_irreversible"
                self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
                self.cpl.linear_constraints.set_coefficients(constraint_name, rid, 1)
                
                # another option is to constrain gamma to be equal to 0, but it is equivalent, I think.
                #self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1)
        
    def ban_futile_cycles(self):
        """
            In order to avoid futile cycles, add a potential for each compound, and constrain the sum on each reaction to be negative
            note that we use S here (not the full sparse reaction as in the thermodynamic constraints, because we want to avoid
            futile cycles in the general sense (i.e. even if one ignores the co-factors).
        """
        
        # TODO: this doesn't work yet, I don't know why
        
        for cid in self.compounds:
            self.cpl.variables.add(names=[cid + "_potential"], lb=[-1e6], ub=[1e6])

        for (rid, sparse) in self.reactions:
            self.cpl.linear_constraints.add(names=[rid + "_futilecycles"], senses='L', rhs=[1e6])
            self.cpl.linear_constraints.set_coefficients(rid + "_futilecycles", rid + "_gamma", 1e6)
            for (cid, coeff) in sparse.iteritems():
                self.cpl.linear_constraints.set_coefficients(rid + "_futilecycles", cid + "_potential", coeff)

    def ban_solution(self, rids):
        if (self.milp == False):
            raise Exception("Cannot ban a solution without the MILP variables")

        constraint_name = "solution_%03d" % self.solution_index
        self.solution_index += 1

        self.cpl.linear_constraints.add(names=[constraint_name], senses='L')
        for rid in rids:
            self.cpl.linear_constraints.set_coefficients(constraint_name, rid + "_gamma", 1)
        self.cpl.linear_constraints.set_rhs(constraint_name, len(rids) - 1)

    def set_flux_objective(self, obj):
        """
            obj must be a list of pairs of (rid, weight)
        """
        self.cpl.objective.set_linear(obj)
        self.cpl.objective.set_sense(self.cpl.objective.sense.maximize)

    def set_pCr_objective(self):
        obj = [("pCr", 1)] # minimize the margin
        self.cpl.objective.set_linear(obj)
        self.cpl.objective.set_sense(self.cpl.objective.sense.minimize)

    def solve(self, export_fname=None):
        self.cpl.solve()
        if (self.milp == False):
            if (self.cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.optimal):
                self.fluxes = self.cpl.solution.get_values([rid for (rid, sparse) in self.reactions])
                self.gammas = []
                for f in self.fluxes:
                    if (f > 1e-6):
                        self.gammas.append(1)
                    else:
                        self.gammas.append(0)
                return True
            else:
                return False
        else:
            if (self.cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.MIP_optimal):
                self.fluxes = self.cpl.solution.get_values([rid for (rid, sparse) in self.reactions])
                self.gammas = self.cpl.solution.get_values([rid + "_gamma" for (rid, sparse) in self.reactions])
                return True
            else:
                return False
    
    def get_fluxes(self):
        fluxes = []
        for r in range(len(self.reactions)):
            if (self.gammas[r] > 0.5 and self.fluxes[r] > 1e-6):
                fluxes.append((r, self.fluxes[r]))
        return fluxes
    
    def get_conc(self):
        cid2conc = {}
        for cid in sorted(self.encountered_cids):
            [conc] = self.cpl.solution.get_values([cid + "concentration"])
            cid2conc[cid] = conc
        return cid2conc

    def get_pCr(self):
        [pCr] = self.cpl.solution.get_values(["pCr"])
        return pCr