import pylab, cplex, sys
from thermodynamics import R

class Stoichiometric_LP():
    
    def __init__(self, name='Stoichiometric_LP', log_file=sys.stderr):
        self.cpl = cplex.Cplex()
        self.cpl.set_problem_name(name)
        self.cpl.set_log_stream(log_file)
        self.cpl.set_results_stream(log_file)
        self.cpl.set_warning_stream(log_file)
        
        self.S = None
        self.weights = None
        self.compounds = None
        self.encountered_cids = set()
        self.reactions = None
        self.solution_index = 0
        self.flux_upper_bound = 100
        self.milp = False
        self.margin = False
        self.use_dG_f = False
        self.target_reaction = None

    def add_stoichiometric_constraints(self, weights, S, compounds, reactions, source, target):
        """
            S is a NCxNR matrix where the rows represent compounds and the columns represent reactions.
            b is the linear constraint vector, and is NCx1
            
            compounds is the list of compounds from KEGG (NC long)
            reactions is a list of pairs of RID and direction (NR long)
        """
        (self.Ncompounds, self.Nreactions) = S.shape
        self.S = S
        self.weights = weights
        self.compounds = compounds
        self.reactions = reactions
        
        # reaction fluxes are the continuous variables
        for r in range(len(self.reactions)):
            self.cpl.variables.add(names=[self.reactions[r].name], lb=[0], ub=[self.flux_upper_bound])
        
        # add a linear constraint on the fluxes for each compound (mass balance)
        for c in range(len(self.compounds)):
            constraint_name = "C%05d_mass_balance" % self.compounds[c].cid
            self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
            for r in pylab.find(self.S[c,:] != 0):
                self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name, self.S[c,r])

        if (source != None):
            self.add_flux("SOURCE", source, lb=1, ub=1)
        
        if (target != None):
            self.add_flux("TARGET", target, lb=-1, ub=-1)

    def add_flux(self, name, sparse, lb=1, ub=1):
        """
            Adds a constraint enforcing the given reaction to have a specific flux.
            The mass-balance constraint will be formulated without taking this flux into consideration.
            This is useful for formulating compound uptake fluxes or biomass reactions.
            
            Note that by changing lb and ub you can add flexibility to this flux.
        """
        self.cpl.variables.add(names=[name], lb=[lb], ub=[ub])
        for (cid, coeff) in sparse.iteritems():
            self.cpl.linear_constraints.set_coefficients("C%05d_mass_balance" % cid, name, coeff)        

    def export(self, fname):
        self.cpl.write(fname, filetype='lp')

    def add_flux_constraint(self, max_flux):
        self.cpl.linear_constraints.add(names=['flux'], senses='L', rhs=[max_flux])
        for (r, weight) in self.weights:
            self.cpl.linear_constraints.set_coefficients('flux', self.reactions[r].name, weight)

    def add_milp_variables(self):
        for r in range(len(self.reactions)):
            # add boolean indicator variables for each reaction, and minimize their sum
            self.cpl.variables.add(names=[self.reactions[r].name + "_gamma"], types='B')

            # add a constrain for each integer variable, so that they will be indicators of the flux variables
            # using v_i - M*gamma_i <= 0
            constraint_name = self.reactions[r].name + "_bound"
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", -self.flux_upper_bound)
        
        self.milp = True

    def add_reaction_num_constraint(self, max_num_reactions):
        if (self.milp == False):
            raise Exception("Cannot add reaction no. constraint without the MILP variables")
        self.cpl.linear_constraints.add(names=['num_reactions'], senses='L', rhs=[max_num_reactions])
        for r in range(len(self.reactions)):
            self.cpl.linear_constraints.set_coefficients('num_reactions', self.reactions[r].name + "_gamma", 1)
    
    def add_dGr_constraints(self, cid2dG0_f, T=300):
        if (self.milp == False):
            raise Exception("Cannot add thermodynamic constraints without the MILP variables")
        
        self.use_dG_f = True

        self.T = T
        for r in range(len(self.reactions)):
            for (cid, coeff) in self.reactions[r].sparse.iteritems():
                self.encountered_cids.add(cid)
        
        for cid in sorted(self.encountered_cids):
            self.cpl.variables.add(names=["C%05d_conc" % cid], lb=[-1e6], ub=[1e6])

        for r in range(len(self.reactions)):
            constraint_name = "%s_thermo" % self.reactions[r].name
            dG0_r = 0
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L')
            for (cid, coeff) in self.reactions[r].sparse.iteritems():
                self.cpl.linear_constraints.set_coefficients(constraint_name, "C%05d_conc" % cid, coeff)
                dG0_r += coeff * cid2dG0_f.get(cid, 0)
            self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1e6)
            self.cpl.linear_constraints.set_rhs(constraint_name, 1e6 - dG0_r/(R*T))

        # if the dG0_f is unknown, the concentration must remain unbound
        # therefore we leave only the known compounds in the encountered_cids set
        self.encountered_cids = self.encountered_cids.intersection(cid2dG0_f.keys())
        
    def add_specific_dGf_constraints(self, cid2bounds):
        """
            add constraints on dG_f for compounds that have a dG0_f and have a specific bound on the concentrations
        """
        bounded_cids = set()
        for cid in sorted(self.encountered_cids):
            (curr_c_min, curr_c_max) = cid2bounds.get(cid, (None, None))
            if (curr_c_min != None):
                self.cpl.variables.set_lower_bounds("C%05d_conc" % cid, pylab.log(curr_c_min))
                bounded_cids.add(cid)
            if (curr_c_max != None):
                self.cpl.variables.set_upper_bounds("C%05d_conc" % cid, pylab.log(curr_c_max))
                bounded_cids.add(cid)
        
        # remove the bounded CIDs from the encountered CID set, to prevent adding global or margin constraints on these compounds
        self.encountered_cids = self.encountered_cids.difference(bounded_cids)

    def add_global_dGf_constraints(self, global_c_range=(1e-6, 1e-2)):
        """
            add global constraints on dG_f for compounds that have a dG0_f (and don't have a specific bound already)
        """
        
        # add constraints on dG_f for compounds that have a dG0_f
        for cid in sorted(self.encountered_cids):
            self.cpl.variables.set_lower_bounds("C%05d_conc" % cid, pylab.log(global_c_range[0]))
            self.cpl.variables.set_upper_bounds("C%05d_conc" % cid, pylab.log(global_c_range[1]))
    
    def add_margin_dGf_constraints(self, c_mid=1e-4):
        """
            Sets the thermodynamic constraints on each of the reactions.
            Note that when using margin optimization, there is no incentive to minimize the number of reactions or the flux,
            and this can cause the emergence of futile cycles in the solutions.
            To avoid this, we add
        """
        if (self.milp == False):
            raise Exception("Cannot add thermodynamic constraints without the MILP variables")
        
        self.margin = True        
        self.cpl.variables.add(names=["pCr"], lb=[0], ub=[1e6])

        for cid in sorted(self.encountered_cids):
            self.cpl.linear_constraints.add(names=["C%05d_conc_minimum" % cid], senses='G', rhs=[pylab.log(c_mid)])
            self.cpl.linear_constraints.set_coefficients("C%05d_conc_minimum" % cid, "C%05d_conc" % cid, 1)
            self.cpl.linear_constraints.set_coefficients("C%05d_conc_minimum" % cid, "pCr", 1)

            self.cpl.linear_constraints.add(names=["C%05d_conc_maximum" % cid], senses='L', rhs=[pylab.log(c_mid)])
            self.cpl.linear_constraints.set_coefficients("C%05d_conc_maximum" % cid, "C%05d_conc" % cid, 1)
            self.cpl.linear_constraints.set_coefficients("C%05d_conc_maximum" % cid, "pCr", -1)

    def add_localized_dGf_constraints(self, cid2dG0_f, cid2bounds, c_range, T=300):
        self.T = T
        for r in range(len(self.reactions)):
            dG0_r = 0
            for (cid, coeff) in self.reactions[r].sparse.iteritems():
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
                    dG0_r += coeff * R*T*pylab.log(curr_c_max)
                else:
                    dG0_r += coeff * R*T*pylab.log(curr_c_min)
            
            if (dG0_r != None and dG0_r > 0):
                # this reaction is a localized bottleneck, add a constraint that its flux = 0
                constraint_name = self.reactions[r].name + "_irreversible"
                self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
                self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name, 1)
                
                # another option is to constrain gamma to be equal to 0, but it is equivalent, I think.
                #self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1)
        
    def ban_futile_cycles(self):
        """
            In order to avoid futile cycles, add a potential for each compound, and constrain the sum on each reaction to be negative
            note that we use S here (not the full sparse reaction as in the thermodynamic constraints, because we want to avoid
            futile cycles in the general sense (i.e. even if one ignores the co-factors).
        """
        for c in range(len(self.compounds)):
            self.cpl.variables.add(names=["C%05d_potential" % self.compounds[c].cid], lb=[-1e6], ub=[1e6])

        for r in range(len(self.reactions)):
            constraint_name = "%s_futilecycles" % self.reactions[r].name
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[1e6])
            self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1e6)
            for c in pylab.find(self.S[:,r] != 0):
                self.cpl.linear_constraints.set_coefficients(constraint_name, "C%05d_potential" % self.compounds[c].cid, self.S[c,r])

    def ban_current_solution(self):
        if (self.milp == False):
            raise Exception("Cannot ban a solution without the MILP variables")

        constraint_name = "solution_%03d" % self.solution_index
        self.solution_index += 1

        self.cpl.linear_constraints.add(names=[constraint_name], senses='L')
        N_active = 0
        for r in range(len(self.reactions)):
            if (self.gammas[r] > 0.5 and self.fluxes[r] > 1e-6):
                self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1)
                N_active += 1
        self.cpl.linear_constraints.set_rhs(constraint_name, N_active - 1)

    def set_objective(self):
        if (self.milp == False): # minimize the total flux of reactions (weighted)
            obj = [(self.reactions[r].name, weight) for (r, weight) in self.weights]
        elif (self.margin): # minimize the margin
            obj = [("pCr", 1)]
        else: # minimize the number of reactions (weighted)
            obj = [(self.reactions[r].name + "_gamma", weight) for (r, weight) in self.weights]

        self.cpl.objective.set_linear(obj)

    def solve(self, export_fname=None):
        self.cpl.solve()
        if (self.milp == False):
            if (self.cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.optimal):
                self.fluxes = self.cpl.solution.get_values([rn.name for rn in self.reactions])
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
                self.fluxes = self.cpl.solution.get_values([rn.name for rn in self.reactions])
                self.gammas = self.cpl.solution.get_values([rn.name + "_gamma" for rn in self.reactions])
                return True
            else:
                return False
    
    def get_total_flux(self):
        return sum(self.fluxes)
    
    def get_fluxes(self):
        fluxes = []
        for r in range(len(self.reactions)):
            if (self.gammas[r] > 0.5 and self.fluxes[r] > 1e-6):
                fluxes.append((r, self.fluxes[r]))
        return fluxes
    
    def get_conc(self):
        if (self.use_dG_f == False):
            return None
        cid2conc = {}
        for cid in sorted(self.encountered_cids):
            [conc] = self.cpl.solution.get_values(["C%05d_conc" % cid])
            cid2conc[cid] = conc
        return cid2conc

    def get_margin(self):
        if (not self.margin):
            return None
        else:
            # the objective is (log_c_max - log_c_min)
            # so its exponent will give c_max/c_min, i.e. the margin
            return pylab.exp(self.cpl.solution.get_objective_value()) 
        
