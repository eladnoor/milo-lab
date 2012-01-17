import pylab, cplex
from pygibbs.thermodynamic_constants import R
from pygibbs.thermodynamics import MissingCompoundFormationEnergy

class Stoichiometric_LP(object):
    
    def __init__(self, name='Stoichiometric_LP', log_file=None):
        self.cpl = cplex.Cplex()
        self.cpl.set_problem_name(name)
        self.cpl.set_log_stream(log_file)
        self.cpl.set_results_stream(log_file)
        self.cpl.set_warning_stream(log_file)
        
        self.S = None
        self.weights = None
        self.compounds = None
        self.cids_with_concentration = set()
        self.reactions = None
        self.solution_index = 0
        self.flux_upper_bound = 100
        self.milp = False
        self.pCr = False
        self.mtdf = False
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
        for r in xrange(len(self.reactions)):
            self.cpl.variables.add(names=[self.reactions[r].name], lb=[0], ub=[self.flux_upper_bound])
        
        # add a linear constraint on the fluxes for each compound (mass balance)
        for c, compound in enumerate(self.compounds):
            cid = compound.cid
                
            constraint_name = "C%05d_mass_balance" % cid
            
            self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
            for r in pylab.find(self.S[c,:]):
                self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name, self.S[c,r])
            
        if source:
            self.add_flux("SOURCE", source, lb=1, ub=1)
        
        if target:
            self.add_flux("TARGET", target, lb=-1, ub=-1)

    def add_flux(self, name, sparse, lb=1, ub=1):
        """
            Adds a constraint enforcing the given reaction to have a specific flux.
            The mass-balance constraint will be formulated without taking this flux into consideration.
            This is useful for formulating compound uptake fluxes or biomass reactions.
            
            Note that by changing lb and ub you can add flexibility to this flux.
        """
        self.cpl.variables.add(names=[name], lb=[lb], ub=[ub])
        for cid, coeff in sparse.iteritems():
            try:
                self.cpl.linear_constraints.set_coefficients("C%05d_mass_balance" % cid, name, coeff)
            except cplex.exceptions.CplexSolverError: #@UndefinedVariable
                print 'CID:', cid

    def export(self, fname):
        self.cpl.write(fname, filetype='lp')

    def add_flux_constraint(self, max_flux):
        self.cpl.linear_constraints.add(names=['flux'], senses='L', rhs=[max_flux])
        for r, weight in self.weights:
            self.cpl.linear_constraints.set_coefficients('flux', self.reactions[r].name, weight)

    def add_milp_variables(self):
        for r in self.reactions:
            # add boolean indicator variables for each reaction, and minimize their sum
            self.cpl.variables.add(names=[r.name + "_gamma"], types='B')

            # add a constrain for each integer variable, so that they will be indicators of the flux variables
            # using v_i - M*gamma_i <= 0
            constraint_name = r.name + "_bound"
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L', rhs=[0])
            self.cpl.linear_constraints.set_coefficients(constraint_name, r.name, 1)
            self.cpl.linear_constraints.set_coefficients(constraint_name, 
                r.name + "_gamma", -self.flux_upper_bound)
        
        self.milp = True

    def add_reaction_num_constraint(self, max_num_reactions):
        if not self.milp:
            raise Exception("Cannot add reaction no. constraint without the MILP variables")
        self.cpl.linear_constraints.add(names=['num_reactions'], senses='L', rhs=[max_num_reactions])
        for r in self.reactions:
            self.cpl.linear_constraints.set_coefficients('num_reactions', r.name + "_gamma", 1)
    
    def add_dGr_constraints(self, thermodynamics, pCr=False, MTDF=False, maximal_dG=0.0):
        """
            Create concentration variables for each CID in the database (at least in one reaction).
            If this compound doesn't have a dG0_f, its concentration will not be constrained.
            If it has a dG0_f, then either it has a specific concentration range (most co-factors do),
            or not - which is the common case. In this case, we either apply global constraints (e.g. 1uM - 10mM),
            or we use the pCr method to minimize the range needed for allowing the pathway thermodynamically.
        """
        
        if not self.milp:
            raise Exception("Cannot add thermodynamic constraints without the MILP variables")
        
        if pCr and MTDF:
            raise Exception("Cannot optimize both the pCr and the MTDF")
        
        self.use_dG_f = True

        if pCr:
            self.pCr = True        
            self.cpl.variables.add(names=["pCr"], lb=[0], ub=[1e6])
        
        if MTDF:
            self.mtdf = True
            self.cpl.variables.add(names=["mtdf"], lb=[-1e6], ub=[1e6])
        
        for r in self.reactions:
            self.cids_with_concentration = self.cids_with_concentration.union(r.sparse.keys())
        
        cids_with_dG0 = thermodynamics.get_all_cids()
        for cid in self.cids_with_concentration:
            self.cpl.variables.add(names=["C%05d_conc" % cid], lb=[-1e6], ub=[1e6])
            if cid not in cids_with_dG0:
                # If the formation energy of this compound cannot be determined,
                # it's concentration cannot be constrained. The value of C%05d_conc 
                # will actually represent the dG_f (both dG0_f and the concentration)
                # and can be used in the thermodynamic constraint of reactions
                continue
                
            (c_min, c_max) = thermodynamics.cid_to_bounds(cid, use_default=False)
                
            if c_min: # override any existing constraint with the specific one for this co-factor
                self.cpl.variables.set_lower_bounds("C%05d_conc" % cid, pylab.log(c_min))
            elif pCr: # use the pCr variable to define the constraint
                self.cpl.linear_constraints.add(names=["C%05d_conc_minimum" % cid],
                    senses='G', rhs=[pylab.log(thermodynamics.c_mid)])
                self.cpl.linear_constraints.set_coefficients(
                    "C%05d_conc_minimum" % cid, "C%05d_conc" % cid, 1)
                self.cpl.linear_constraints.set_coefficients(
                    "C%05d_conc_minimum" % cid, "pCr", 1)
            else: # otherwise, use the global concentration bounds
                self.cpl.variables.set_lower_bounds("C%05d_conc" % cid, 
                    pylab.log(thermodynamics.c_range[0]))
            
            if c_max: # override any existing constraint with the specific one for this co-factor
                self.cpl.variables.set_upper_bounds("C%05d_conc" % cid, pylab.log(c_max))
            elif pCr: # use the pCr variable to define the constraint
                self.cpl.linear_constraints.add(names=["C%05d_conc_maximum" % cid],
                    senses='L', rhs=[pylab.log(thermodynamics.c_mid)])
                self.cpl.linear_constraints.set_coefficients(
                    "C%05d_conc_maximum" % cid, "C%05d_conc" % cid, 1)
                self.cpl.linear_constraints.set_coefficients(
                    "C%05d_conc_maximum" % cid, "pCr", -1)                
            else: # otherwise, use the global concentration bounds
                self.cpl.variables.set_upper_bounds(
                    "C%05d_conc" % cid, pylab.log(thermodynamics.c_range[1]))

        for r in self.reactions:
            constraint_name = "%s_thermo" % r.name
            dG0_r = 0
            self.cpl.linear_constraints.add(names=[constraint_name], senses='L')
            for cid, coeff in r.sparse.iteritems():
                self.cpl.linear_constraints.set_coefficients(constraint_name, "C%05d_conc" % cid, coeff)
                try:
                    dG0_r += coeff * thermodynamics.GetTransformedFormationEnergies(cid)
                except MissingCompoundFormationEnergy:
                    # if this CID is not in cid2dG0, it means its formation energy is 
                    # part of its concentration variable, and therefore it doesn't contribute to dG0_r
                    continue
            
            gamma_factor = 1e6
            self.cpl.linear_constraints.set_coefficients(constraint_name, r.name + "_gamma", gamma_factor)
            if self.mtdf:
                self.cpl.linear_constraints.set_coefficients(constraint_name, "mtdf", -1)
            
            rhs = gamma_factor - (dG0_r - maximal_dG)/(R*thermodynamics.T) 
            self.cpl.linear_constraints.set_rhs(constraint_name, rhs) 

    def add_localized_dGf_constraints(self, thermodynamics):
        for r in self.reactions:
            dG0_r = 0
            for cid, coeff in r.sparse.iteritems():
                try:
                    dG0_r += coeff * thermodynamics.GetTransformedFormationEnergies(cid)
                except MissingCompoundFormationEnergy:
                    dG0_r = None
                    break
                
                (curr_c_min, curr_c_max) = thermodynamics.cid_to_bounds(cid)

                if (coeff < 0):
                    dG0_r += coeff * R * thermodynamics.T * pylab.log(curr_c_max)
                else:
                    dG0_r += coeff * R * thermodynamics.T * pylab.log(curr_c_min)
            
            if (dG0_r != None and dG0_r > 0):
                # this reaction is a localized bottleneck, add a constraint that its flux = 0
                constraint_name = r.name + "_irreversible"
                self.cpl.linear_constraints.add(names=[constraint_name], senses='E', rhs=[0])
                self.cpl.linear_constraints.set_coefficients(constraint_name, r.name, 1)
                
                # another option is to constrain gamma to be equal to 0, but it is equivalent, I think.
                #self.cpl.linear_constraints.set_coefficients(constraint_name, self.reactions[r].name + "_gamma", 1)
        
    def ban_futile_cycles(self):
        """
            In order to avoid futile cycles, add a potential for each compound, and constrain the sum on each reaction to be negative
            note that we use S here (not the full sparse reaction as in the thermodynamic constraints, because we want to avoid
            futile cycles in the general sense (i.e. even if one ignores the co-factors).
        """
        for c in self.compounds:
            self.cpl.variables.add(names=["C%05d_potential" % c.cid], lb=[-1e6], ub=[1e6])

        for r, reaction in enumerate(self.reactions):
            constraint_name = "%s_futilecycles" % reaction.name
            
            self.cpl.linear_constraints.add(names=[constraint_name], 
                senses='L', rhs=[1e6])
            self.cpl.linear_constraints.set_coefficients(constraint_name, 
                reaction.name + "_gamma", 1e6)
            
            for c in pylab.find(self.S[:,r] != 0):
                self.cpl.linear_constraints.set_coefficients(constraint_name, 
                    "C%05d_potential" % self.compounds[c].cid, self.S[c,r])

    def set_objective(self):
        if (self.milp == False): # minimize the total flux of reactions (weighted)
            obj = [(self.reactions[r].name, weight) for (r, weight) in self.weights]
        elif (self.pCr): # minimize the pCr
            obj = [("pCr", 1)]
        elif (self.mtdf):
            obj = [("mtdf", 1)]
        else: # minimize the number of reactions (weighted)
            obj = [(self.reactions[r].name + "_gamma", weight) for (r, weight) in self.weights]

        self.cpl.objective.set_linear(obj)

    def solve(self, export_fname=None):
        try:
            self.cpl.solve()
        except cplex.exceptions.CplexSolverError: #@UndefinedVariable
            return False
            
        if not self.milp:
            if (self.cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.optimal):
                self.fluxes = self.cpl.solution.get_values([rn.name for rn in self.reactions])
                self.gammas = []
                for f in self.fluxes:
                    if (f > 1e-6):
                        self.gammas.append(1)
                    else:
                        self.gammas.append(0)
                self.log_concentrations = None
                return True
            else:
                return False
        else:
            if (self.cpl.solution.get_status() == cplex.callbacks.SolveCallback.status.MIP_optimal):
                self.fluxes = self.cpl.solution.get_values([rn.name for rn in self.reactions])
                self.gammas = self.cpl.solution.get_values([rn.name + "_gamma" for rn in self.reactions])
                self.log_concentrations = []
                if self.use_dG_f:
                    for c in self.compounds:
                        if (c.cid in self.cids_with_concentration):
                            self.log_concentrations += self.cpl.solution.get_values(["C%05d_conc" % c.cid])
                        else:
                            self.log_concentrations += [pylab.NaN]
                return True
            else:
                return False

    def ban_current_solution(self):
        if not self.milp:
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
    
    def get_total_flux(self):
        return sum(self.fluxes)
    
    def get_fluxes(self):
        reactions = []
        fluxes = []
        for r in xrange(len(self.reactions)):
            if self.gammas[r] > 0.5 and self.fluxes[r] > 1e-6:
                reactions.append(self.reactions[r])
                fluxes.append(self.fluxes[r])
        return (reactions, fluxes)
    
    def get_conc(self):
        if not self.use_dG_f:
            return None
        
        # first create a set that contains all the CIDs that participate in active reactions
        active_cids = set()
        for r in range(len(self.reactions)):
            if (self.gammas[r] > 0.5 and self.fluxes[r] > 1e-6):
                active_cids = active_cids.union(self.reactions[r].sparse.keys())
        
        # get the concentrations for these CIDs in the solution
        cids = []
        concentrations = []
        for c in xrange(len(self.compounds)):
            if self.compounds[c].cid in active_cids:
                cids.append(self.compounds[c].cid)
                concentrations.append(pylab.exp(self.log_concentrations[c]))
        
        # return a pair with two lists, one of the CIDs and the other of the concentrations
        return (cids, concentrations)

    def get_margin(self):
        if not self.pCr:
            return None
        else:
            # the objective is (log_c_max - log_c_min)
            # so its exponent will give c_max/c_min, i.e. the margin
            return pylab.exp(self.cpl.solution.get_objective_value()) 
        
