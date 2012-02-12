import pulp
import numpy as np
from pygibbs.thermodynamic_constants import R

class StoichiometricSolution(object):
    
    def __init__(self):
        self.reactions = None
        self.fluxes = None
        self.dG_r = None
        self.dG0_r = None
        
        self.compounds = None
        self.concentrations = None
        self.dG_f = None
        self.dG0_f = None

class Stoichiometric_LP(object):
    
    def __init__(self, name='Stoichiometric_LP', log_file=None):
        self.prob = pulp.LpProblem(name, pulp.LpMinimize)
        self.prob.solver = pulp.solvers.CPLEX(msg=False,
                                              logfilename='../res/cplex.log')
        
        self.S = None
        self.weights = None
        self.compounds = None
        self.reactions = None
        self.solution_index = 0
        self.flux_upper_bound = 100
        self.gamma_vars = None
        self.formation_vars = None
        self.pCr = None
        self.mtdf = None
        
        self.cids = None
        self.dG0_f = None
        self.RT = None
        

    def export(self, fname):
        self.prob.writeLP(fname)

    def add_stoichiometric_constraints(self, weights, S, compounds, reactions,
                                       net_reaction):
        """
            S is a NCxNR matrix where the rows represent compounds and the columns represent reactions.
            b is the linear constraint vector, and is NCx1
            
            compounds is the list of compounds from KEGG (NC long)
            reactions is a list of pairs of RID and direction (NR long)
        """
        self.weights = weights
        self.S = S
        self.compounds = compounds
        self.reactions = reactions
        
        assert S.shape[0] == len(self.compounds)
        assert S.shape[1] == len(self.reactions)

        self.cids = [c.cid for c in self.compounds]
        
        # reaction fluxes are the continuous variables
        self.flux_vars = pulp.LpVariable.dicts("Flux",
                                               [r.name for r in self.reactions],
                                               lowBound=0,
                                               upBound=self.flux_upper_bound,
                                               cat=pulp.LpContinuous)

        # add a linear constraint on the fluxes for each compound (mass balance)
        for c, compound in enumerate(self.compounds):
            reactions_with_c = np.nonzero(self.S[c, :])[0]
            if len(reactions_with_c) == 0:
                continue
            
            # If the compound participates in the net reaction, force the mass
            # balance to be -1 times the desired net reaction coefficient
            if compound.cid in net_reaction.sparse:
                mass_balance = net_reaction.sparse[compound.cid]
            else:
                mass_balance = 0.0
                
            # Sum of all fluxes involving this compounds times the stoichiometry
            # of their reactions
            cid_sum = pulp.lpSum([self.S[c, r] * self.flux_vars[self.reactions[r].name]
                                  for r in reactions_with_c])
            
            self.prob.addConstraint(cid_sum == mass_balance,
                                    "C%05d_mass_balance" % compound.cid)
        
        obj = pulp.lpSum([self.flux_vars[self.reactions[r].name]*weight
                          for (r, weight) in self.weights])
        self.prob.setObjective(obj)
        
    def add_milp_variables(self):
        # add boolean indicator variables for each reaction, and minimize their sum
        self.gamma_vars = pulp.LpVariable.dicts("Gamma",
                                                [r.name for r in self.reactions],
                                                cat=pulp.LpBinary)
        
        for r in self.reactions:
            # add a constrain for each integer variable, so that they will be indicators of the flux variables
            # using v_i - M*gamma_i <= 0
            f = self.flux_vars[r.name]
            g = self.gamma_vars[r.name]
            self.prob.addConstraint(f - self.flux_upper_bound * g <= 0,
                                    r.name + "_bound")
        
        obj = pulp.lpSum([self.flux_vars[self.reactions[r].name]*weight
              for (r, weight) in self.weights])
        self.prob.setObjective(obj)
        
    def add_reaction_num_constraint(self, max_num_reactions):
        if self.gamma_vars is None:
            raise Exception("Cannot add reaction no. constraint without the MILP variables")
        sum_gammas = pulp.lpSum(self.gamma_vars.values())
        self.prob.addConstraint(sum_gammas <= max_num_reactions, 'num_reactions')
    
    def add_dGr_constraints(self, thermo, pCr=False, MTDF=False, maximal_dG=0.0):
        """
            Create concentration variables for each CID in the database (at least in one reaction).
            If this compound doesn't have a dG0_f, its concentration will not be constrained.
            If it has a dG0_f, then either it has a specific concentration range (most co-factors do),
            or not - which is the common case. In this case, we either apply global constraints (e.g. 1uM - 10mM),
            or we use the pCr method to minimize the range needed for allowing the pathway thermodynamically.
        """
        
        if self.gamma_vars is None:
            raise Exception("Cannot add thermodynamic constraints without the MILP variables")
        
        if pCr and MTDF:
            raise Exception("Cannot optimize both the pCr and the MTDF")
        
        if pCr:
            self.pCr = pulp.LpVariable("pCr", lowBound=0, upBound=1e6,
                                       cat=pulp.LpContinuous)
            self.prob.setObjective(self.pCr)
        elif MTDF:
            self.mtdf = pulp.LpVariable("mtdf", lowBound=-1e6, upBound=1e6,
                                        cat=pulp.LpContinuous)
            self.prob.setObjective(self.mtdf)
        
        self.formation_vars = pulp.LpVariable.dicts("Formation",
            ["C%05d" % cid for cid in self.cids],
            lowBound=1e-6,
            upBound=1e6,
            cat=pulp.LpContinuous) 

        self.dG0_f = thermo.GetTransformedFormationEnergies(self.cids)[:, 0]
        self.RT = R * thermo.T
        
        for c, cid in enumerate(self.cids):
            dG0 = self.dG0_f[c]
            if np.isnan(dG0):
                # If the standard formation energy of this compound cannot be
                # determined, it's formation energy variable cannot be
                # constrained.
                continue
            
            c_var = self.formation_vars["C%05d" % cid]
            c_min, c_max = thermo.cid_to_bounds(cid, use_default=False)
                
            dG0_mid = dG0 + self.RT * np.log(thermo.c_mid)
            if self.pCr is not None and c_min is None:
                self.prob.addConstraint(c_var >= dG0_mid - self.RT * self.pCr,
                                        "C%05d_conc_minimum" % cid)
            else:
                c_var.lowBound = dG0 + self.RT * np.log(c_min or thermo.c_range[0])
                
            if self.pCr is not None and c_max is None:
                self.prob.addConstraint(c_var <= dG0_mid + self.RT * self.pCr,
                                        "C%05d_conc_maximum" % cid)
            else:
                c_var.upBound = dG0 + self.RT * np.log(c_max or thermo.c_range[1])

        for r in self.reactions:
            g = self.gamma_vars[r.name]
            dG_r = pulp.lpSum([self.formation_vars['C%05d' % cid]*coeff
                               for cid, coeff in r.sparse.iteritems()])
            
            if self.mtdf is not None:
                self.prob.addConstraint(dG_r <= self.mtdf + 1e6*(1 - g),
                                        r.name + "_thermo")
            else:
                self.prob.addConstraint(dG_r <= maximal_dG + 1e6*(1 - g),
                                        r.name + "_thermo")
            
    def add_localized_dGf_constraints(self, thermo):
        for r in self.reactions:
            dG0_r = thermo.GetTransfromedKeggReactionEnergies([r])[0, 0]
            if np.isnan(dG0_r): # no constraints on reactions with unknown dG0'
                continue
            for cid, coeff in r.sparse.iteritems():
                c_min, c_max = thermo.cid_to_bounds(cid, use_default=True)

                if coeff < 0:
                    dG0_r += coeff * R * thermo.T * np.log(c_max)
                else:
                    dG0_r += coeff * R * thermo.T * np.log(c_min)
            
            if dG0_r > 0:
                self.prob.addConstraint(self.flux_vars[r.name] == 0,
                                        r.name + "_irreversible")

    def solve(self, export_fname=None):
        self.prob.solve()
            
        if self.prob.status != pulp.LpStatusOptimal:
            return False

        self.fluxes = np.array([self.flux_vars[r.name].varValue
                                for r in self.reactions], ndmin=2).T
        
        if self.gamma_vars is not None:
            self.gammas = np.array([self.gamma_vars[r.name].varValue
                                    for r in self.reactions], ndmin=2).T
            self.gammas = np.logical_and(self.gammas, self.fluxes > 1e-6)
        else:
            self.gammas = self.fluxes > 1e-6

        self.concentrations = None
        if self.formation_vars is not None:
            self.dG_f = np.array([self.formation_vars["C%05d" % cid].varValue or np.nan
                                   for cid in self.cids])

            # if dG0_f is NaN for a compound, then so will the concentration be NaN
            self.concentrations = np.exp((self.dG_f - self.dG0_f) / self.RT)
            
            # multiply columns in S one-by-one, since there are NaN values
            # in 'formations' which will otherwise cause the result to be
            # all NaNs.
            self.dG_r = np.zeros((len(self.reactions)))
            self.dG0_r = np.zeros((len(self.reactions)))
            for r, reaction in enumerate(self.reactions):
                for cid, coeff in reaction.iteritems():
                    c = self.cids.index(cid)
                    self.dG0_r[r] += coeff * self.dG0_f[c]
                    self.dG_r[r] += coeff * self.dG_f[c]
        
        return True

    def ban_current_solution(self):
        """
            Use the binary reaction indicators (gammas) to form a constraint
            that will exclude the last solution from the MILP solution space.
            
            The trick is to add a constraint on the sum of the indicators of 
            the reactions participating in the solution, so that it will be 
            less than the number of reactions in that solution.
            
            In order to eliminate solutions which are equivalent (i.e. use 
            a reaction which differs only in the co-factors it uses) we need
            to find all rows in S which are identical to the rows in the solution
            and add their indicators to the sum as well.
        """
        if self.gamma_vars is None:
            raise Exception("Cannot ban a solution without the MILP variables")

        equivalence_set = []
        for r in np.where(self.gammas)[0]:
            for i in xrange(self.S.shape[1]):
                if np.all(self.S[:, i] == self.S[:, r]):
                    equivalence_set.append(self.gamma_vars[self.reactions[i].name])
                    
        N = np.sum(self.gammas)
        self.prob.addConstraint(pulp.lpSum(equivalence_set) <= N - 1,
                                "solution_%03d" % self.solution_index)
        self.solution_index += 1
    
    def get_margin(self):
        if self.pCr is None:
            return None
        else:
            # the objective is (log_c_max - log_c_min)
            # so its exponent will give c_max/c_min, i.e. the margin
            return np.exp(2*self.pCr.varValue) 

    def get_mtdf(self):
        if self.mtdf is None:
            return None
        else:
            return self.mtdf.varValue

    def get_active_reaction_data(self):
        solution = StoichiometricSolution()

        active_r = np.where(self.gammas)[0]
        solution.reactions = [self.reactions[r] for r in active_r]

        active_cids = set()
        for reaction in solution.reactions:
            active_cids.update(reaction.get_cids())
        active_c = [self.cids.index(cid) for cid in sorted(active_cids)]

        solution.fluxes = self.fluxes[active_r, :]
        solution.compounds = [self.compounds[c] for c in active_c]
        if self.formation_vars is not None:
            solution.dG0_f = self.dG0_f[active_c]
            solution.dG_f = self.dG_f[active_c]
            solution.concentrations = self.concentrations[active_c]
            solution.dG0_r = self.dG0_r[active_r]
            solution.dG_r = self.dG_r[active_r]
        
        return solution

