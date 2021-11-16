import numpy as np
import sympy

class States:

    def __init__(self, J, P : int, nmax : int, lmax : int, Lmax : int):
        """
        Class for constructing possible wavefunction basis for
        a specific J^P baryon state

        Parameters
        ----------
        J : float
            Total angular momentum of baryon can be 1/2, 3/2 etc.
        P : int
            Parity of baryon J^P can be either +1 or -1
        nmax : int
            maximum value of n (radial quantum number)
        lmax : int 
            maximum value of orbital angular momenta lp and ly, 
            these couple to give total orbital angular momentum L
            lp + ly = L (vector coupling). There is no bound on
            lp and ly, so we need to constrain it.
        Lmax : int
            maximum value of L (where L = lp + ly)  
        """
        
        # Initial class variables
        self.J = sympy.Rational(J)
        self.P = P
        self.lmax = lmax
        self.Lmax = Lmax
        self.nmax = nmax

        """Lmax is equal to J + 3/2 because the highest L still needs to 
        give us J, therefore, L_highest - 3/2 = J, turning this 
        equation around we get L_highest = J + 3/2. The 3/2 comes
        from the highest value of the total spin that a baryon can
        possess.
        """
        self.physical_Lmax = self.J + sympy.Rational(3,2)
        if self.Lmax > self.physical_Lmax:
            self.Lmax = self.physical_Lmax

    
    def generate_all(self):

        # List for collecting possible combinations of [[lp,ly,L]]
        llist = []

        # Loop for populating llist 
        for lp in range(self.lmax + 1):
            for ly in range(self.lmax + 1):

                # Obtaining possible L's using vector sum lp + ly = L
                coupled_L_min = abs(lp - ly)
                coupled_L_max = lp + ly

                # Constraining L to Lmax
                if coupled_L_max > self.Lmax:
                    coupled_L_max = self.Lmax

                """Generate list of positive parity states, given
                that parity P, that is lp + ly is even. So, all 
                combinations that give lp + ly even, from lp and ly
                0 to lp and ly lmax for all possible couplings to L
                up to Lmax. 
                [lp, ly, |lp - ly|], [lp, ly, |lp - ly| + 1], ...
                [lp, ly, lp + ly if <= Lmax]
                """
                if self.P == 1:
                    if (lp + ly) % 2 == 0:
                        for i in np.arange(coupled_L_min, 
                                coupled_L_max+1):
                            llist.append([lp, ly, i])

                """Generate list of negative parity states, given
                that parity P, that is lp + ly is odd. So, all 
                combinations that give lp + ly odd, from lp and ly
                0 to lp and ly lmax for all possible couplings to L
                up to Lmax. 
                [lp, ly, |lp - ly|], [lp, ly, |lp - ly| + 1], ...
                [lp, ly, lp + ly if <= Lmax]
                """
                if self.P == -1:
                    if (lp + ly) % 2 != 0:
                        for i in np.arange(coupled_L_min, 
                                coupled_L_max+1):
                            llist.append([lp, ly, i])
        return llist

    def find_spin(self, I : int):
        """Determine the parity of the the spin x space product
        Based on Isospin, 0 means flavour = sym, 1 means
        flavour is antisymmetric, therefore, if I = 0, then
        spin x space must be symmetric and if I = 1 then
        spin x space must be antisymmetric.
        """
        if I % 2 == 0:
            spinSpace_parity = +1   # Symmetric
        else:
            spinSpace_parity = -1   # Anti symmetric 

        # Generate list of all possible [[lp, ly, L]'s]
        llist = self.generate_all()

        # Bin for storing possible spin states
        state_bin = []

        # For each element in llist...
        for el in llist:
            # Determine if space part of w.f. (lp) is sym or antisym
            lp = el[0]
            if lp % 2 == 0:
                space_parity = +1   # even = sym
            else:
                space_parity = -1   # odd  = antisym

            # Grab L from el
            L = el[2]

            # Calculate possible J's for |L-1/2|, .., L+1/2
            Jmin_12 = abs(L-sympy.Rational(1,2))
            Jmax_12 = L + sympy.Rational(1,2)
            Jlist12 = np.arange(Jmin_12, Jmax_12 + 1)

            # Calculate possible J's for |L-3/2|, .., L+3/2
            Jmin_32 = abs(L-sympy.Rational(3,2))
            Jmax_32 = L + sympy.Rational(3,2)
            Jlist32= np.arange(Jmin_32, Jmax_32 + 1)
            
            # determine parity of spin
            # spinSpace_parity = spin_parit x space_parity
            spin_parity = spinSpace_parity * space_parity

            # Collectors for the spin basis |s_12, s>
            A012 = 0    # antisymmetric |0, 1/2>
            S112 = 0    # symmetric |1, 1/2>
            S132 = 0    # symmetric |1, 3/2>
            
            # Check if we can have |0, 1/2> 
            if spin_parity == -1:
                if self.J in Jlist12:
                    A012 += 1

            # Check if we can have |1, 1/2> 
            if spin_parity == +1:
                if self.J in Jlist12:
                    S112 += 1

            # Check if we can have |1, 3/2> 
            if spin_parity == +1:
                if self.J in Jlist32:
                    S132 += 1

            # Append possibilities into state_bin
            state_bin.append([el, [A012, S112, S132]])
        return state_bin

    def build_basis(self, I):
        basis_list   = self.find_spin(I)
        r12          = sympy.Rational(1,2)
        r32          = sympy.Rational(3,2)
        basis_states = [[0,r12], [1,r12], [1,r32]]
        full_states  = []

        for el in basis_list:
            orbital = el[0]
            spin    = el[1]
            for i, state in enumerate(spin):
                if state == 1:
                    full_states.append([orbital, basis_states[i]])


        return full_states

    def add_radial(self, I):
        angular_part = self.build_basis(I)
        full_basis = []
        for state in angular_part:
            for nRho in range(self.nmax + 1):
                for nLam in range(self.nmax + 1):
                    lp  = state[0][0]
                    ly  = state[0][1]
                    L   = state[0][2]
                    s12 = state[1][0]
                    s   = state[1][1]
                    full_basis.append([nRho, lp, nLam, ly, L, s12, s])
        return full_basis

if __name__ == '__main__':
    s1 = States(1/2, +1, 1, 2, 2)
    print(s1.J)
    print(s1.P)
    print(s1.lmax)
    print(s1.Lmax)
    print(s1.generate_all())
    print(len(s1.generate_all()))
    print(s1.add_radial(0))
    print(len(s1.add_radial(0)))
