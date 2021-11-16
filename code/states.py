import numpy as np

class StateMatrix:

    def __init__(self, gammaMax: int, Lambda: int, remove_ls=False):
        self.gammaMax   = gammaMax
        self.Lambda     = Lambda
        self.nlist      = self.n_list()
        self.llist      = self.l_list()
        self.all_states = self.all_states()
        if remove_ls == True:
            self.remove_ls()

    # Create list of possible n's
    def n_list(self):
        # bin to collect possible values of n
        temp_nlist = []
        # start with n = 0
        n = 0
        # while loop
        while 2*n <= self.gammaMax:
            # while the condition is true append n to bin
            temp_nlist.append(n)
            # increment n by 1
            n += 1
        # return the list of possible n's
        return temp_nlist

    # Create list of possible l's
    def l_list(self):
        # bin to collect possible values of l
        temp_llist = []
        # start with l = 0
        l = 0
        # while loop
        while l <= self.gammaMax:
            # while the condition is true append l to bin
            temp_llist.append(l)
            # increment l by 1
            l += 1
        # return the list of possible l's
        return temp_llist

    # Create all possible combinations with allowed energy (nonphysical)
    def all_states(self):
        list_of_states = []
        for nrho in self.nlist:
            for lrho in self.llist:
                for nlam in self.nlist:
                    for llam in self.llist:
                        # Calculate energy of state
                        energy = 2*nrho + lrho + 2*nlam + llam
                        # Obtain possible total angular momenta
                        Lambdas = [i for i in range(abs(lrho-llam),
                            lrho+llam+1)]
                        # If energy of state is equal or below gammaMax
                        if energy <= self.gammaMax:
                            # And if l's can couple to Lambda, append
                            if self.Lambda in Lambdas:
                                state = [nrho, lrho, nlam, llam]
                                list_of_states.append(state)
        return list_of_states

    # remove all states with non-zero l's
    #TODO remove according to key, e.g. (0,0) would remove all lrho=0
    # and llam=0 states etc.
    def remove_ls(self):
        states = []
        for state in self.all_states:
            if state[1] == 0 and state[3] == 0:
                states.append(state)
        self.all_states = states

    def make_matrix(self):
        state_array = np.zeros
        for state1 in self.all_states:
            for state2 in self.all_states:




if __name__ == '__main__':

    sm1 = StateMatrix(6,0)
    print("All states")
    for state in sm1.all_states:
        print(state)
    sm1.remove_ls()
    print("States with l = 0")
    for state in sm1.all_states:
        print(state)
