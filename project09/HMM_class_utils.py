import numpy as np
from numbers import Number
from collections.abc import Iterable
from math import log

class State:
    """Hidden state for HMM"""
    def __init__(self, name:str, emissions: list, probabilities: list[Number], transitions:dict[str:float]):
        self.name = name
        self.emissions = set(emissions) # Create a set of emission labels for faster future checks
        self.emission_probs = dict(zip(emissions, probabilities))  # Create state emission: probabilities pairs from input
        self.transitions = transitions
        self.total_emission_prob = sum(probabilities)  # Save sum of emission probabilities
        if self.total_emission_prob != 1:  # Confirm emission probabilities sum to 1, else raise error
            raise ValueError("Emission probabilities do not sum to 1")
        
    def __repr__(self):
        return self.name

    def add_emission(self, emission, probability: Number):
        """
        Add emission:probability pair to State emissions dictionary

        Args:
            name (str): name of State
            emissions (list): observations as seen in data
            probabilities (list): probabilities of emissions

        Returns:
            None:
        """
        self.emissions.add(emission)
        self.emission_probs[emission] = probability  # Create new emission: probability pair
        self.total_emission_prob += probability # Update sum of emission probabilities

        # Check if emission probabilities are still equal to 1, else get relative probabilities
        if self.total_emission_prob > 1: 
            print(f"\nWARNING: sum of State: {self.name} emission probabilities has exceeded 0.\n"
                  f"Refactoring to maintain relative probabilities with sum of 1\n")
            for emit in self.emissions:
                self.emission_probs[emit] = self.emission_probs[emit] / self.total_emission_prob
    
class HMM():
    def __init__(self, name: str, betas: dict[State:float], emissions: set, states: list[State]):
        self.name = name  # Name for identification
        self.states = states  # List of state objects in HMM
        self.t_mat = self.build_transition_mat_from_states()  # Transition matrix of states -> states in HMM
        self.emissions = set(emissions) # Set of emissions in HMM
        self.betas = betas

    def __repr__(self):
        return(f"{self.name}\n"
               f"{self.emissions}\n"
               f"{self.states}\n"
               f"{self.t_mat}\n"
               f"{self.betas}")
    
    def build_transition_mat_from_states(self):
        """
        Build a transition matrix from the transition probabilities of states in HMM

        Returns:
            np.array: Matrix of transitions from states transition dicts (transition[from state][to state])
        """
        t_mat = [[] for state in self.states]  # Initialize empty list of lists for transition matrix

        # Fill transition matrix list with transition probabilities from states in HMM.states
        for row, state in enumerate(self.states):
            for trans in self.states:
                t_mat[row].append(state.transitions[trans.name])

        return np.array(t_mat)

    def add_state(self, name:str = None, emissions: list = [], probabilities: list = [], state:State = None):
        """
        Add a State to the HMM either by providing a State or arguments for a state

        Args:
            name (str): name of new State if not using existing State
            emissions (list): list of State emissions
            probabilities (list): list of probabilities corresponding to list of emissions
            state (State): existing State object to add instead of creating new from other arguments

        Examples:
            HMM.add_state(name = "my_state", emissions = ["A", "B", "C"], probabilities = [0.2, 0.3, 0.5])  
            HMM.add_state(state = existing_state)
        """
        if state is None: # Check if pre-existing State given

            # Check if name provided to create new State, else raise ValueError
            if name is None:
                raise ValueError("No name given to create new State object from arguments")
            state = State(name=name, emissions=emissions, probabilities=probabilities)

        self.states.append(state) # Add state to HMM list of States

        # Compare state emission to HMM emissions and add missing state emissions to HMM
        for emission in state.emissions:
            if not emission in self.emissions:
                self.emissions.add(emission)

        # Compare HMM emissions to state emissions and add missing HMM emissions with probability 0
        for emission in self.emissions:
            for state in self.states:
                if not emission in state.emissions:
                    state.add_emission(emission=emission, probability=0)

    def _get_prev_state_options(self, obs:int, observations: Iterable, mat:np.ndarray, mat_row:int) -> list:
        """
        Create a list of previous path * transition * emission probabilities for each potential state transition

        Args:
            obs (int): observation sequence index
            observations (Iterable): Iterable object of observations for model
            mat (np.ndarray): matrix up to current observation
            mat_row (int):row index for transitioning in to

        Returns:
            list: list of potential path probabilities
        """

        options = [state_row[obs-1] + log(self.t_mat[trans_state][mat_row]) + log(self.states[mat_row].emission_probs[observations[obs]])
                    for trans_state, state_row in enumerate(mat)]
        
        return options
    
    def viterbi(self, observations: Iterable) -> list[State]:
        """
        Calculate most likely state at each observation

        Args:
            observations (Iterable): iterable of observations

        Returns:
            list: list of which state is most likely at each observation
        """
        # Initialize empty arrays to hold path probabilities and traceback
        vit = np.zeros((len(self.states), len(observations)))
        traceback = np.zeros((len(self.states), len(observations)), dtype=int)

        # Set 0 index of trace array rows to index not in states list to recognize as stop signal
        for row in traceback:
            row[0] = len(self.states)

        # Set 0 index of probability array rows to beta probability * emission probability
        for state, row in enumerate(vit):
            row[0]=log(self.betas[self.states[state]]) + log(self.states[state].emission_probs[observations[0]])

        # Populate each row at each position with max probability of cumulative prob * transition prob * emission prob
        for obs in range(1, len(observations)):
            for state, row in enumerate(vit):
                options = self._get_prev_state_options(obs=obs, observations=observations, mat=vit, mat_row=state)
                row[obs] = max(options)  # Take highest probability of possible paths to this emission
                traceback[state][obs] = np.argmax(options) #.index(row[obs]) # Save state path extended from as index in states list

        end_states = [row[-1] for row in vit]  # List of final probabilities in viterbi probability matrix
        max_path = np.argmax(end_states)  # Get row with highest final probability

        # Start from end of traceback matrix at highest probability path and follow traceback to beginning
        index = -1
        trace = max_path
        print("trace=", trace)
        print("index =", index)
        state_path = [max_path]
        while traceback[trace][index] != len(self.states):
            trace = traceback[trace][index]
            state_path.append(trace)
            index -= 1

        for ind, code in enumerate(state_path):
            state_path[ind] = self.states[code]

        state_path = list(reversed(state_path))

        print(traceback)
        print(vit)
        return state_path
    
    def forward(self, observations: Iterable) -> np.ndarray:
        forward_mat = np.ndarray((len(self.states), len(observations)))

        for state, row in enumerate(forward_mat):
            beta = self.betas[self.states[state]]
            emission = self.states[state].emission_probs[observations[0]]
            row[0] = log(beta) + log(emission)

        for obs_ind in range(1, len(observations)):
            for state, row in enumerate(forward_mat):
                options = self._get_prev_state_options(obs=obs_ind, observations=observations, mat=forward_mat, mat_row=state)
                total_prob = np.logaddexp.reduce(options)
                row[obs_ind] = total_prob

        return forward_mat
    
    def backward(self, observations: Iterable) -> np.ndarray:
        pass

if __name__ == "__main__":
    observations = "ABC"

    my_name = "my_state"
    my_emissions = ["A", "B", "C"]
    my_probs = [0.3, 0.2, 0.5]
    my_state = State(name="my_state", emissions=my_emissions, probabilities=my_probs, transitions={"my_state":0.7, "my_state2":0.3})
    my_state2 = State(name = "my_state2", emissions=["A", "B", "C"], probabilities=[0.2, 0.7, 0.1], transitions={"my_state2":0.9, "my_state":0.1})
    my_HMM = HMM(name="My_HMM", betas={my_state:0.5, my_state2:0.5}, emissions={"A", "B", "C"}, states=[my_state, my_state2])
    print(my_HMM)

    print(my_HMM.viterbi(observations=observations))

    fwmat = my_HMM.forward(observations=observations)
    print(fwmat[0][-1]+fwmat[1][-1])