import numpy as np
from numbers import Number
from collections.abc import Iterable
from math import log

#from project09.HMM_notebook import backward_mat


class State:
    """Hidden state for HMM"""
    def __init__(self, name: str, emissions: list, probabilities: list[Number], transitions: dict[str:float]):
        self.name = name
        self.emissions = set(emissions)  # Create a set of emission labels for faster future checks
        self.emission_probs = dict(zip(emissions, probabilities))  # Create state emission: probabilities pairs from input
        self.transitions = transitions
        self.total_emission_prob = sum(probabilities)  # Save sum of emission probabilities
        self.total_emission_prob = sum(probabilities)
        # Confirm emission probabilities sum to 1, else raise error
        if not np.isclose(self.total_emission_prob, 1.0):
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
        self.total_emission_prob += probability  # Update sum of emission probabilities

        # Check if emission probabilities are still equal to 1, else get relative probabilities
        if self.total_emission_prob > 1:
            print(f"\nWARNING: sum of State: {self.name} emission probabilities has exceeded 0.\n"
                  f"Refactoring to maintain relative probabilities with sum of 1\n")
            for emit in self.emissions:
                self.emission_probs[emit] = self.emission_probs[emit] / self.total_emission_prob
            self.total_emission_prob = sum(self.emission_probs.values())


class HMM():
    def __init__(self, name: str, betas: dict[State:float], emissions: set, states: list[State]):
        self.name = name  # Name for identification
        self.states = states  # List of state objects in HMM
        self.t_mat = self.build_transition_mat_from_states()  # Transition matrix of states -> states in HMM
        self.emissions = set(emissions)  # Set of emissions in HMM
        self.betas = betas
        self._synchronize_emission()

    def __repr__(self):
        return (f"{self.name}\n{self.emissions}\n{self.states}\n{self.t_mat}\n{self.betas}")


    @staticmethod
    def _safe_log(x:float) -> float:
        """
        Return log(x) if x > 0, else -inf
        """
        if x <= 0:
            return float("-inf")
        return log(x)

    def _synchronize_emission(self):
        """
        Ensure every state has every emission in the MM.
        Missing emission are added with probability 0
        """
        for state in self.states:
            for emission in state.emissions:
                self.emissions.add(emission)

        for emission in self.emissions:
            for state in self.states:
                if emission not in state.emissions:
                    state.emissions.add(emission)
                    state.emission_probs[emission] = 0.0

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

        return np.array(t_mat, dtype=float)

    def add_state(self, name: str = None, emissions: list = None, probabilities: list = None, transitions: dict[str, float] = None, state: State = None):
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
        if emissions is None:
            emissions = []
        if probabilities is None:
            probabilities = []

        if state is None:
            if name is None:
                raise ValueError("No name given to create new State object from argument")
            if transitions is None:
                raise ValueError("No name given to create new State object from argument")
            state = State(name=name, emissions=emissions, probabilities=probabilities, transitions=transitions)

        self.states.append(state)
        self._synchronize_emission()
        self.t_mat = self.build_transition_mat_from_states()

    def _get_prev_state_options(self, obs: int, observations: Iterable, mat: np.ndarray, mat_row: int) -> list:
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
        # For each possible path probability into the given state at the current observation store that cumulative probability in a list of paths
        options = [
            mat[trans_state][obs - 1]
            + self._safe_log(self.t_mat[trans_state][mat_row])
            + self._safe_log(self.states[mat_row].emission_probs[observations[obs]])
            for trans_state in range(len(self.states))
        ]

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
        vit = np.full((len(self.states), len(observations)), float("-inf"))
        traceback = np.zeros((len(self.states), len(observations)), dtype=int)

        # Set 0 index of trace array rows to index not in states list to recognize as stop signal
        for row in traceback:
            row[0] = len(self.states)

        # Set 0 index of probability array rows to beta probability * emission probability
        for state, row in enumerate(vit):
            row[0] = self._safe_log(self.betas[self.states[state]]) + self._safe_log(self.states[state].emission_probs[observations[0]])

        # Populate each row at each position with max probability of cumulative prob * transition prob * emission prob
        for obs in range(1, len(observations)):
            for state, row in enumerate(vit):
                options = self._get_prev_state_options(obs=obs, observations=observations, mat=vit, mat_row=state)
                row[obs] = max(options)  # Take highest probability of possible paths to this emission
                traceback[state][obs] = np.argmax(options)  # Save state path extended from as index in states list

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
        """
        Create forward matrix from sequence of observations

        Args:
            observations (Iterable): sequence of observations

        Returns:
            np.ndarray: matrix of cumulative path probabilities for each state at each position
        """
        # Initialize empty forward matrix with a row for each state and column for each observation
        forward_mat = np.full((len(self.states), len(observations)), float("-inf"))

        # Set each row in first column equal to log(beta * emission) for corresponding state
        for state, row in enumerate(forward_mat):
            beta = self.betas[self.states[state]]
            emission = self.states[state].emission_probs[observations[0]]
            row[0] = self._safe_log(beta) + self._safe_log(emission)

        # For each observation and each possible state at that observation get the cumulative path probabilities into that state
        for obs_ind in range(1, len(observations)):
            for state, row in enumerate(forward_mat):
                options = self._get_prev_state_options(obs=obs_ind, observations=observations, mat=forward_mat, mat_row=state)
                row[obs_ind] = np.logaddexp.reduce(options)

        return forward_mat

    def _get_future_options(self, obs_ind: int, observations: Iterable, mat: np.ndarray, mat_row: int) -> list:
        """
        Create a list of cumulative path probabilities for a backward matrix

        Args:
            obs_ind (int): index of observation in observation sequence
            observations (Iterable): sequence of observations
            mat (np.ndarray): current backward matrix
            mat_row (int): current row of matrix being calculated for position obs_ind

        Returns:
            list: list of cumulative path probabilities up to current row and position
        """
        # For each possible sum of current path probability, transition probability, and emission probability, create a corresponding element in a list of path options
        options = [
            self._safe_log(self.t_mat[mat_row][next_state])
            + self._safe_log(self.states[next_state].emission_probs[observations[obs_ind]])
            + mat[next_state][obs_ind]
            for next_state in range(len(self.states))
        ]

        return options

    def backward(self, observations: Iterable) -> np.ndarray:
        """
        Create backward matrix from sequence of observations

        Args:
            observations (Iterable): sequence of observations

        Returns:
            np.ndarray: backward matrix
        """
        # Initialize backward matrix
        backward_mat = np.full((len(self.states), len(observations)), float("-inf"))

        # Set last column values as 0 (log space equivalent of 1)
        for row in backward_mat:
            row[-1] = 0.0

        # Iterate from back of observation sequence to front, populating array
        for obs_ind in range(len(observations)-1, 0, -1):
            for state, row in enumerate(backward_mat):
                options = self._get_future_options(observations=observations, obs_ind=obs_ind, mat=backward_mat, mat_row=state)
                row[obs_ind-1] = np.logaddexp.reduce(options)

        return backward_mat

    def posterier_decoding(self, observations: Iterable):
        """
        Calculate the most likely observation for a given state using the forward-backward algorithm

        Args:
            observations (Iterable): sequence of observations

        Returns:
            (np.ndarray, list): posterior decoding matrix and list of most likely state at each observation position
        """
        # Create the forward and backward matrices
        forward_mat = self.forward(observations)
        backward_mat = self.backward(observations)

        # Calculate the log
        log_prob = np.logaddexp.reduce(forward_mat[:, -1])

        # Create posterior matrix
        posterior_mat = forward_mat + backward_mat - log_prob

        # Determine the most likely state at a given position
        state_idx = np.argmax(posterior_mat, axis=0)
        state_path = [self.states[i] for i in state_idx]

        return posterior_mat, state_path

    def _expectation_values(self, observations: Iterable):
        "Compute shared E-step quantities for one sequence."

        observations = list(observations)
        n_states = len(self.states)
        T = len(observations)

        forward_mat = self.forward(observations)
        backward_mat = self.backward(observations)
        log_prob = np.logaddexp.reduce(forward_mat[:, -1])

        gamma = np.exp(forward_mat + backward_mat - log_prob)

        xi = np.zeros((T-1, n_states, n_states), dtype=float)
        for t in range (T-1):
            for i in range(n_states):
                for j in range(n_states):
                    log_xi = (
                            forward_mat[i, t]
                            + self._safe_log(self.t_mat[i, j])
                            + self._safe_log(self.states[j].emission_probs[observations[t+1]])
                            + backward_mat[j, t+1]
                            - log_prob)
                    xi[t, i, j] = np.exp(log_xi)

        return forward_mat, backward_mat, log_prob, gamma, xi

    def sequence_log_likelihood(self, observations: Iterable) -> float:
        " Compute log-likelihood of a sequence of observations"
        _, _, log_prob, _, _ = self._expectation_values(observations)
        return log_prob

    def gamma_matrix(self, observations: Iterable) -> np.ndarray:
        """return gamma matrix"""
        _, _, _, gamma, _ = self._expectation_values(observations)
        return gamma

    def xi_tensor(self, observations: Iterable) -> np.ndarray:
        """Return xi matrix"""
        _, _, _, _,xi = self._expectation_values(observations)
        return xi

    def posterior_decoding(self, observations: Iterable): #new approach
        "Compute posterior decoding using gamma matrix"

        gamma = self.gamma_matrix(observations)
        posterior_mat = np.where(gamma > 0, np.log(gamma), float("-inf"))

        state_idx = np.argmax(posterior_mat, axis=0)
        state_path = [self.states[i] for i in state_idx]

        return posterior_mat, state_path

    def baum_welch(self, sequences, max_iter=100, tol=1e-4, pseudocount=1e-6):

        n_states = len(self.states)

        emission_symbols = sorted(self.emissions)
        emission_to_idx = {}
        for i, symbol in enumerate(emission_symbols):
            emission_to_idx[symbol] = i

        history = []
        prev_loglik = float("-inf")

        for iteration in range(max_iter):

            beta_counts = np.zeros(n_states)
            trans_counts = np.zeros((n_states, n_states))
            trans_denom = np.zeros(n_states)
            emit_counts = np.zeros((n_states, len(emission_symbols)))
            emit_denom = np.zeros(n_states)

            total_loglik = 0.0

            # -------- E-STEP --------
            for observations in sequences:

                observations = list(observations)

                _, _, log_prob, gamma, xi = self._expectation_values(observations)

                total_loglik += log_prob

                beta_counts += gamma[:, 0]

                trans_counts += np.sum(xi, axis=0)
                trans_denom += np.sum(gamma[:, :-1], axis=1)

                for t in range(len(observations)):
                    obs = observations[t]
                    obs_idx = emission_to_idx[obs]
                    emit_counts[:, obs_idx] += gamma[:, t]

                emit_denom += np.sum(gamma, axis=1)

            # -------- M-STEP --------

            # update initial probs
            beta_total = np.sum(beta_counts) + pseudocount * n_states
            for i in range(n_states):
                self.betas[self.states[i]] = (beta_counts[i] + pseudocount) / beta_total

            # update transitions
            for i in range(n_states):
                denom = trans_denom[i] + pseudocount * n_states
                for j in range(n_states):
                    prob = (trans_counts[i, j] + pseudocount) / denom
                    self.t_mat[i, j] = prob
                    self.states[i].transitions[self.states[j].name] = prob

            # update emissions
            for i in range(n_states):
                denom = emit_denom[i] + pseudocount * len(emission_symbols)
                for symbol in emission_symbols:
                    k = emission_to_idx[symbol]
                    prob = (emit_counts[i, k] + pseudocount) / denom
                    self.states[i].emission_probs[symbol] = prob

            history.append(total_loglik)

            # convergence check
            if len(history) > 1:
                if abs(total_loglik - prev_loglik) < tol:
                    break

            prev_loglik = total_loglik

        return history


if __name__ == "__main__":
    observations = "ABACB"

    my_name = "my_state"
    my_emissions = ["A", "B", "C"]
    my_probs = [0.3, 0.2, 0.5]
    my_state = State(name="my_state", emissions=my_emissions, probabilities=my_probs, transitions={"my_state": 0.7, "my_state2": 0.3})
    my_state2 = State(name="my_state2", emissions=["A", "B", "C"], probabilities=[0.2, 0.7, 0.1], transitions={"my_state2": 0.9, "my_state": 0.1})
    my_HMM = HMM(name="My_HMM", betas={my_state: 0.5, my_state2: 0.5}, emissions={"A", "B", "C"}, states=[my_state, my_state2])
    print(my_HMM)

    print(my_HMM.viterbi(observations=observations))

    fwmat = my_HMM.forward(observations=observations)
    bkmat = my_HMM.backward(observations=observations)
    print(fwmat)
    print("np.logaddexp(fwmat[0][-1], fwmat[1][-1])", np.logaddexp(fwmat[0][-1], fwmat[1][-1]))
    print(bkmat)
    print("np.logaddexp(bkmat[0][-1], bkmat[1][-1])", np.logaddexp(bkmat[0][0], bkmat[1][0]))
    post_mat, post_path = my_HMM.posterier_decoding(observations=observations)
    print(post_mat)
    print(post_path)

    log_prob_fw = np.logaddexp.reduce(fwmat[:, -1])
    log_prob_bk = np.logaddexp.reduce(bkmat[:, 0])
    print("log_P from forward: ", log_prob_fw)
    print("log_P from backward:", log_prob_bk)

    # --- Assignment example for Baum-Welch ---
    obs = ["GGCACTGAA", "ATGCAATGC", "AATGCCTGA"]
    sequences = [list(seq) for seq in obs]

    state_H = State(
        name="H",
        emissions=["A", "C", "G", "T"],
        probabilities=[0.2, 0.3, 0.3, 0.2],
        transitions={"H": 0.6, "L": 0.4}
    )

    state_L = State(
        name="L",
        emissions=["A", "C", "G", "T"],
        probabilities=[0.3, 0.2, 0.2, 0.3],
        transitions={"H": 0.3, "L": 0.7}
    )

    my_HMM = HMM(
        name="GC_Model",
        betas={state_H: 0.5, state_L: 0.5},
        emissions={"A", "C", "G", "T"},
        states=[state_H, state_L]
    )

    print("\n--- Initial Parameters ---")
    print("Initial betas:")
    for state in my_HMM.states:
        print(state.name, my_HMM.betas[state])

    print("\nInitial transition matrix:")
    print(my_HMM.t_mat)

    print("\nInitial emission probabilities:")
    for state in my_HMM.states:
        print(state.name, state.emission_probs)

    history = my_HMM.baum_welch(sequences, max_iter=20, tol=1e-4, pseudocount=1e-6)

    print("\n--- Baum-Welch Log-Likelihood History ---")
    print(history)

    print("\n--- Updated Parameters After Baum-Welch ---")
    print("Updated betas:")
    for state in my_HMM.states:
        print(state.name, my_HMM.betas[state])

    print("\nUpdated transition matrix:")
    print(my_HMM.t_mat)

    print("\nUpdated emission probabilities:")
    for state in my_HMM.states:
        print(state.name, state.emission_probs)

    print(sum(my_HMM.states[0].emission_probs.values()))
    print(sum(my_HMM.states[1].emission_probs.values()))

