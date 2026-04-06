# Introduction
# Hidden Markov Model — Part 2: Forward, Backward, and Posterior Decoding

## Introduction

This project is the second part of a two-part implementation of a Hidden Markov Model (HMM) in Python. Part 1 covered the **Viterbi algorithm**, which finds the single most likely sequence of hidden states given a sequence of observations. This part builds on that foundation by implementing three closely related algorithms: the **Forward algorithm**, the **Backward algorithm**, and **Posterior Decoding**.

## Background

A Hidden Markov Model is a statistical model used to represent systems where the true underlying states are not directly observable, but produce observable outputs. An HMM is defined by:

- **Hidden states** — the unobservable conditions of the system
- **Observations** — the symbols or events that can be seen
- **Initial probabilities (betas)** — the probability of starting in each state
- **Transition probabilities** — the probability of moving from one state to another
- **Emission probabilities** — the probability of a state producing a given observation

## What This Part Implements

### Forward Algorithm
The Forward algorithm computes the probability of observing a given sequence of emissions up to time `t`, while being in a particular hidden state at time `t`. It works left to right through the observation sequence, accumulating probabilities at each step. The total probability of the entire observation sequence can be derived from its final column.

### Backward Algorithm
The Backward algorithm computes the probability of observing the remaining emissions from time `t+1` to the end of the sequence, given a particular hidden state at time `t`. It mirrors the Forward algorithm but works right to left through the observation sequence.

### Posterior Decoding
Posterior Decoding combines the Forward and Backward matrices to compute, at each position in the sequence, the probability of being in each hidden state given **all** observations — both past and future. The most likely state at each position is then selected independently. This differs from the Viterbi algorithm, which finds the globally most likely path through all states jointly.

## Implementation Notes

All three algorithms are implemented as methods of the `HMM` class established in Part 1. To avoid floating-point underflow from multiplying many small probabilities together, all computations are performed in **log space**. Multiplications become additions, and summations of probabilities are handled using `numpy`'s `logaddexp` function, which computes `log(exp(a) + exp(b))` in a numerically stable way.
# Pseudocode
```
FORWARD(observations):

    # Initialize matrix: rows = states, cols = observations
    forward_mat = empty matrix of shape (num_states, len(observations))

    # Fill first column using betas (initial probs) + emission
    for each state in self.states:
        forward_mat[state][0] = log(self.betas[state]) + log(state.emission_probs[observations[0]])

    # Fill remaining columns left to right
    for obs_ind from 1 to len(observations)-1:
        for each state in self.states:

            # Call _get_prev_state_options() to get list of:
            # forward_mat[prev_state][obs_ind-1] + log(t_mat[prev_state][state]) + log(emission)
            options = _get_prev_state_options(obs=obs_ind, observations, mat=forward_mat, mat_row=state)

            # SUM all options together in log space (logaddexp = log(sum(exp(options))))
            forward_mat[state][obs_ind] = logaddexp.reduce(options)

    return forward_mat

Helper function for Backward
_get_future_options(obs_ind, observations, mat, mat_row):

    # For each next state, compute:
    # log(t_mat[current_state → next_state]) + log(emission of next_state at obs_ind) + mat[next_state][obs_ind]
    options = [ log(t_mat[mat_row][next_state])
                + log(next_state.emission_probs[observations[obs_ind]])
                + mat[next_state][obs_ind]
                for each next_state ]

    return options                              
                              
BACKWARD(observations):

    # Initialize matrix: rows = states, cols = observations
    backward_mat = empty matrix of shape (num_states, len(observations))

    # Fill last column with 0 (which is log(1), meaning no future observations)
    for each row in backward_mat:
        row[-1] = 0

    # Fill remaining columns right to left
    for obs_ind from len(observations)-1 down to 1:
        for each state in self.states:

            # Call _get_future_options() to get list of:
            # log(t_mat[state → next]) + log(emission of next at obs_ind) + backward_mat[next][obs_ind]
            options = _get_future_options(obs_ind=obs_ind, observations, mat=backward_mat, mat_row=state)

            # SUM all options in log space, store one column to the LEFT
            backward_mat[state][obs_ind - 1] = logaddexp.reduce(options)

    return backward_mat


POSTERIOR_DECODING(observations):

    # Run both algorithms
    forward_mat  = FORWARD(observations)     # shape: (num_states, T)
    backward_mat = BACKWARD(observations)    # shape: (num_states, T)

    # Get log P(observations) by summing last column of forward_mat in log space
    log_prob_obs = logaddexp.reduce(forward_mat[:, -1])

    # Compute gamma for every state at every position
    # In log space: gamma = forward + backward - log_prob_obs
    gamma = forward_mat + backward_mat - log_prob_obs

    # At each time step, pick the state with the highest gamma value
    best_state_indices = argmax(gamma, axis=0)   # one index per time step

    # Map indices back to State objects using self.states
    return [ self.states[i] for i in best_state_indices ]                              
```

# Successes
One of our biggest strengths was that our team developed a strong understanding of the underlying logic behind the Forward, Backward, and Posterior Decoding algorithms. Building on our work from the Viterbi implementation, we were able to clearly visualize how the probabilities are processed through the model in both directions and how combining these perspectives provides a more complete picture of state likelihoods. We became comfortable working in log space and understood why functions like logaddexp are necessary for numerical stability. We also maintained the level of communication we had last week, which helped us all stay up to date with progress any of us made despite the busy schedules we all had. 

# Struggles
One of the main challenges we faced was coordinating time to meet as a group. The busy nature of our scheudles made it hard to find a time where we all could meet, but we still did our best through frequent updates and messages on teams. Additionally, translating the mathematical concepts and pseudocode into working Python code proved to be more confusing than we expected. While we understood the concepts and had working code that we checked through careful debugging, we found ourselves stumped trying to make sense of the results. Eventually, we did verify our results and found it was working as intented so that wasn't too big of a hurdle for us.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

## Shameem Shahib
Compared to the first part of the HMM, I found the concepts for the forward and backward algorithms to be a little easier to visualize and understand. Posterier decoding was interesting to tackle as well, and I found it to be enriching to put the forward and backward together. Working with my group continued to be a big strength. Being able to talk through the algorithm and troubleshoot together made a big difference in my learning. Overall, I'm glad with how progress is going, and the final code we put together for this week.

# Generative AI Appendix
Claude Sonnet 4.6 was used to check our logic in our code and ensure our results matched the logic behind the algorithm.
