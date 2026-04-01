# Introduction
Hidden Markov Models (HMMs) can provide probabilistic characterization for a sequence of observations that would otherwise be hard to parse.
For example predicting the seasons of the year based on weather patterns. For this project we will be using Viterbi's algorithm to select the most likely hidden states underlying a sequence of observations.

# Pseudocode
```
**class State(name, emissions)**
Initialize State:
    INPUT: name, list of emissions, list of probabilities, transition dictionary
    SET state name
    STORE emissions as a set
    CREATE dictionary mapping emissions → probabilities
    STORE transition probabilities
    
    CALCULATE total emission probability
    
    IF total emission probability ≠ 1:
        RAISE error
Add Emission
    Input: emission, probability
    
    Add emission to emission set
    Add (emission -> probability) to dictionary
    update total emission probability
    
    If total emission probabily > 1:
        print warning
        normalize all emission probabilities so sum = 1
        
**class HMM(emissions, states)**
Initialize HMM:
    Input: name, initial probabilities, emissions, list of states
    
    Store name
    store states
    store emissions
    store initial state probabilities
    
    built transtion matrix from state transition dictionaries.
    
Build transition matrix:
    create empty HMM with list of emissions
    for each state i: 
        for each state j:
            add transition probability from state i to state j
    return matrix
    
Add State:
    Input: either (state object) or (name, emissions, probabilities)
 
    If no state object provided:
        create new state
        
    Add state to HMM
    
    For each emission in state:
        if emission not in HMM emissions:
            add it
    
    For each emission in HMM:
        for each state:
            if state does not contain emission:
                  add emission with probability 0

**Viterbi Algorithm:**
Initialize
    Input: observations
    
    create matrix V (state * observation) -> store probabilities
    create matrix traceback (state * observations) -> store paths
    
    for each state:
        V[state][0] = log(initial probability) + log(emission probability of first observation)
        traceback[state][0] = STOP marker
        
Dynamic Programming Step
    For each observation t from 1 to end:
        calculate all possible paths:
            previous probability
            +log(transition probability)
            +log(emssion probability)
         
        select maximum value
        store in V[s][t]
        
        store index of best previous state in traceback[s][t]

Termination:
    Find state with maximum probability at last observation
    set this as final state
    
Backtracking:
    Initialize path with final state
    
    while previous state exists:
        follow traceback matrix backward
        add states to path
    
    reverse path to correct order
    
Return most likely sequence of states
    
```

# Successes
Description of the team's learning points

# Struggles
One of the main challenges we faced was aligning our understanding of the algorithm and ensuring consistency across different parts of the implementation. Since multiple people were working on related components, there were moments where assumptions about data structures or function behavior did not fully match.

Another difficulty was debugging the Viterbi algorithm, especially when dealing with indexing, transition probabilities, and traceback logic. Small mistakes in these areas could lead to incorrect outputs, which required careful step-by-step verification.

We also encountered challenges related to integrating code from different team members. Differences in coding style and structure required additional time to standardize and ensure compatibility across the project.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Tien Nguyen
This project was more challenging compared to previous assignments because we were not provided with a notebook or step-by-step instructions to follow. Instead, we had to design our own algorithm and develop the functions from scratch. This required a deeper level of understanding and independent thinking.

Our group chose to use a class-based structure, as suggested in class. While this approach was powerful, I initially found it difficult to fully understand how the classes interacted and how to use them effectively. Without example outputs or reference implementations, deciding on the appropriate data structures for each component was also a challenge.

Despite these difficulties, the project was a valuable learning experience. It helped me gain a clearer understanding of Hidden Markov Models, especially the Viterbi algorithm and how it is used to find the most optimal path of hidden states. Overall, this project improved both my problem-solving skills and my understanding of probabilistic modeling.
# Generative AI Appendix
Anthropic. (2026). Claude (claude-sonnet-4-6) [Large language model]. https://claude.ai

AI was utilized to assist in debugging.
