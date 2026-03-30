# Introduction
Hidden Markov Models (HMMs) can provide probabilistic characterization for a sequence of observations that would otherwise be hard to parse.
For example predicting the seasons of the year based on weather patterns. For this project we will be using Viterbi's algorithm to select the most likely hidden states underlying a sequence of observations.

# Pseudocode
Put pseudocode in this box:

```
class State(name, emissions)
    give state a name
    give state dict of emit: probability

class HMM(emissions, states)
    create empty HMM with list of emissions
    give empty list of state or build from exisiting list

    def add_state(state)
        add a state to HMM
        update state emission dict to include emissions missing from emissions list
        make sure all emissions in state are in HMM emissions list
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
Anthropic. (2026). Claude (claude-sonnet-4-6) [Large language model]. https://claude.ai

AI was utilized to assist in debugging.
