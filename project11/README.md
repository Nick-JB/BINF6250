# Introduction
This project focuses on implementing a Profile HMM for protein motif detection by extending the provided HMM.py framework. Starting from a multiple sequence alignment (MSA) in FASTA format, the model is constructed by classifying alignment columns into match or insertion positions and building a structured, left-to-right state topology consisting of match, insertion, and deletion states, along with begin and end states. Using this structure, the model estimates transition and emission probabilities from the aligned sequences. The profile HMM leverages existing algorithms, such as Forward and Viterbi, to analyze new, unaligned protein sequences. This enables both probabilistic motif detection and sequence-to-profile alignment, demonstrating how probabilistic models can be applied to biological sequence analysis and pattern recognition.


# Pseudocode
```
read_fasta(filepath):
Parses an MSA FASTA file into a dictionary mapping headers to sequences.
Initialize an empty dictionary, a current header string, and a current sequence string

Read file line by line
If line starts with >, treat it as a new header
If a previous header exists, store it with its sequence in the dictionary
Reset sequence accumulation for the new entry

Otherwise append cleaned line to current sequence

After reading all lines, store the final sequence
Return dictionary with all sequences uppercased

ProfileHMM class:

Extends a base HMM while fixing multi-character state handling

Inherits forward, backward, forward_backward, and baum_welch without changes

Overrides hidden state and probability setters to store states as a list instead of a string
This prevents corruption of multi-character states like M1 or D3

Maintains all original probability validation rules

from_msa(filepath):

Builds a profile HMM from a multiple sequence alignment

Column classification:

Compute gap fraction per column
Columns below threshold become match columns
Others become insert columns
Number of match columns defines profile length L

State structure:

Start with B
For each position add match, insert, and delete states
End with IL and E
Delete and terminal states are silent

Sequence processing:

Walk each sequence through alignment
Match columns emit from match states
Gaps in match columns create delete states
Insert columns emit from insert states or update background counts
Build state path ending in E
Record transitions between states

Emissions:

Match states use column specific counts
Insert states use their own counts or fallback background distribution
Silent states use background distribution as placeholder

Transitions:

Use Plan7 topology
Normalize observed transitions with pseudocounts
Ensure valid row stochastic matrix

Initialization:

Only B has probability 1
All others are 0

Output:

Return fully constructed HMM

_build_allowed_transitions(L)

Defines Plan7 topology:

B connects to M1, I0, D1

For each position k
Match states connect forward, to insert, and delete states
Insert states loop and connect forward
Delete states skip emission and move forward

End state E is absorbing

viterbi(sequence):

Implements corrected Viterbi decoding for multi character states

Uses log probabilities with safe handling of zero values

Initialize scores from start probabilities and emissions

For each sequence position compute best predecessor state
Store traceback pointers instead of string paths

At end select best final state
Reconstruct path by backtracking pointers

Return best score and state sequence

score_sequence(sequence):

Removes gaps
Runs forward algorithm
Returns total probability

log_odds_score(sequence):

Compares model likelihood to uniform background model

Computes forward probability
Computes random expectation based on alphabet size
Returns log ratio

Positive score indicates motif-like sequence

export_params:

Serializes model parameters using existing to_json function
Prints initial, transition, and emission probabilities 
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Tien Nguyen
I struggled to understand the provided HMM.py, especially how it was structured and how different components interacted. One of the main challenges was that hidden_states were represented as a string instead of a list, which made it less intuitive to work with, particularly for multi-character states like M1 or D2. I decided to modify that part of HMM.py to store hidden_states as a list, which made the implementation clearer and easier to debug.

I also found it difficult to translate the concept of a Profile HMM into code. Mapping alignment columns to match, insertion, and deletion states, and enforcing the correct transition structure required careful thinking and debugging.

Overall, I feel relieved that we were able to figure it out to this point. We were not able to implement Baum–Welch due to time constraints, but I gained a much better understanding of how HMMs are constructed and applied in practice.

# Generative AI Appendix
ChatGPT was used to explain the given HMM.py.
