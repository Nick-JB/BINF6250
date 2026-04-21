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
We were able to understand the provided HMM.py framework and develop a clear plan for extending it into a working Profile HMM. Interpreting the model and its required steps took longer than expected, but once we established a solid foundation through well-structured pseudocode, the implementation itself was relatively straightforward. We ensured that the model correctly classifies alignment columns into match and insert states and constructs a full Plan7-style topology, including match, insertion, deletion, begin, and end states derived directly from the input alignments. We were able to compute forward probabilities and generate Viterbi paths that pass through match states for sequences consistent with the training motif. We also resolved a key compatibility issue between the Profile HMM’s multi-character state names and BaseHMM’s original string-based state representation. By storing hidden states as a list instead of a single concatenated string, we preserved compatibility with the existing forward, backward, and forward-backward algorithms without needing to modify those methods. Most importantly, despite our hectic schedules with the end of the semester, we kept each other updated with progress and took time to explain any confusion to make sure we were all on the same page.

# Struggles
Prior to any scripting, one of the biggest challenges we faced was with the foundational pseudocode.  Mapping alignment columns to match, insert, and delete states, handling gaps appropriately, and constructing emission distributions required careful attention to detail. Understanding the algorithm and translating our ideas from concepts to pseudocode took us a lot longer than we had planned, but it also helped make the script writing much faster later on once we got to that step. Reading through HMM.py also proved to be quite time-consuming, and we had to go back several times to make sure our profile HMM didn't cause any errors with the logic in the baseHMM. One of the first challenges we faced with this logic stemmed from a mismatch between the design of BaseHMM and the requirements of a Profile HMM. BaseHMM represents hidden states as a single concatenated string and iterates over it character by character. While this works for single-character state names, it breaks for multi-character states such as M1 or D3, since a string is interpreted as individual characters rather than complete state labels. Once we got past the initial hurdles with our pseudocode and logic, the rest was not as complicated. However, with the time constraints and scheduling conflicts, we had to allocate most of our collaboration to be asynchronous.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Tien Nguyen
I struggled to understand the provided HMM.py, especially how it was structured and how different components interacted. One of the main challenges was that hidden_states were represented as a string instead of a list, which made it less intuitive to work with, particularly for multi-character states like M1 or D2. I decided to modify that part of HMM.py to store hidden_states as a list, which made the implementation clearer and easier to debug.

I also found it difficult to translate the concept of a Profile HMM into code. Mapping alignment columns to match, insertion, and deletion states, and enforcing the correct transition structure required careful thinking and debugging.

Overall, I feel relieved that we were able to figure it out to this point. We were not able to implement Baum–Welch due to time constraints, but I gained a much better understanding of how HMMs are constructed and applied in practice.

## Shameem Shahib
I found reading through HMM.py and understanding its logic to be the most time-consuming part, with the pseudocode being the second most time-consuming. I realized early on that it would not be as simple as just using some methods and finishing the implementation, and a lot of time was put into making sure I actually understood how a profile HMM works and what parts of the baseHMM I could use and what needed to be modified. Thankfully, talking with my group helped clear up some of my confusion, and once we met to finalize the pseudocode, the scripting did not take as long as the beginning steps did. I personally found the first half working through the pseudocode to be the most enriching since it helped me grasp what makes a profile HMM different to a regular HMM. Overall, I am glad with the progress we made and the results we were able to come up with.

# Generative AI Appendix
ChatGPT was used to explain the given HMM.py.
