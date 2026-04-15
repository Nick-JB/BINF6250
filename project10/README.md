# Introduction
This project was completed as a group and tasked us with implementing the Baum-Welch algorithm for parameter estimation in Hidden Markov Models (HMMs). The Baum-Welch algorithm is an Expectation-Maximization (EM) approach that allows an HMM to learn its own parameters — initial state probabilities, transition probabilities, and emission probabilities — directly from unlabeled observation sequences, without needing any ground-truth state labels.

Our implementation takes in four data structures (observation sequences, initial probabilities, transition probabilities, and emission probabilities) and uses them to iteratively refine the model until it converges on a locally optimal set of parameters. We made no assumptions about the number of hidden states, the size of the observation alphabet, or the domain of the data, keeping our solution fully general.

To structure our work, we used Object-Oriented Programming (OOP). Having already implemented Viterbi, Forward, Backward, and Posterior decoding in our existing HMM class, the Baum-Welch algorithm was a natural extension — the E-step draws directly on our Forward and Backward methods to compute expected counts, and the M-step uses those counts to re-estimate parameters, completing the full unsupervised learning pipeline.
# Pseudocode

```
## Pseudocode

### `_safe_log(x)`
- Input: probability value `x`
- If `x <= 0`
  - Return `-∞`
- Else
  - Return `log(x)`

---

### `_synchronize_emission()`
- For each state in the HMM:
  - For each emission in that state:
    - Add emission to global emission set

- For each emission in global emission set:
  - For each state:
    - If state does not contain emission:
      - Add emission to state
      - Set emission probability to `0`

---

### `_expectation_values(observations)`
- Convert observations to list
- Let `n_states = number of hidden states`
- Let `T = length of observations`

- Compute `forward_mat`
- Compute `backward_mat`

- Compute sequence log-likelihood:
  - `log_prob = log-sum-exp(last column of forward_mat)`

- Compute gamma:
  - `gamma = exp(forward_mat + backward_mat - log_prob)`

- Initialize `xi` as zeros of shape `(T-1, n_states, n_states)`

- For each time step `t` from `0` to `T-2`:
  - For each state `i`:
    - For each state `j`:
      - Compute:
        - forward probability at `(i, t)`
        - transition probability `i → j`
        - emission probability of next symbol from state `j`
        - backward probability at `(j, t+1)`
      - Combine into log-xi
      - Convert to normal probability
      - Store in `xi[t, i, j]`

- Return:
  - forward matrix
  - backward matrix
  - log-likelihood
  - gamma
  - xi

---

### `sequence_log_likelihood(observations)`
- Call `_expectation_values(observations)`
- Extract `log_prob`
- Return `log_prob`

---

### `gamma_matrix(observations)`
- Call `_expectation_values(observations)`
- Extract `gamma`
- Return `gamma`

---

### `xi_tensor(observations)`
- Call `_expectation_values(observations)`
- Extract `xi`
- Return `xi`

---

### `posterior_decoding(observations)`
- Compute `gamma = gamma_matrix(observations)`

- Convert gamma to log-space:
  - If `gamma > 0` → `log(gamma)`
  - Else → `-∞`

- For each position:
  - Select state with highest posterior probability

- Convert indices to state objects

- Return:
  - posterior matrix
  - state path

---

### `baum_welch(sequences, max_iter, tol, pseudocount)`

- Let `n_states = number of hidden states`
- Create sorted list of emission symbols
- Create dictionary mapping symbol → index

- Initialize empty list `history`
- Set `prev_loglik = -∞`

- For each iteration up to `max_iter`:

  - Initialize:
    - `beta_counts = 0`
    - `trans_counts = 0`
    - `trans_denom = 0`
    - `emit_counts = 0`
    - `emit_denom = 0`
    - `total_loglik = 0`

  - For each sequence:
    - Convert sequence to list
    - Compute `(log_prob, gamma, xi)`

    - Add `log_prob` to total log-likelihood

    - Update initial state counts:
      - `beta_counts += gamma[:, 0]`

    - Update transition counts:
      - `trans_counts += sum(xi over time)`
      - `trans_denom += sum(gamma except last column)`

    - For each position `t`:
      - Identify observed symbol
      - Find its index
      - Add `gamma[:, t]` to emission counts

    - Update emission denominator:
      - `emit_denom += sum(gamma)`

  - Update initial probabilities:
    - Normalize `beta_counts` with pseudocount

  - Update transition probabilities:
    - Normalize `trans_counts` with pseudocount

  - Update emission probabilities:
    - Normalize `emit_counts` with pseudocount

  - Append `total_loglik` to history

  - If improvement < `tol`:
    - Stop training

  - Update `prev_loglik`

- Return `history`
```

# Successes
One of the major successes of this project was successfully integrating the Baum-Welch algorithm into our existing HMM framework. Because we had already implemented the Forward and Backward algorithms, we were able to directly integrate these components in the E-step, which reinforced our understanding of how these algorithms connect within the broader HMM pipeline. Another key success was our thorough learning process, which helped us learn from one another and set a solid foundation of understanding for us to work from. Additionally, despite our tight schedules as the end of the semester gets hectic, we were able to maintain frequent communication to keep each other up to date on any progress we made. Most importantly, we were able to effectively help each other learn and cross-validate each other’s logic. 

# Struggles
One of the primary challenges we faced was understanding the underlying concepts before implementation. We found that we needed to spend significant time working through the theory to fully grasp how the Forward and Backward algorithms contribute to expected counts. Rather than rushing into coding, we deliberately slowed down to ensure a solid conceptual foundation, which ultimately improved the correctness of our implementation but extended development time. Translating our concepts from formula to code also proved to be a confusing endeavor at first, but through our debugging, we were able to make sense of what needed to happen and how. Once we got over the conceptual hurdle, the code and testing were not as difficult.

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

## Shameem Shahib
This project helped solidify my understanding of HMMs and how the different algorithms work together. Implementing Baum-Welch made the connection between probability theory and practical machine learning much clearer, and I found it to be enriching to work alongside my group to reinforce my learning. Collaborating with the team was valuable, as it exposed me to different problem-solving approaches and helped clarify challenging concepts through discussion. Despite the initial confusion and conceptual roadblocks, I am happy with the final result we were able to come up with.

# Generative AI Appendix
As per the syllabus
