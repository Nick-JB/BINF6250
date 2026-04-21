from HMM import HMM

class ProfileHMM(HMM):
    def __init__(self, filepath, alphabet, pseudocount=1.0, gap_threshold=0.5):
        self.filepath = filepath
        self.pseudocount = pseudocount
        self.gap_threshold = gap_threshold

        self.headers, self.msa = self.read_fasta_msa(filepath)
        self.match_columns = self.get_match_columns(self.msa, gap_threshold)
        self.states, self.num_match_states = self.build_profile_states(self.match_columns)
        self.allowed_transitions = self.build_allowed_transitions(self.num_match_states)

        init_counts, trans_counts, emit_counts = self.count_transitions_and_emissions(
            self.msa,
            self.match_columns,
            self.states,
            alphabet,
            self.allowed_transitions,
            pseudocount
        )

        init_probs, trans_probs, emit_probs = self.normalize_counts(
            init_counts,
            trans_counts,
            emit_counts
        )

        background_probs = self.get_background_distribution(self.msa, alphabet, pseudocount)
        emit_probs = self.apply_background_to_insertions(emit_probs, self.states, background_probs)

        super().__init__(
            alphabet=alphabet,
            hidden_states=self.states,
            init_probs=init_probs,
            trans_probs=trans_probs,
            emit_probs=emit_probs
        )

    def read_fasta_msa(self, filepath):
        """
        Read an aligned FASTA file and return the sequences.

        This function assumes the FASTA file is a multiple sequence alignment.
        So all sequences must have the same length, including gaps.
        """
        headers = []
        sequences = []
        current_header = None
        current_seq = []

        with open(filepath, "r") as infile:
            for line in infile:
                line = line.strip()

                # skip empty lines
                if not line:
                    continue

                # new FASTA header
                if line.startswith(">"):
                    # save the previous sequence before starting a new one
                    if current_header is not None:
                        headers.append(current_header)
                        sequences.append("".join(current_seq).upper())
                    current_header = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
        # Save the final sequence after the loop ends:
        if current_header is not None:
            headers.append(current_header)
            sequences.append("".join(current_seq).upper())

        if not sequences:
            raise ValueError("No sequences were found.")

        # All aligned sequences must have the same length
        seq_length = len(sequences[0])
        for seq in sequences:
            if len(seq) != seq_length:
                raise ValueError("All sequences must have the same length.")

        return headers, sequences

    def get_match_columns(self, msa, gap_threshold=0.5):
        """Determine which columns is the MSA are match columns"""
        num_sequences = len(msa)
        aligned_length = len(msa[0])

        match_columns = []

        # iterate over each column index
        for col_idx in range(aligned_length):
            gap_count = 0

            # count how many sequences have a gap in this column
            for seq in msa:
                if seq[col_idx] == '-':
                    gap_count += 1

            gap_fraction = gap_count / num_sequences

            # apply threshold rule
            if gap_fraction < gap_threshold:
                match_columns.append(True)  # match column
            else:
                match_columns.append(False)  # insertion column

        return match_columns

    def build_profile_states(self, match_columns):
        """
        Build the ordered list of profile HMM states.
        """
        num_match_states = sum(match_columns)

        states = []
        states.append("B")
        states.append("I0")

        for k in range(1, num_match_states + 1):
            states.append("M" + str(k))
            states.append("D" + str(k))
            states.append("I" + str(k))

        states.append("E")

        return states, num_match_states

    def build_allowed_transitions(self, num_match_states):
        """
        Build the allowed transition structure for a profile HMM
        """
        allowed = {}

        # begin
        allowed['B'] = ['M1', 'I0', 'D1']

        # I0
        allowed['I0'] = ['I0', 'M1', 'D1']

        # internal positions
        for k in range(1, num_match_states + 1):
            m = "M" + str(k)
            d = "D" + str(k)
            i = "I" + str(k)

            if k < num_match_states:
                next_m = "M" + str(k + 1)
                next_d = "D" + str(k + 1)

                allowed[m] = [next_m, i, next_d]
                allowed[d] = [next_m, i, next_d]
                allowed[i] = [next_m, i, next_d]
            else:
                # final position goes to E
                allowed[m] = [i, "E"]
                allowed[i] = [i, "E"]
                allowed[d] = [i, "E"]
        # end
        allowed['E'] = ["E"]

        return allowed

    def get_background_distribution(self, msa, alphabet, pseudocount=1.0):
        """
        Estimate a background residue distribution from the MSA

        This background distribution will be used for all insertion states including the initial insertion
        Gaps are ignored because '-' is not an emitted residue
        """

        # make sure alphabet is easy to iterate over
        alphabet = list(alphabet)

        # start each residue count with a pseudocount
        counts = {}
        for residue in alphabet:
            counts[residue] = pseudocount

        # count residues across the whole alignment, ignoring gap
        for seq in msa:
            for char in seq:
                if char in counts:
                    counts[char] += 1

        # convert counts to probabilities
        total = sum(counts.values())
        bg_probs = {}
        for residue in alphabet:
            bg_probs[residue] = counts[residue] / total

        return bg_probs

    def apply_background_to_insertions(self, emit_probs, states, background_probs):
        """
        Force all insertion states to share the same background emission distribution.
        """
        for state in states:
            if state.startswith("I"):
                emit_probs[state] = {}

                for residue in background_probs:
                    emit_probs[state][residue] = background_probs[residue]

        return emit_probs

    def get_state(self, symbol, is_match_column, match_index):
        """
        Determine the HMM state for a given symbol and column type.
        """
        if is_match_column:
            match_index += 1

            if symbol == "-":
                return "D" + str(match_index), match_index
            else:
                return "M" + str(match_index), match_index

        else:
            if symbol != "-":
                return "I" + str(match_index), match_index
            else:
                return None, match_index

    def aligned_to_state_path(self, aligned_seq, match_columns):
        """
        Convert one aligned sequence into a profile HMM state path
        """
        path = []
        path.append("B")
        match_index = 0

        for col_idx in range(len(aligned_seq)):
            symbol = aligned_seq[col_idx]
            is_match_column = match_columns[col_idx]

            state, match_index = self.get_state(symbol, is_match_column, match_index)

            if state is not None:
                path.append(state)

        # add ending for path
        path.append("E")
        return path

    def count_transitions_and_emissions(self, msa, match_columns, states, alphabet, allowed_transitions, pseudocount=1.0):
        """
        Count initial states, transitions, and emissions from an aligned training MSA"
        """
        alphabet = list(alphabet)

        # 1. Initialize initial state-counts
        init_counts = {}
        for state in states:
            init_counts[state] = 0.0
        init_counts["B"] = 1.0

        # 2. Initialize transition counts
        trans_counts = {}
        for from_state in states:
            trans_counts[from_state] = {}
            for to_state in states:
                trans_counts[from_state][to_state] = 0.0

        # add pseudocounts only to allowed transitions
        for from_state in allowed_transitions:
            for to_state in allowed_transitions[from_state]:
                trans_counts[from_state][to_state] = pseudocount

        # 3. Initialize emission counts
        emit_counts = {}
        for state in states:
            emit_counts[state] = {}
            for residue in alphabet:
                emit_counts[state][residue] = pseudocount

        # 4. Process each aligned sequence
        for aligned_seq in msa:
            path = self.aligned_to_state_path(aligned_seq, match_columns)

            # . count transitions from the path
            for i in range(len(path) - 1):
                from_state = path[i]
                to_state = path[i + 1]
                if to_state in self.allowed_transitions[from_state]:
                    trans_counts[from_state][to_state] += 1

            # 6 count emission by waking through the aligned sequence again
            match_index = 0

            for col_idx in range(len(aligned_seq)):
                symbol = aligned_seq[col_idx]
                is_match_column = match_columns[col_idx]

                state, match_index = self.get_state(symbol, is_match_column, match_index)

                if state is not None:
                    if state.startswith("D"):
                        continue
                    emit_counts[state][symbol] += 1

        return init_counts, trans_counts, emit_counts

    def normalize_counts(self, init_counts, trans_counts, emit_counts):
        """"
        Convert count tables into probabilites tables
        """
        # 1. normalize initial counts
        init_probs = {}

        init_total = 0.0

        for state in init_counts:
            init_total += init_counts[state]

        for state in init_counts:
            init_probs[state] = init_counts[state] / init_total

        # 2 normalize transition counts
        trans_probs = {}

        for from_state in trans_counts:
            trans_probs[from_state] = {}

            row_total = 0.0
            for to_state in trans_counts[from_state]:
                row_total += trans_counts[from_state][to_state]

            for to_state in trans_counts[from_state]:
                if row_total == 0:
                    trans_probs[from_state][to_state] = 0.0
                else:
                    trans_probs[from_state][to_state] = trans_counts[from_state][to_state] / row_total

        # 3 normalize emission counts
        emit_probs = {}

        for state in emit_counts:
            emit_probs[state] = {}
            row_total = 0.0
            for residue in emit_counts[state]:
                row_total += emit_counts[state][residue]

            for residue in emit_counts[state]:
                if row_total == 0:
                    emit_probs[state][residue] = 0.0
                else:
                    emit_probs[state][residue] = emit_counts[state][residue] / row_total

        return init_probs, trans_probs, emit_probs

def main():
    train_file = "C:/Users/botto/Documents/NU_26_Spring/BINF6250/11BINF6250/BINF6250/project11/data/phmm_train_motif1.fasta"
    test_file = "C:/Users/botto/Documents/NU_26_Spring/BINF6250/11BINF6250/BINF6250/project11/data/phmm_test_sequences.fasta"
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    model = ProfileHMM(train_file, alphabet)

    print("States:")
    print(model.hidden_states)
    print()

    print("Match Columns:")
    print(model.match_columns)
    print()

    print("Initial Probabilities:")
    print(model.init_probs)
    print()

    print("Transition Probabilities:")
    print(model.trans_probs)
    print()

    print("Emission Probabilities:")
    print(model.emit_probs)
    print()

    test_headers, test_sequences = model.read_fasta_msa(test_file)

    for i in range(len(test_sequences)):
        seq = test_sequences[i].replace("-", "")
        name = test_headers[i]

        print("Sequence:", name)
        print("Ungapped Sequence:", seq)

        forward_prob, _ = model.forward(seq)
        print("Forward Probability:", forward_prob)

        viterbi_path = model.viterbi(seq)
        print("Viterbi Path:", viterbi_path)
        print("-" * 50)


if __name__ == "__main__":
    main()
