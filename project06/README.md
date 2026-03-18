# Introduction
This project aim was to create an algorithm for implementing the neighbor joining method to build a phylogenetic tree. FASTA files are read in, Smith-Waterman scores are calculated and converted to distance scores, and a Newick string is produced to represent relatedness. The program then produces a graphical representation of the unrooted tree.

# Pseudocode
```
read_fasta(filename) -> Dict
    read in fasta file and return dictionary of 
    headers (keys) and sequences (values)

cal_score(matrix, seq1, seq1, i, j, match, mismatch, gap)
    calculate score of position in SW matrix
    (taken from assignment 5)

smith_waterman(seq1, seq2, match, mismatch, gap)
    create smith-waterman alignment matrix
    (adapted from assignment 5)

    select max score and divide by max possible score, then subtract from
    one to create distance score

    return normalized distance

build_distance_matrix(sequences from read_fasta)
    use smith waterman to find distance between each sequence from the fasta
    
    initialize empty distance matrix, shape = len(seqs) x len(seqs)

    create a sequence id list using sequences key

    Loop through sequences dict to populate distance matrix

    return distance matrix

build_Q
    Uses the raw distance matrix from S-W to create Q matrix
    calculate Q distance between each node in raw distance matrix using Q formula 
    Q(i,j) = (n − 2) · d(i,j) − Σd(i) − Σd(j)

    initialize empty q matrix, with shape = distance matrix shape
    initialize empty row sums list

    collect row sums from distance matrix

    calculate Q value for each node pair and populate Q matrix

    return Q matrix

Node (Class)
    init(name, children, distance to parent)
    name is node name
    children is list of child nodes
    distance to parent is branch length to parent internal node

calculate_distances(dij, ri, rj, n)
    Calculates the distance from node i and j to their parent internal node

    dij is distance from node i to node j

    ri is sum of node distances from node i and the same is true for rj with respect to node j

    n is the number of nodes in the distance matrix

    di = (dij + (ri - rj) / (n - 2)) / 2
    dj = dij - di

create_newick(top)
    Uses Node pointed to as top for start of recursion to dig to leaves and report branch lengths from the ground up

    base case = node.children is None
        return node name and distance to parent

    join child nodes

    set self name outside tuple of children

    add distance to parent if not top

    return newick string

neighbor_joining(distance matrix, seq ids)

    initialize empty list to hold nodes

    add leaf nodes to node list using sequence ids

    Iteration: Condense 
        Base case: 2 nodes left
            choose 1 to be parent of other and return pointer to top node

        create q matrix from distance matrix

        use q matrix to select neigbor nodes in distance matrix

        calculate distance to closest internal node

        remove child nodes from raw distance matrix

        create new node k and calculate distances to each remaining node

    create newick string from tree via recursion and return

plot_tree(newick tree)
    Plot phylo tree from newick string
```

# Successes
We met several times for brainstorming and group programming, and ingestigated both iterative and recursive methods which lent us a greater understanding of the algorithm. Our first main success came with producing a raw distance matrix and creating a Q-matrix from it. The next big step was successfully condensing and trimming the nodes from the distance matrix. Finally, we determined how to build the tree structure from the information trimmed during each iteration of the tree building.

# Struggles
Our first stumbling block was that we noticed that the lecture didn't mention Q-matrices which were mentioned in the notebook as well as external resources about neighbor joining, so we spent a lot of time trying to understand how Q-matrices interacted with the methods described in the lecture. 
We also spent a good amount of time trying to decide on Node object structure and how best to point between child and parent nodes, and how that would affect our ability to build the Newick string.

# Personal Reflections
## Nicholas Bottomley 
This project was more difficult than I expected it to be. Originally, we tried to implement a tree building method recursively utilizing balding and trimming, but found ourselves confused about how to construct the tree once we had reached a 2x2 matrix. After multiple group meetings we were able to implement an iterative logic for condensing the distance matrix. Overall, this project was not too difficult, but there were definitely a couple spots where it felt like the answer was right there, but I was just short of it. It was extremely rewarding when we were able to finally produce a phyogenetic tree however.

## Hongyuan Deng 
Reflecting on this module, our team spent a significant amount of time wrestling with the algorithmic design of the Neighbor-Joining (NJ) function, specifically deciding between a recursive and an iterative approach.However, through deep diving into the implementation,we utilized an iterative "string-building" technique. By dynamically updating the numpy distance matrix and formatting the merged nodes directly into Newick substrings on the fly, we bypassed the deep recursion limits in Python.

## Victoria Van Berlo
This project seemed simple, but proved rather confounding the deeper we went. Recursion is a weakness of mine, but my group members helped me to understand and having both methods solidified things for me.

# Generative AI Appendix
Anthropic. (2026). Claude (claude-sonnet-4-6) [Large language model]. https://claude.ai

AI was used to assist with a Q matrix equation formula, debugging, and assistance with numpy and pandas syntax.
