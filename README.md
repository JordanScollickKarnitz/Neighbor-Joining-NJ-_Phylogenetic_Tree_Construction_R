## Neighbor-Joining (NJ)

This implementation constructs a phylogenetic tree using the Neighbor-Joining
algorithm with Jukes–Cantor distances.

**Features**
- Explicit Q-matrix calculation
- Branch length estimation
- Tree topology reconstruction

**Time Complexity:** O(n³)  
**Space Complexity:** O(n²)

**Language:** R

Algorithmic Complexity
Let n = number of sequences.
Distance Matrix Construction

Time: O(n² · L) (L = sequence length)
Space: O(n²)

Neighbor-Joining Iterations

Each iteration:
Q-matrix computation: O(n²)
Minimum search: O(n²)
Distance matrix update: O(n)
Number of iterations: n − 2

➡️ Total NJ runtime:
𝑂(𝑛3)
O(n3)
Space Complexity
Distance matrix: O(n²)
Q-matrix: O(n²)
Tree storage: O(n)

➡️ Total space:
𝑂(𝑛2)
O(n2)

Interpretation (Bioinformatics Context):
Neighbor-Joining scales better than UPGMA/WPGMA for unequal evolutionary rates
Cubic time makes NJ suitable for small–medium datasets
Practical NJ implementations use heuristics or optimizations for large phylogenies
