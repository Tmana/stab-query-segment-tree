# stab-query-segment-tree
A heap based implementation of a segment tree, optimized to make cheap stab queries in O(log(n)) time. 

Implemented design from Lei Wang and Xiaodong Wang (2018) paper "A Simple and Space Efficient Segment Tree Implementation" - https://arxiv.org/pdf/1807.05356.pdf


## Problem:

The exercise is to write a computer program in Python3 (only using the standard library) to calculate the coverage for each of the positions in loci.csv, based on the reads in reads.csv.

### Data

#### Input - reads.csv
two columns corresponding to the
start position and length of the read.

start,length

#### Output - loci.csv

this file contains 1000 positions of interest and a blank "coverage" column. In the course of this exercise,
populate the coverage column with the number of reads overlapping that position.

position,coverage


N = # lines in reads.csv, aka segments = ~2 million
n = # of unique endpoints in N
k = # lines of loci (queries) = 1000


## Approaches:

### Naive
I first considered the performance of the naive solution, checking each point in loci against every interval in reads, which would take O(N*k). since loci is relatively short, this might not be completely terrible despite being the brute force solution.


### Simple
Then I thought about pre-processing a sorted list of segments, and binary search tree. That would be O(N*log(N)) one time for the sort, and O(N*log(N)) for each query. again, loci being so short, this might not be bad. In the case where we would have lots more than 1000 of queries (likely), it seems worth seeking a more efficient way to make queries cheap.


### Optimized
This problem seems to be a stabbing query specific sub-problem of the interval and segment tree problems. Specifically, we want to minimize the number of reads we check in order to count the 
number of overlapping segments given a integer position. 

Using a segment tree (binary tree that stores segment information in its inherent structure) should allow us to count the number of overlapping intervals in O(k + log(n)) time, where n is the number of distinct intervals (smaller than our N segments) and k is the number of retrieved intervals during a query. In this case, the height of our tree which is log(n), giving us O(2log(n)) -> O(log(n)) cost per query. We would also have to store each node object’s attributes for left and right children and counts and intervals.

However, we can do a little better by using a heap structure, we can avoid storing any structure and just navigate our tree by some “Fun” Heap Math! (jazz hands). Instead we can store the values in an arbitrary data structure for later access based on its node index. we will use an array. By doing this, we can reduce the storage cost of the tree down to a 2n-1 size array to store counts, an n size array for the points values, and each stab query down to O(log(n)), I believe is the optimal possible query cost.



## Other Considerations - pre-processing

Once the tree has been populated, it has statically encoded the intervals, since I only implemented insert and no delete, update or other heap functions. In order to avoid retraining (the longest part by far), I pickled the populated segment tree for quick loading if we want to make new stab queries not in loci.csv. Similarly for the points, I saved them to a text file to avoid reprocessing it every time. 

I also plotted the values in reads and loci to get an idea of their distributions, and found that most segments are around 150 in length (very short), there are many repeat segments, and there are large ranges that are totally uncovered.  Loci.csv however, having only 1000 values, is sparsely distributed by comparison. 

What follows are optimizations in preprocessing that could improve this specific exercise, but may not be wanted in most applications. If we really really didn’t care about maintaining all reads, and our segment tree being generalizable outside of the exact set of given loci (we probably would care in a pipeline), we could reduce our N reads and n points by removing all segments outside some margin outside the furthest left  (worst case max segment length)  and right points in loci.csv. Similarly, we could remove any segment too far away from each point in loci, to encode the smallest segment tree necessary to answer the queries in loci.csv.
