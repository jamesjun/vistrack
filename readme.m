Data from

Signatures of value comparison in ventral striatum neurons

by
Caleb E. Strait, Brianna J. Sleezer, and Benjamin Y. Hayden

VSgamblingdata.mat contains 124-item cell array, and each cell contains the following two data structures for one of the 124 VS neurons:

vars - A [T x 11] matrix where each row is one of T trials and each of 11 columns is a variable:

1.  [1st option] probability to win
2.  [1st option] reward size
3.  [1st option] expected value
4.  [2nd option] probability to win
5.  [2nd option] reward size
6.  [2nd option] expected value
7.  Side of first (1 = L; 0 = R)
8.  Choice (Left or Right)
9.  Choice (1st or 2nd)
10. Experienced reward
11. Valid trial (1 = use for analyses; 2 = invalid trial)

psth - A [T x 750] matrix, where T is the number of trials recorded for that neuron. Each row is a trial split into 750 20ms bins, each with the number of times that neuron fired in that bin.
    