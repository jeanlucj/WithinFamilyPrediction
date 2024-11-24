# WithinFamilyPrediction
Within a general population, find segments that are IBD to the two parents to
create a kernel that (hopefully) will predict progeny within that family better
than a general population wide kernel

assessIndMrkSimilarityToPar.R
Functions
calcMrkSimIndPar
Count the number of sequential markers that have the same states as the parent
Pass both the individual and the parent as Formal class pop
Return a list with one, two or four vectors depending on whether both, one
or neither are inbred
Each vector contains the number of markers in the sequence, and its start position

findPosMaxSeqSame
Find the longest number of contiguous markers that are IBS

calcProbSameSame
Moving from left to right along a chromosome, if two haplotypes are the same
at the first locus, what probability are they the same at the second locus?

findPosMinProbSame
Find the sequence of contiguous markers that are IIS with the lowest probability
Return two values: -log10 of the probability, the sequence start position

calcAllPosProbSame
This function calculates the probabilities of all streaks of contiguous
markers that are IIS. It returns the position of the streak, its length, and
its probability for all such streaks on a chromosome.

partitionMrkForBiparPred
Function to split markers into one matrix if they appear to be IBD with either
parent of a biparental and a separate matrix if they they are not.


A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr
