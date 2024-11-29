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
24 Nov 2024 I think partitionMrkForBiparPred is the full deal, but I don't have
an example using it yet.  Something like

1. Make founder population.  Say 1000 inbreds.  
2. Split off 200 for test.  
3. X number of times.  
4. Choose two parents and make progeny.  
5. With remaining population determine IBD and use two kernel prediction.  
6. Try different thresholds and different weights for the two kernels.  
7. Choose the combo of thresholds and weights that gave the best accuracy.  
8. Try that combo on some biparentals from the test population.  

genomicPredictionFunction.R  
genPred  
Use RRBLUP to run genomic prediction  

genPredMultiChr  
it's assumed that you enter with a list of training populations that has as many
training populations as there are chromosomes  

genPredPartitionedSnps  
There are two SNP matrices, one for portions of the genome that are presumed IBD
with the parents and one for portions that are presumed not IBD This function
returns the SNP effects from the ibd matrix

Other thoughts
Decreasing the weight put on individuals as the error associated with their true
degree of relatedness increases would make sense.
Hill, W.G., and B.S. Weir. 2011. Variation in actual relationship as a consequence of Mendelian sampling and linkage. Genet. Res.  93(1): 47–64. doi: 10.1017/S0016672310000480.

A [workflowr][] project.
[workflowr]: https://github.com/workflowr/workflowr
