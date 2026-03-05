from pna_calc.tm import pna_tm
from pna_calc.tm import dna_nn_tm
from pna_calc.metrics import sequence_metrics

seq = "CGTTGCATAACG"

print("DNAnn Tm:", dna_nn_tm(seq))
print("PNA Tm:", pna_tm(seq))
print(sequence_metrics(seq))
