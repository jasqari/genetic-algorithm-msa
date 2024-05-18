# Multiple Sequence Alignment using Genetic Algorithm
<p align="center">
 <img src="https://github.com/jasqari/genetic-algorithm-msa/assets/44480584/c95b0f64-bf6c-4f20-8a5f-42d2a5350aa5"/>
</p>

## Overview
This repository incorporates the concepts of sequence alignment and genetic algorithm to perform the task of multiple sequence alignment.
* `msa.py` is comprised of different functions useful for MSA,
* `blosum62.py` and `nuc44.py` are protein and nucleotide substitution matrices used for calculating alignment scores,
* `genetic_algorithm.py` is an implementation of the genetic algorithm appropriate for the task,
* `main.py` is the integration of all the modules.

## Requirements
```
pip3 install -r requirements.txt
```

## Usage
Take a closer look at the arguments and options involved:
```
python3 main.py -h
```

Select the type of the input sequences and specify the location of the file containing them:
```
python3 main.py protein "Sample Inputs/sample_input_protein.txt"
```

Steer the search process towards better solutions by tuning the available parameters:
```
python3 main.py dna "Sample Inputs/sample_input_dna.txt" --pop_size 3000 --num_gens 500 --epsilon 1e-5
```

Output a log of the process:
```
python3 main.py protein "Sample Inputs/sample_input_protein.txt" --log
```

Use the visualization functionality (only with IDLE shell):
```
python3 main.py protein "Sample Inputs/sample_input_protein.txt" --log --display
```

Solutions are saved to `optimal_alignment.txt`.