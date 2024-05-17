import numpy as np


def genetic_algorithm(
    _input, pop_size, num_gens, rand_generator, fitness_func, epsilon, seq_type, log
):
    """The Core function of the genetic algorithm"""
    # Create an initial population
    population = [rand_generator(_input) for _ in range(pop_size)]

    # Keep track of convergence
    last_average_score = sum([fitness_func(chrom, _input, seq_type) for chrom in population]) / len(
        population
    )
    if log:
        print("Gen 1 average score: ", last_average_score)

    # Run the main loop
    for generation in range(num_gens):
        next_population_candidates = []
        while len(next_population_candidates) < pop_size:
            children = crossover(_input, population, fitness_func, seq_type)
            mutated_children = [mutation(_input, c) for c in children]
            next_population_candidates += mutated_children
        next_population_candidates += population
        population = selection_for_survival(
            next_population_candidates, _input, pop_size, fitness_func, seq_type
        )

        average_score = sum([fitness_func(chrom, _input, seq_type) for chrom in population]) / len(
            population
        )
        if log and (generation + 1) % 5 == 0:
            print("Gen {} average score: {}".format(generation + 1, last_average_score))

        # Check for convergence
        if abs(last_average_score - average_score) < epsilon:
            if log:
                print("Gen {} average score: {}".format(generation + 1, last_average_score))
            break
        else:
            last_average_score = average_score
    return population[np.argmax([fitness_func(chroms, _input, seq_type) for chroms in population])]


def selection_for_variation(_input, population, fitness_func, seq_type, method="tournament", K=2):
    """Function to perform K-tournament selection with a chance for duplicates"""
    if method == "tournament":
        participants = [population[np.random.choice(len(population))] for _ in range(K)]
        return participants[
            np.argmax([fitness_func(participant, _input, seq_type) for participant in participants])
        ]


def crossover(_input, population, fitness_func, seq_type):
    """Function to perform vertical crossover"""
    parent1 = selection_for_variation(_input, population, fitness_func, seq_type)
    parent2 = selection_for_variation(_input, population, fitness_func, seq_type)

    # Find the valid positions
    valid_crossover_positions = []
    max_seq_length = max([len(seq) for seq in _input])
    for seq in range(len(_input)):
        seq1 = parent1[seq]
        seq2 = parent2[seq]
        valids = []
        for position in range(max_seq_length):
            if (
                len([gap_pos for gap_pos in seq1 if gap_pos <= position])
                - len([gap_pos for gap_pos in seq2 if gap_pos <= position])
                == 0
            ):
                valids.append(position)
        valid_crossover_positions.append(valids)

    # Choose the vertical points at random
    crossover_points = [np.random.choice(ps) for ps in valid_crossover_positions]

    # Generate the children
    child1 = []
    child2 = []
    for gene in range(len(parent1)):
        gap_pos1 = []
        gap_pos2 = []
        for gap_position in parent1[gene]:
            if gap_position <= crossover_points[gene]:
                gap_pos1.append(gap_position)
            else:
                gap_pos2.append(gap_position)
        for gap_position in parent2[gene]:
            if gap_position > crossover_points[gene]:
                gap_pos1.append(gap_position)
            else:
                gap_pos2.append(gap_position)
        child1.append(sorted(gap_pos1))
        child2.append(sorted(gap_pos2))
    return [child1, child2]


def mutation(_input, chromosome):
    """Function to perform the MergeSpace mutation"""
    random_seq = 0
    distant_gap_positions = []

    # Try to merge space in a sequence, if not possible, try another sequence
    for it in range(len(chromosome)):
        if len(distant_gap_positions) > 1:
            break
        distant_gap_positions = []
        random_seq = np.random.choice(len(chromosome))
        seq = chromosome[random_seq]

        # Sequences with no gaps are substituted
        while not seq:
            random_seq = np.random.choice(len(chromosome))
            seq = chromosome[random_seq]

        # Find the single gap positions and check the first and last positions
        # to handle potential errors
        if not abs(seq[0] - seq[1]) <= 1:
            distant_gap_positions.append(seq[0])
        if not abs(seq[-1] - seq[-2]) <= 1:
            distant_gap_positions.append(seq[-1])

        # Then, check the rest of the positions
        for pos in range(1, len(seq) - 1):
            if not (abs(seq[pos] - seq[pos - 1]) <= 1 or abs(seq[pos] - seq[pos + 1]) <= 1):
                distant_gap_positions.append(seq[pos])
    distant_gap_positions = sorted(distant_gap_positions)

    # Not enough possible positions for a mutation, hence, return the same
    if len(distant_gap_positions) <= 1:
        return chromosome

    # Randomly choose the positions to be merged
    pos1 = np.random.choice(distant_gap_positions)
    distant_gap_positions.remove(pos1)
    pos2 = np.random.choice(distant_gap_positions)

    # Randomly select target positions
    possible_positions = [
        i for i in range(max([len(seq) for seq in _input])) if not i in chromosome[random_seq]
    ]
    possible_pairs = []
    for i in range(len(possible_positions) - 1):
        if possible_positions[i + 1] - possible_positions[i] == 1:
            possible_pairs.append((possible_positions[i], possible_positions[i + 1]))
    target_position = possible_pairs[np.random.choice(len(possible_pairs))]

    # Mutate the chromosomes
    mutated_chromosome = [
        chromosome[i] if i != random_seq else [pos for pos in target_position] + chromosome[i]
        for i in range(len(chromosome))
    ]
    mutated_chromosome[random_seq].remove(pos1)
    mutated_chromosome[random_seq].remove(pos2)
    mutated_chromosome[random_seq] = sorted(mutated_chromosome[random_seq])
    return mutated_chromosome


def selection_for_survival(next_population_candidates, _input, pop_size, fitness_func, seq_type):
    """Function to perform selection for survival to the next evolution"""
    # The best individuals from the union of parents and children are selected
    next_population = sorted(
        next_population_candidates, key=lambda x: fitness_func(x, _input, seq_type), reverse=True
    )
    return next_population[:pop_size]
