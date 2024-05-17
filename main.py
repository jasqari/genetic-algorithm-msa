import argparse
from msa import *
from genetic_algorithm import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "seq_type",
        type=str,
        choices=["protein", "dna"],
        help="Type of the sequences (protein or dna)",
    )
    parser.add_argument(
        "seqs_file_path",
        type=str,
        help="Path to the text file that contains 2 or more protein or dna sequences",
    )
    parser.add_argument(
        "--pop_size",
        type=int,
        default=3000,
        help="Population size of each generation of the evolutionary algorithm",
    )
    parser.add_argument(
        "--num_gens",
        type=int,
        default=500,
        help="Number of generations as the stopping criterion of the iterative algorithm",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1e-5,
        help="Tolerance value as the stopping criterion of the iterative algorithm",
    )
    parser.add_argument(
        "--log",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Log the iterative optimization process or not",
    )
    parser.add_argument(
        "--display",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Display the converged solution on a mesh grid or not",
    )
    args = parser.parse_args()

    seqs = load_sequences(args.seqs_file_path)

    converged_solution = genetic_algorithm(
        seqs,
        args.pop_size,
        args.num_gens,
        random_encoding,
        sum_of_pairs_score,
        args.epsilon,
        args.seq_type,
        args.log,
    )
    with open("optimal_alignment.txt", "w") as file:
        for seq in decode(converged_solution, seqs):
            file.write("{}\n".format("".join(seq)))

    if args.display:
        display_msa(decode(converged_solution, seqs), args.seq_type)
