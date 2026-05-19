import msprime
import argparse
import random
import numpy as np  # Import numpy for growth rate calculations

def run_simulation(sample_size, seq_length, Ne,N1,N2,T1,T2, mut_rate, rec_rate, output_vcf, seed):
    """
    Runs a single neutral simulation with msprime and outputs a VCF file.
    This version uses the modern msprime API (>= 1.0) and includes a
    recent contraction followed by an ancient expansion.

    :param sample_size: Number of diploid individuals to simulate.
    :param seq_length: The length of the sequence (chromosome) in base pairs.
    :param Ne: The effective population size *at the present day*.
    :param mut_rate: The mutation rate per base pair per generation.
    :param rec_rate: The recombination rate per base pair per generation.
    :param output_vcf: Path to the output VCF file.
    :param seed: A random seed for reproducibility.
    """
    print(f"Starting simulation with seed {seed}:")
    print(f"  Sample size: {sample_size} individuals ({sample_size * 2} haplotypes)")
    print(f"  Sequence length: {seq_length} bp")
    print(f"  Effective population size (Ne at present): {Ne}")
    print(f"  Mutation rate: {mut_rate}")
    print(f"  Recombination rate: {rec_rate}")

    N0_now = Ne  # Present-day effective population size
    N1_3k = N1  # Size at 3,000 generations ago
    N2_100k = N2  # Size at 100,000 generations ago

    t_recent_end = T1
    t_ancient_end = T2

    alpha_0_3k   = np.log(N0_now / N1_3k) / t_recent_end
    alpha_3k_100k = -np.log(N2_100k / N1_3k) / (t_ancient_end-t_recent_end)

    demography = msprime.Demography()
    
    demography.add_population(name="POP", initial_size=Ne, growth_rate=alpha_0_3k)

    # # 4. Add demographic events (time moves backward)
    demography.add_population_parameters_change(
        time=t_recent_end,
        initial_size=N1_3k,
        growth_rate=alpha_3k_100k,
        population="POP"
    )
    demography.add_population_parameters_change(
        time=t_ancient_end,
        initial_size=N2_100k,
        growth_rate=0,  # Stop growth
        population="POP"
    )
    


    # # Event 2: At 20,000 generations ago, the population enters the second growth phase.
    # #          At this point, size is 0.2 * Ne and starts growing at rate_ancient.
    # demography.add_population_parameters_change(
    #     time=t_recent_end,
    #     initial_size=size_at_t_recent,
    #     growth_rate=rate_ancient,
    #     population="POP"
    # )

    
    print("\n--- Demography Model (using modern API) ---")
    print(f"  - 0 -> {t_recent_end} gens ago: Size changes from {N0_now} to {N1_3k}")
    print(f"  - {t_recent_end} -> {t_ancient_end} gens ago: Size changes from {N1_3k} to {N2_100k}")
    print(f"  - > {t_ancient_end} gens ago: Ne = {N2_100k}")
    print("-------------------------------------------\n")

    ts = msprime.sim_ancestry(
        samples={"POP": sample_size},
        demography=demography,
        sequence_length=seq_length,
        recombination_rate=rec_rate,
        random_seed=seed,
        model="hudson"
    )

    mts = msprime.sim_mutations(ts, rate=mut_rate, random_seed=seed)

    with open(output_vcf, "w") as vcf_file:
        mts.write_vcf(vcf_file)

    print(f"Simulation complete. VCF written to {output_vcf}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a neutral simulation with msprime.")
    parser.add_argument("--sample-size", type=int, required=True, help="Number of diploid individuals.")
    parser.add_argument("--seq-length", type=int, required=True, help="Length of the chromosome in bp.")
    parser.add_argument("--Ne", type=float, default=10000, help="Effective population size.")
    parser.add_argument("--N1", type=float, default=10000, help="Effective population size.")
    parser.add_argument("--N2", type=float, default=10000, help="Effective population size.")
    parser.add_argument("--T1", type=float, default=10000, help="Effective population size.")
    parser.add_argument("--T2", type=float, default=10000, help="Effective population size.")
    parser.add_argument("--mut-rate", type=float, default=1e-8, help="Mutation rate (e.g., 1e-8).")
    parser.add_argument("--rec-rate", type=float, default=1e-8, help="Recombination rate (e.g., 1e-8).")
    parser.add_argument("--output-vcf", type=str, required=True, help="Path for the output VCF file.")
    parser.add_argument("--seed", type=int, help="Random seed for the simulation.", default=random.randint(1, 2**32 - 1))

    args = parser.parse_args()

    run_simulation(
        sample_size=args.sample_size,
        seq_length=args.seq_length,
        Ne=args.Ne,
        N1=args.N1,
        N2=args.N2,
        T1=args.T1,
        T2=args.T2,
        mut_rate=args.mut_rate,
        rec_rate=args.rec_rate,
        output_vcf=args.output_vcf,
        seed=args.seed
    )
