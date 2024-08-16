import subprocess
import re
import functools
from operator import mul

#mersenne_factors[i] lists the factors of ith mersenne number (2^i - 1)
mersenne_factors = [[]]

def is_power_of_2(n):                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    return (n == 0) or (n & (n - 1) == 0)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    

def sord_2(n):  
    m = 1
    while True:
        two_powered  = 2 ** (m + 1)
        two_n_plus_1 = 2 * n + 1

        remainder = two_powered % two_n_plus_1

        if remainder == 1 or remainder == 2 * n:
            return m + 1

        m += 1

def max_dividing_power_2(n):
    if n == 1:
        return 0
    elif is_power_of_2(n):
        return n.bit_length() - 2
    else:
        return (n & -n).bit_length() - 1

def read_mersenne_prime_factors():
    global mersenne_factors

    mersenne_factor_counts = []

    with open("b046051.txt", "r") as f:
        for l in f:
            line_match = re.match("(\\d+)\\s+(\\d+)\\s*", l)
            if line_match is not None:
                mersenne_factor_counts.append(int(line_match[2]))

    with open("b001265.txt", "r") as f:
        current_number = 2 #Start from 2^2 - 1
        while True:
            if current_number > len(mersenne_factor_counts):
                break

            mersenne_factors.append([])
            for _ in range(0, mersenne_factor_counts[current_number - 1]):
                line = f.readline()
                if len(line) == 0:
                    break

                line_match = re.match("(\\d+)\\s+(\\d+)\\s*", line)
                if line_match is not None:
                    mersenne_factors[-1].append(int(line_match[2]))

            current_number += 1

def find_cycle_length_cyclic_rule_150(n, verifier_path):
    topology = "torus"
    rule = "--rule-150"

    # From Stephen Wolfram in https://content.wolfram.com/sw-publications/2020/07/algebraic-properties-cellular-automata.pdf:
    # For even n, the cycle length П(n) is 2*П(n/2)
    # For odd n, the cycle length П(n) divides 2^sord(n) - 1
    # In other words, П(n) divides 2^k * (2^sord(n/2^k) - 1), where k is the highest power of 2 dividing n
    two_exponent = max_dividing_power_2(n)
    sord_number  = n // (2 ** two_exponent)

    sord_exponent = sord_2((sord_number - 1) // 2) #Sord_2(n) is computed for [(n - 1)/2]th odd number

    power_two_factor = 2 ** two_exponent
    mersenne_factor  = 2 ** sord_exponent - 1

    cycle_length_multiple = power_two_factor * mersenne_factor

    cycle_proc = subprocess.run(f"{verifier_path} --topology {topology} {rule} --verify_cycle_length {cycle_length_multiple} --size {n} --offset {cycle_length_multiple * 2}", capture_output=True)
    cycle_str = str(cycle_proc.stdout, encoding="utf-8")

    if not "success" in cycle_str:
        print(f"Cycle length verification error for n = {n}: {cycle_length_multiple}", cycle_proc)
        return 0, 0

    cycle_offset_proc = subprocess.run(f"{verifier_path} --topology {topology} {rule} --calc_cycle_offset {cycle_length_multiple} --size {n}", capture_output=True)
    cycle_offset_str  = str(cycle_offset_proc.stdout, encoding="utf-8")

    cycle_offset_match = re.match("Cycle length offset is (\\d+).\\s*", cycle_offset_str)
    if cycle_offset_match is None:
        print(f"Error validating cycle offset for n = {n}: {cycle_offset_str}")
        return 0, 0

    cycle_offset = int(cycle_offset_match[1])

    cycle_proc_with_offset = subprocess.run(f"{verifier_path} --topology {topology} {rule} --verify_cycle_length {cycle_length_multiple} --size {n} --offset {cycle_offset}", capture_output=True)
    cycle_str_with_offset  = str(cycle_proc_with_offset.stdout, encoding="utf-8")

    if "success" not in cycle_str_with_offset:
        print(f"cycle offset calculated incorrectly for n = {n}: {cycle_offset}")
        return 0, 0

    cycle_length_divisors = [2] * two_exponent + mersenne_factors[sord_exponent - 1]
    assert functools.reduce(mul, cycle_length_divisors, 1) == cycle_length_multiple

    #Try cycle length without each of the divisors; remove unnecessary divisors
    is_cycle_minimum = False
    reduced_cycle_divisors = cycle_length_divisors
    while not is_cycle_minimum:
        reduced_cycle_length = functools.reduce(mul, reduced_cycle_divisors, 1)

        for index, p in enumerate(reduced_cycle_divisors):
            reduced_cycle_multiple = reduced_cycle_length // p

            reduced_cycle_proc = subprocess.run(f"{verifier_path} --topology {topology} {rule} --verify_cycle_length {reduced_cycle_multiple} --size {n} --offset {cycle_offset}", capture_output=True)
            reduced_cycle_str = str(reduced_cycle_proc.stdout, encoding="utf-8")

            if "success" in reduced_cycle_str:
                reduced_cycle_divisors = reduced_cycle_divisors[:index] + reduced_cycle_divisors[index + 1:]
                break

        else:
            is_cycle_minimum = True

    reduced_cycle_length = functools.reduce(mul, reduced_cycle_divisors, 1)
    return reduced_cycle_length, cycle_length_multiple

if __name__ == "__main__":
    read_mersenne_prime_factors()

    verifier_path = "../out/build/x64-Release/CycleLengthVerifier.exe"

    with open("b085588.txt", "w") as f:
        for n in range(1, 1000 + 1):    
            cycle_length, cycle_length_heuristic = find_cycle_length_cyclic_rule_150(n, verifier_path)

            if cycle_length == 0:
                break

            if cycle_length != cycle_length_heuristic:
                print(f"The maximum cycle length for rule 150 of width {n} is {cycle_length} (differs from the initial heuristic by a factor of {cycle_length_heuristic // cycle_length})")
            else:
                print(f"The maximum cycle length for rule 150 of width {n} is {cycle_length}.")

            f.write(f"{n} {cycle_length}\n")

        f.close()