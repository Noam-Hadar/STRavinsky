def reSTRucture(sequence):
    sequences = []
    for i in range(len(sequence)):
        shuffled = sequence[i:] + sequence[:i]
        sequences.append(shuffled)
        sequences.append(shuffled[::-1])
        complementary = shuffled.replace('A','t').replace('T','a').replace('G','c').replace('C','g').upper()
        sequences.append(complementary)
        sequences.append(complementary[::-1])
    reSTRuctured = sorted(sequences)[-3]
    return reSTRuctured
