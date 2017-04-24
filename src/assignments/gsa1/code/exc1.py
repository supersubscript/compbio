#!/usr/bin/env python

import numpy as np
import random
import math
import matplotlib.pyplot as plt
from operator import itemgetter
from string import maketrans


#################################
# PROGRAM 1
#################################
# Get index for letter (symbol) occurring.
def get_ind(probs):
    rand = random.random()
    for ii in xrange(len(probs)):
        if rand < probs[ii]:
            return ii


# Read in transition matrix and initial distribution from files.
# If sumrows == TRUE, sum the rows
def read_and_output(mat_file, idist_file, seq_length, sumrows, transpose):
    # Read in data
    trans_mat = np.loadtxt(open(mat_file, "r"))
    init_dist = np.loadtxt(open(idist_file, "r"))

    # Act on the optional things
    if transpose:
        trans_mat = np.transpose(trans_mat)
    if sumrows:
        trans_mat = np.cumsum(trans_mat, 1)

    # Fill out the output probabilistically
    output = seq_length * [None]
    output[0] = np.random.choice(np.arange(0, len(trans_mat)), p=list(init_dist))
    for ii in xrange(1, seq_length):
        output[ii] = get_ind(trans_mat[output[ii - 1], :])

    "".join(str(s) for s in output)
    return output


mat_file = "/home/henrik/compbio/src/assignments/gsa1/report/mat_file.dat"
idist_file = "/home/henrik/compbio/src/assignments/gsa1/report/idist_file.dat"
output = read_and_output(mat_file, idist_file, 115, True, False)


# print output


#################################
# PROGRAM 2
#################################
# Infer the transition matrix and initial distribution from a sequence.
def infer_mat_and_init(sequence):
    # Get some data and create the matrices we need.
    elements = np.unique(list(sequence))
    trans_mat = np.zeros((len(elements), len(elements)))
    init_dist = np.zeros(len(elements))  # This is specific to this assignment
    init_dist[int(sequence[0])] = 1  # Count the first one

    # Count the occurrences of all events
    for ii in xrange(1, len(sequence)):
        trans_mat[int(sequence[ii - 1]), int(sequence[ii])] += 1

    # What are the probabilities?
    for ss in xrange(len(trans_mat)):
        trans_mat[ss, :] /= np.sum(trans_mat[ss,])

    return trans_mat, init_dist


trans_mat, init_dist = infer_mat_and_init("010101010101111010010101100101001010101010101010010101")

#################################
# PROGRAM 3
#################################
S = np.array([0, 1])
V = np.array([0, 1, 2, 3, 4])
A = np.matrix("0.8, 0.2;"
              "0.1, 0.9")
B = np.matrix("0.2, 0.5, 0.2, 0.1, 0.0;"
              "0.0, 0.1, 0.4, 0.4, 0.1")
mu0 = np.array([.5, .5])


# Emit a sequence given some transition and emission matrices, as well as
# an initial distribution.
def emit_sequence(seq_length, trans_matrix, emiss_matrix, init_dist):
    c_trans_mat = np.cumsum(trans_matrix, 1)
    c_emiss_mat = np.cumsum(emiss_matrix, 1)

    emit = np.zeros(seq_length)
    hidden = np.zeros(seq_length)
    hidden[0] = np.random.choice(np.arange(0, len(trans_matrix)), p=list(init_dist))
    emit[0] = get_ind((c_emiss_mat[int(hidden[0]), :]).tolist()[0])

    # Parse sequence
    for ii in xrange(1, seq_length):
        hidden[ii] = get_ind((c_trans_mat[int(hidden[ii - 1]), :]).tolist()[0])
        emit[ii] = get_ind((c_emiss_mat[int(hidden[ii]), :]).tolist()[0])

    emit = "".join(str(int(x)) for x in emit)
    hidden = "".join(str(int(x)) for x in hidden)
    return emit, hidden


emit, hidden = emit_sequence(115, A, B, mu0)


# random.seed(2)
# plt.interactive(True)
# plt.yticks([])
# pl_emit, = plt.step(xrange(len(emit)), list(emit), label = "Emit")
# pl_hidd, = plt.step(xrange(len(hidden)), list(int(s) +4 for s in hidden), label = "Hidden")
# plt.legend(handles=[pl_emit, pl_hidd])
# plt.show()
# print(emit)a


#################################
# PROGRAM 4
#################################
# Read in sequence from file.
def read_sequence(file, skipFirst):
    seq = ""
    with open(file, 'r') as in_file:
        if skipFirst:
            in_file.readline()
        seq = in_file.read().replace('\n', '')
    return seq


# Use the forward algorithm to determine (log) likelihood of sequence
def forward(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist):
    # Function for getting normalisation constant
    def get_constant(a):
        sum = np.sum(a)
        return 1.0 / len(trans_matrix) if sum == 0 else 1.0 / sum

    # Pre-allocate
    fw = np.zeros((len(hidden_states), len(sequence)))
    constants = []

    # Initialise
    for ss in xrange(len(hidden_states)):
        fw[ss, 0] = init_dist[ss] * emiss_matrix[ss, int(sequence[0])]
    const = get_constant(fw[0, :])
    constants.append(const)

    # Scale accordingly
    for ss in xrange(len(hidden_states)):
        fw[ss, 0] *= const

    # Now go through the rest of the sequence
    probs = np.zeros(len(hidden_states))
    for ii in xrange(1, len(sequence)):
        for ss in xrange(len(hidden_states)):
            probs[ss] = np.sum(
                (fw[substate, ii - 1] * trans_matrix[substate, ss] * emiss_matrix[ss, int(sequence[ii])]) for
                substate in xrange(len(hidden_states)))

        # Scale again
        const = get_constant(probs)
        constants.append(const)
        for ss in xrange(len(hidden_states)):
            fw[ss, ii] = const * probs[ss]

    # Sum up the probabilities; disregard if 0
    ln_prob = -np.sum([math.log(const) if const != 0 else 0 for const in constants])
    return ln_prob, constants, fw


out_sequence, __ = emit_sequence(115, A, B, mu0)
seq_file = '/home/henrik/compbio/src/assignments/gsa1/random_output_sequence.dat'
with open(seq_file, "w") as f_out:
    f_out.write(out_sequence)
    f_out.close()
that_same_sequence = read_sequence(seq_file, False)
prob, consts, fw = forward(that_same_sequence, A, B, S, mu0)


# print(math.exp(prob))


#################################
# PROGRAM 5
#################################
# We cheat a little and put the remnants in a bin of their own, i.e. if the last bin only has e.g. 20 elements, the GC
# content in this will be "the number of GC bases in those 20 bases"/bin_size, which is wrong, but we'll have to
# live with it.
def calc_gc(sequence, bin_size):
    gc_cont = []
    bin_size = float(bin_size)
    for ii in xrange(int(math.ceil(float(len(sequence)) / bin_size))):
        ceil = int((ii + 1) * bin_size)
        if ceil > len(sequence):
            ceil = len(sequence)
        bin = sequence[int(ii * bin_size): ceil]
        G = bin.count("G")
        C = bin.count("C")
        gc_cont.append((G + C) / bin_size)
    return gc_cont


# Read genome
sc_file = '/home/henrik/compbio/src/assignments/gsa1/sc_gen.fa'
sc_gen = read_sequence(sc_file, True)
gc = calc_gc(sc_gen, 100)


# Relabel according to predefined cuts
def relabel(seq, cuts):
    relab = []
    for ii in seq:
        for cut in xrange(len(cuts)):
            if ii <= cuts[cut]:
                relab.append(cut)
                break

    relab = "".join(str(x) for x in relab)
    return relab


# Compare the two strings
relab = relabel(gc, [.285, .342, .40, .5, 1.])  # Split by percentages
model3_seq, model3_hidden = emit_sequence(len(relab), A, B, mu0)
log_p, const, fw_seq = forward(relab, A, B, S, mu0)


# plt.hist(list(int(s) for s in model3_seq), bins=5, alpha=0.5, label='x')
# plt.hist(list(int(s) for s in relab), bins = 5, alpha=0.5, label='y')
# # plt.show()


#################################
# PROGRAM 6
#################################
# We only use this for Baum-Welch, so we use the same constants as in our forward case
def backward(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist, constants):
    # Preallocate for our backwards sequence
    bw = np.zeros((len(hidden_states), len(sequence)))

    # Initialise
    for ss in xrange(len(hidden_states)):
        bw[ss, len(sequence) - 1] = 1.
        bw[ss, len(sequence) - 1] = constants[len(sequence) - 1]

    # Go through the rest of the sequence, continue from the back end
    for ii in xrange(len(sequence) - 2, -1, -1):
        b = np.zeros(len(hidden_states))
        for ss in xrange(len(hidden_states)):
            b[ss] = np.sum([(bw[substate, ii + 1] * trans_matrix[ss, substate] *
                             emiss_matrix[substate, int(sequence[ii + 1])]) for substate in xrange(len(hidden_states))])

        # Again, scale
        for ss in xrange(len(hidden_states)):
            bw[ss, ii] = constants[ii] * b[ss]

    ln_prob = -np.sum([math.log(c) if c != 0 else 0 for c in constants])
    return ln_prob, bw


def baum_welch(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist):
    # Get those lists we need
    p1, constants, fw = forward(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist)
    p1, bw = backward(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist, constants)

    # Initial distribution
    for ss in xrange(len(hidden_states)):
        init_dist[ss] = fw[ss, 0] * constants[0] * bw[ss, 0]
    init_dist /= np.sum(init_dist)

    # New state matrix
    for ss in xrange(len(hidden_states)):
        for substate in xrange(len(hidden_states)):
            # Sum up the estimate, normalise and assign
            new_est = np.sum([(fw[ss, ii] * bw[substate, ii + 1] * trans_matrix[ss, substate] * emiss_matrix[
                substate, int(sequence[ii + 1])]) for ii in xrange(len(sequence) - 1)])
            norm = np.sum([(fw[ss, ii] * bw[ss, ii] / constants[ii]) for ii in xrange(len(sequence) - 1)])
            trans_matrix[ss, substate] = new_est / norm if norm != 0 else 1. / len(hidden_states)

    # New emission matrix (structure same as above)
    for ss in xrange(len(hidden_states)):
        for jj in xrange(emiss_matrix.shape[1]):
            new_est = np.sum([(fw[ss, ii] * bw[ss, ii] / constants[ii] if int(sequence[ii]) == jj else 0) for ii in
                              xrange(len(sequence))])
            norm = np.sum([(fw[ss, ii] * bw[ss, ii] / constants[ii]) for ii in xrange(len(sequence))])
            emiss_matrix[ss, jj] = new_est / norm if norm != 0 else 1. / emiss_matrix.shape[1]
    return p1


#################################
# PROGRAM 7
#################################
def viterbi(sequence, trans_matrix, emiss_matrix, hidden_states, init_dist):
    # Define some things we're gonna use.
    constants = np.zeros(len(sequence))
    solution = np.zeros(len(sequence))
    dyn_mat = np.zeros((len(hidden_states), len(sequence)))
    store_path = np.zeros((len(hidden_states), len(sequence)))

    # Initialise rows
    dyn_mat[:, 0] = init_dist * emiss_matrix[:, int(sequence[0])]
    constants[0] = 1.0 / np.sum(dyn_mat[:, 0])
    dyn_mat[:, 0] *= constants[0]

    # Fill out the table
    for ii in xrange(1, len(sequence)):
        for ss in xrange(len(hidden_states)):
            probabilities = np.diag(dyn_mat[:, ii - 1]) * trans_matrix[:, ss]  # Probability of coming from either state
            store_path[ss, ii], dyn_mat[ss, ii] = max(enumerate(probabilities), key=itemgetter(1))
            dyn_mat[ss, ii] *= emiss_matrix[ss, int(sequence[ii])]

        # Now we scale
        constants[ii] = 1.0 / np.sum([dyn_mat[ss, ii] for ss in xrange(len(hidden_states))])
        dyn_mat[:, ii] *= constants[ii]

    # Backtrack
    solution[len(sequence) - 1] = dyn_mat[:, len(sequence) - 1].argmax()  # last state
    for ii in xrange(len(sequence) - 1, 0, -1):
        solution[ii - 1] = store_path[int(solution[ii]), ii]

    return solution


model7_seq = sc_gen
gc = calc_gc(sc_gen, 100)

intab = "ATCG"
outtab = "0123"
trantab = maketrans(intab, outtab)
model7_seq = model7_seq.translate(trantab)


prev = -99999999999999
change = -prev
while change > 0.0001:
    prob = baum_welch(relab, A, B, S, mu0)
    change = (prob - prev)
    prev = prob

solution = viterbi(relab, A, B, S, mu0)
# out_len = 200
# solution = solution[0:out_len]
# plt.interactive(True)
# plt.yticks([])
# plt.step(xrange(len(solution)), solution + 5, label="Viterbi")
# plt.plot(xrange(out_len), [(x + 4) for x in gc[0:out_len]], label="GC content")
# plt.step(xrange(out_len), list(relab[0:out_len]), label="Categorized GC content")
# plt.legend()
# plt.show()

plt.interactive(True)
# plt.yticks([])
plt.hist(gc, bins=100, width=1./100. - 0.006, label="GC")
plt.legend()
fig = plt.gcf()
fig.set_size_inches(6, 4)
fig.savefig("/home/henrik/compbio/src/assignments/gsa1/report/gc_hist.png")
plt.show()
