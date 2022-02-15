######################################################################
# set_cover_greedy.py
# Implementation of a greedy approximation to the minimum set cover
# problem.
#
# The greedy algorithm is analyzed in Kleinberg and Tardos, and is
# shown to approximate the optimal solution up to a factor of log(n).
######################################################################

from data_structures.bucket_queue.bucket_queue import BucketQueue
import sys
import random

MIN_ELEMENT = 0
MAX_ELEMENT = int(10e6)
NOISE_PERCENTAGE = 0.10


def main():
    print(single_randomized_test(10, 25, 5, 10, True, True, set_cover_greedy))


def set_cover_greedy(sets, elements):
    """
    Use a greedy algorithm to find a set cover of a given collection of
    elements. Note that we assume that the set of elements is a subset
    of the union of the input sets.

    Inputs:
        sets (list[set()]) - an iterable of set objects.
        elements (set()) - a set of the elements to be covered.
    Outputs:
        set_cover (list[set()]) - a subset of the input sets that covers
                                  all of the elements. This is meant to be
                                  a greedy approximation to the minimum
                                  set cover.
    """
    # Convert set into a hash table.
    elements = {ele: False for ele in elements}
    n_false = len(elements)

    # The set cover to be returned.
    cover = list()

    # Add the sets to a bucket queue
    bq = BucketQueue(C=len(elements), L=0)
    sets = [(s, None) for s in sets]
    for i, (s, _) in enumerate(sets):
        # calculate size of intersection between this set and elements
        priority = 0
        for ele in s:
            if ele in elements:
                priority += 1
        # don't include sets with no relevant elements.
        if priority > 0:
            bq.insert(i, priority)
            sets[i] = (s, priority)

    # Unused sets - only include sets that intersect with elements.
    uncovered = set([i for i, x in enumerate(sets) if x[1]])

    while n_false > 0:
        # make sure at least one element is covered each round, otherwise break
        n_false_start = n_false

        # remove the set with the currently greatest number of uncovered
        # elements.
        set_added_index = bq.extract_max()

        # If set_added_index is None, or the maximum priority of the sets is zero
        # there are no more sets left, but some elements are still uncovered.
        if set_added_index is None or sets[set_added_index] == 0:
            print("Sets do not fully cover elements!", file=sys.stderr)
            break
        else:
            set_added, _ = sets[set_added_index]
            for ele in set_added:
                # want to skip elements that are irrelevant or already added
                if not elements.get(ele, True):
                    elements[ele] = True
                    n_false -= 1
            cover.append(set_added)
            uncovered.remove(set_added_index)

            # update the priorities of all of the other sets
            for remaining_set_index in uncovered:
                remaining_set, remaining_set_priority = sets[remaining_set_index]
                new_priority = 0
                for ele in remaining_set:
                    # want to skip elements that are irrelevant or already added
                    if not elements.get(ele, True):
                        new_priority += 1
                # minimum possible priority is zero
                bq.change_priority(remaining_set_index, remaining_set_priority,
                                   new_priority)
                sets[remaining_set_index] = (remaining_set, new_priority)

            # make sure at least one element is covered each round, otherwise break
            if n_false == n_false_start:
                print("Sets do not fully cover elements!", file=sys.stderr)
                break

    return cover


def run_randomized_tests(ntests, nmin_elements=50, nmax_elements=int(10e6),
                         nmin_sets=1, nmax_sets=1000, noise=False, cover_all=True,
                         algorithm=set_cover_greedy, **kwargs):
    """
    Creates a random set of elements in the range [nmin_elements, nmax_elements] 
    to be covered. Randomly generates sets so that the union of these sets covers 
    all or a subset of these elements (depending on boolean value of cover_all).

    If noise is True, some elements not needing to be covered will be added.

    Returns a boolean value indicating whether or not all of the elements are 
    covered by the returned greedy algorithm cover.  
    """
    assert nmin_sets >= 1, "Error: must have at least one set!"


def single_randomized_test(nmin_elements, nmax_elements, nmin_sets, nmax_sets,
                           noise, cover_all, algorithm, **kwargs):
    """
    Creates a random set of elements in the range [nmin_elements, nmax_elements] 
    to be covered. Randomly generates sets so that the union of these sets covers 
    all or a subset of these elements (depending on boolean value of cover_all).

    If noise is True, some elements not needing to be covered will be added to 
    the sets.

    Returns a boolean value indicating whether or not all of the elements are 
    covered by the returned greedy algorithm cover.  
    """
    nsets = random.randint(nmin_sets, nmax_sets)
    nelements = random.randint(nmin_elements, nmax_elements)

    set_of_all_elements = set()
    for _ in range(nelements):
        set_of_all_elements.add(random.randint(nmin_elements, nmax_elements))

    sets = generate_sets(set_of_all_elements, nsets, noise=True)

    cover = algorithm(sets, set_of_all_elements)

    set_of_all_elements = {x: True for x in set_of_all_elements}
    ct = len(set_of_all_elements)
    for included_set in cover:
        for x in included_set:
            if set_of_all_elements.get(x, False):
                ct -= 1
                set_of_all_elements[x] = False

    print(set(set_of_all_elements.keys()))
    print(cover)
    print(sets)

    return ct == 0


def generate_sets(universe, nsets, noise):
    """
    Randomly generates nsets sets from a universe of elements, such that the 
    union of the sets is equal to the universe. 
    """
    # the minimum and maximum elements in the set
    min_element = min(universe)
    max_element = max(universe)

    min_size = 1
    max_size = len(universe)

    universe = list(universe)
    added = [False for _ in range(len(universe))]
    sets = []
    for _ in range(nsets-1):
        size = random.randint(min_size, max_size)
        new_set = set()
        for _ in range(size):
            new_element_index = random.randint(0, len(universe)-1)
            new_element = universe[new_element_index]
            added[new_element_index] = True
            new_set.add(new_element)
            # add noise element with specified probability
            if random.random() < NOISE_PERCENTAGE:
                # if random.random() < 0.50:
                #     new_set.add(random.randint(
                #         max_element + 1, MAX_ELEMENT + 1000))
                # else:
                #     new_set.add(random.randint(
                #         MIN_ELEMENT - 1000, min_element - 1))
                new_set.add(random.randint(max_element+1, max_element+1000))
        sets.append(new_set)

    last_set = {x for i, x in enumerate(universe) if not added[i]}
    additional_size = random.randint(min_size - len(last_set),
                                     max_size - len(last_set))
    for i in range(0, additional_size):
        new_element_index = random.randint(0, len(universe)-1)
        new_element = universe[new_element_index]
        last_set.add(new_element)
    sets.append(last_set)

    return sets


if __name__ == '__main__':
    main()
