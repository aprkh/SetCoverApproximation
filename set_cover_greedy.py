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


def main():
    elements = set([1, 3, 5])
    sets = [set([1, 3]), set([5])]
    set_cover_greedy(sets, elements)

    elements = set([1, 3, 5, 6, 7, 8, 9, 10, 21, 48, 59, 62, 71, 80, 99, 100])
    set1 = set([1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 48])
    set2 = set([59, 62, 71, 80, 81, 99])
    set3 = set([101, 102])
    sets = [set1, set2, set3]
    set_cover_greedy(sets, elements)


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
    bq = BucketQueue(C=len(elements), L=1)
    sets = [(s, None) for s in sets]
    for i, (s, _) in enumerate(sets):
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

        # If set_added_index is None, there are no more sets left, but some
        # elements are still uncovered.
        if set_added_index is None:
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
                new_priority = remaining_set_priority
                for ele in remaining_set:
                    # want to skip elements that are irrelevant or already added
                    if elements.get(ele, False):
                        remaining_set_priority -= 1
                bq.change_priority(remaining_set_index,
                                   remaining_set_priority, new_priority)
                sets[remaining_set_index] = (remaining_set, new_priority)

            # make sure at least one element is covered each round, otherwise break
            if n_false == n_false_start:
                print("Sets do not fully cover elements!", file=sys.stderr)
                break

        print("Number of remaining items:", n_false,
              "\nCurrent cover:", cover, file=sys.stderr)

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
        set_of_all_elements.add(random.randint(MIN_ELEMENT, MAX_ELEMENT))

    sets = generate_sets(set_of_all_elements, nsets, noise=False)

    cover = algorithm(sets, set_of_all_elements)

    set_of_all_elements = {x: False for x in set_of_all_elements}
    ct = len(set_of_all_elements)
    for included_set in cover:
        for x in included_set:
            if set_of_all_elements.get(x, False):
                ct -= 1
                set_of_all_elements[x] = True

    return ct == 0


def generate_sets(universe, nsets, noise):
    """
    Randomly generates nsets sets from a universe of elements, such that the 
    union of the sets is equal to the universe. 
    """
    raise NotImplementedError


if __name__ == '__main__':
    main()
