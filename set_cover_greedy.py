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


def main():
    elements = set([1, 3, 5])
    sets = [set([1, 3]), set([5])]
    set_cover_greedy(sets, elements)

    elements = set([1, 3, 5, 6, 7, 8, 9, 10, 21, 48, 59, 62, 71, 80, 99, 100])
    set1 = set([1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 48])
    set2 = set([59, 62, 71, 80, 81, 99])
    # set3 = set([101, 102])
    sets = [set1, set2]
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
        bq.insert(i, priority)
        sets[i] = (s, priority)

    # Unused sets
    uncovered = set([i for i in range(len(sets))])

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


def randomized_test(ntests, nmin_elements=50, nmax_elements=int(10e6),
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
    raise NotImplementedError


if __name__ == '__main__':
    main()
