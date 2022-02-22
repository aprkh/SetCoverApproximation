######################################################################
# set_cover_greedy.py
# Implementation of a greedy approximation to the minimum set cover
# problem.
#
# The greedy algorithm is analyzed in Kleinberg and Tardos, and is
# shown to approximate the optimal solution up to a factor of log(n).
######################################################################

from data_structures.bucket_queue.bucket_queue import BucketQueue
from data_structures.eQTL_genotype_interface.People import People
import numpy as np
import sys
import random

MIN_ELEMENT = 0
MAX_ELEMENT = int(10e6)
NOISE_PERCENTAGE = 0.10
REMOVE_PERCENTAGE = 0.20


def main():
    # print(run_randomized_tests(25, 10, 25, 5, 10, True, False, set_cover_greedy))

    maternal_file = "data_structures\\eQTL_genotype_interface\\Xm1.txt"
    paternal_file = "data_structures\\eQTL_genotype_interface\\Xp1.txt"

    loci_to_cover = set([1, 4, 18, 12, 10, 19, 18, 15])
    people = People(np.loadtxt(maternal_file)[
                    1:, :], np.loadtxt(paternal_file)[1:, :])

    het_sets = people.get_heterozygote_loci()

    print(het_sets)
    print(set_cover_greedy(het_sets, loci_to_cover))


def set_cover_greedy(sets, elements):
    """
    Use a greedy algorithm to find a set cover of a given collection of
    elements. Note that we assume that the set of elements is a subset
    of the union of the input sets.

    Inputs:
        sets (list[set()]) - an iterable of set objects.
        elements (set()) - a set of the elements to be covered.
    Outputs:
        cover (list[set()]) - a subset of the input sets that covers
                                  all of the elements. This is meant to be
                                  a greedy approximation to the minimum
                                  set cover.
        cover_index (list[int]) - the indices corresponding to which sets 
                                  are included in the cover. 
    """
    # Convert set into a hash table.
    elements = {ele: False for ele in elements}
    n_false = len(elements)

    # The set cover to be returned.
    cover = list()
    cover_index = list()

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
            cover_index.append(set_added_index)
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

    return cover, cover_index


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
    results = [None for _ in range(ntests)]
    for i in range(ntests):
        results[i] = single_randomized_test(nmin_elements, nmax_elements,
                                            nmin_sets, nmax_sets, noise,
                                            cover_all, algorithm, **kwargs)
    return all(results)


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

    sets = generate_sets(set_of_all_elements, nsets,
                         noise, cover_all)

    cover, cover_index = algorithm(sets, set_of_all_elements)

    # If we want everything covered, then there should be no uncovered elements
    # left in set_of_all_elements after we count all the elements in the returned
    # cover. Alternatively, if we have some missing elements from the set, we should
    # still cover every element in the intersection of set_of_all_elements and sets
    # in our returned cover.
    result = None
    if cover_all:
        set_of_all_elements = {x: True for x in set_of_all_elements}
        ct = len(set_of_all_elements)
        for included_set in cover:
            for x in included_set:
                if set_of_all_elements.get(x, False):
                    ct -= 1
                    set_of_all_elements[x] = False
        result = (ct == 0)
    else:
        possible_elements = set()
        for included_set in sets:
            for element in included_set:
                if element in set_of_all_elements:
                    possible_elements.add(element)
        covered_elements = set()
        for included_set in cover:
            for element in included_set:
                if element in set_of_all_elements:
                    covered_elements.add(element)
        result = (covered_elements == possible_elements)

    print(set_of_all_elements)
    print(cover, cover_index)
    print(sets)

    return result


def generate_sets(universe, nsets, noise, cover_all):
    """
    Randomly generates nsets sets from a universe of elements, such that the 
    union of the sets is equal to the universe. 
    """
    # the minimum and maximum elements in the set
    min_element = min(universe)
    max_element = max(universe)

    min_size = 1
    max_size = len(universe)

    # may or may not cover all elements
    universe_list = []
    for x in universe:
        if cover_all:
            universe_list.append(x)
        else:
            # randomly remove elements with pre-specified probability
            if random.random() > REMOVE_PERCENTAGE:
                universe_list.append(x)
    universe = universe_list
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
    for _ in range(0, additional_size):
        new_element_index = random.randint(0, len(universe)-1)
        new_element = universe[new_element_index]
        last_set.add(new_element)
    sets.append(last_set)

    return sets


if __name__ == '__main__':
    main()
