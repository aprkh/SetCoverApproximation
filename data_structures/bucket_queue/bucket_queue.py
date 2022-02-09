######################################################################
# bucket_queue.py
# Preliminary implementation of a bucket queue, used for implementation
# of the greedy algorithm used in the set cover problem.
######################################################################

import sys
import numpy as np


def main():
    bq = BucketQueue(50)
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    data = [(15, 25), (39, 40), (30, 15)]
    for x, p in data:
        bq.insert(x, p)
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)
    bq.change_priority(39, 40, 1)
    bq.change_priority(30, 15, 1)
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)
    bq.insert(2, 2)
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)
    print("removed a minimum element:", bq.extract_min())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a maximum element:", bq.extract_max())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a maximum element:", bq.extract_max())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a maximum element:", bq.extract_max())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a maximum element:", bq.extract_max())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("#" * 60)
    bq = BucketQueue(50)
    data = [(15, 25), (39, 40), (30, 15)]
    for x, p in data:
        bq.insert(x, p)

    bq.change_priority(39, 40, 50)
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a minimum element:", bq.extract_min())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a minimum element:", bq.extract_min())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a minimum element:", bq.extract_min())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)

    print("removed a minimum element:", bq.extract_min())
    print("lower bound:", bq.lowest)
    print("upper bound:", bq.greatest)
    print("array:", bq.A)


class BucketQueue():

    def __init__(self, C, L=1, data_structure=set, **kwargs):
        """
        C is the highest possible priority, while L is the lowest possible priority.
        data_structure is the type of data structure to use as the "buckets"
            - the data_structure object must have the following operations:
                - add(x) - add an element x
                - remove(x) - remove an element x from the collection
                - len(data_structure) - returns the number of elements in the data structure. 
        """
        self.C = C
        self.L = L
        self.A = [data_structure(**kwargs) for _ in range(self.L, self.C+1)]
        # lower and upper bounds of priorities with non-empty containers
        self.lowest = None
        self.greatest = None

    def insert(self, x, p):
        """
        Insert an element x with priority p.
        """
        assert self.L <= p <= self.C, "invalid value for p"
        p = self.priority_to_index(p)
        self.A[p].add(x)
        self.update_bounds_eager(p)

    def change_priority(self, x, p, q):
        """
        Change the priority of an element x from priority p to priority q.
        Note: x's current priority p must be specified.
        """
        assert self.L <= p <= self.C, "invalid value for p"
        assert self.L <= q <= self.C, "invalid value for q"
        p = self.priority_to_index(p)
        q = self.priority_to_index(q)
        self.A[p].remove(x)
        self.A[q].add(x)
        self.update_bounds_eager(q, p)  # we removed from p, moved to q

    def extract_max(self):
        """
        Extract the maximum method.
        """
        self.find_highest()
        if self.greatest is None:
            return None
        else:
            return self.A[self.greatest].pop()

    def extract_min(self):
        self.find_lowest()
        if self.lowest is None:
            return None
        else:
            return self.A[self.lowest].pop()

    #### Helper methods ####
    def priority_to_index(self, index):
        """
        Maps a given priority to an array index. 
        """
        return index - self.L

    def find_lowest(self):
        if self.lowest is None:
            assert self.greatest is None
            return
        for L in range(self.lowest, self.greatest+1):
            if len(self.A[L]) != 0:
                self.lowest = L
                return
        # No non-empty container found
        self.lowest = None
        self.greatest = None

    def find_highest(self):
        if self.greatest is None:
            assert self.lowest is None
            return
        for U in range(self.greatest, self.lowest-1, -1):
            if len(self.A[U]) != 0:
                self.greatest = U
                return
        # No non-empty container found
        self.lowest = None
        self.greatest = None

    def update_bounds_eager(self, p, q=None):
        """
        Update the lower and upper bounds of the non-empty priority buckets 
        given we've inserted an element with priority p. 
        Optionally, we may have removed an element with priority q, and will 
        want to take this into account. 
        """
        # update lower and upper bounds of priorities with non-empty containers
        if self.lowest is None:
            assert self.greatest is None,\
                "Error: cannot have no lowest priority if bucket queue is non-empty."
            self.lowest = p
            self.greatest = p
        else:
            assert self.greatest is not None,\
                "Error: cannot have no greatest priority if bucket queue is non-empty."
            self.lowest = min(p, self.lowest)
            self.greatest = max(p, self.greatest)
        if q is None or len(self.A[q]) != 0:
            return
        if q == self.greatest:
            while q >= self.lowest and len(self.A[q]) == 0:
                q -= 1
            if q < self.lowest:
                # then there are no more items in the bucket queue
                self.lowest = None
                self.greatest = None
            else:
                self.greatest = q
        elif q == self.lowest:
            while q < len(self.A) and len(self.A[q]) == 0:
                q += 1
            if q > self.greatest:
                # then there are no more elements in the bucket queue
                self.lowest = None
                self.greatest = None
            else:
                self.lowest = q

    def __str__(self):
        return str(self.A)


if __name__ == '__main__':
    main()
