#!/usr/bin/python
# This code is from python.org's sandbox for collections.
# http://svn.python.org/projects/sandbox/trunk/collections/pairing_heap.py
# It is a copy, completely not written by me.
"""

Note: this implementation is not currently thread-safe.

The code in this file implements a pairing heap.  By default it
operates as a min-heap, but it can also work as a max-heap via a flag
passed to the heap's __init__ routing.

The important papers on pairing heaps are these:

[1] Michael L. Fredman, Robert Sedjewick, Daniel D. Sleator, and
Robert E. Tarjan.  The Pairing Heap: A New Form of Self-Adjusting
Heap.  Algorithmica 1:111-129, 1986.

[2] John T. Stasko and Jeffrey Scott Vitter.  Pairing heaps:
experiments and analysis.  Communications of the ACM, Volume 30, Issue
3, Pages 234-249, 1987.

[3] Michael L. Fredman.  On the efficiency of pairing heaps and
related data structures.  Journal of the ACM, Volume 46, Issue 4,
Pages 473-501, 1999.

There are actually several varieties of pairing heaps.  Five are
described in [1], and two more are described in [2].  This
implementation is the basic pairing heaps described in [1], also known
as the twopass variant.  There isn't much published evidence that any
of the variants are significantly better than the others.

The twopass pairing heap works as follows.

A pairing heap is a tree that maintains the heap property: each node
has a lower value than all its children.  The pairing heap tree need
not be binary nor balanced.  The good amortized efficiency of the
pairing heap results from clever algorithms for .adjust_key() and
.delete().

The basic operation in the pairing heap is the comparison-link.  Two
nodes are compared, and the largest is linked as a child of the
smaller.  To represent the tree, each node maintains a pointer to the
leftmost child and to its right sibling.  Each node also maintains a
pointer to its parent.  When a comparison-link is performed, the new
child becomes the leftmost child.

[1] suggests a memory optimization where only the right-most sibling
stores a pointer to the parent, and each node maintains a single bit
indicating if it is that rightmost sibling.  This theoretically saves
the size of a pointer in each node, but only if we can store a bit for
free.  However, storing a bit takes just as much space as storing a
pointer, due to padding.  So, we bite the bullet and store the
pointer.  Also, [3] shows that the amortized cost of .adjust_key()
is at least Omega(log log n) if the parent pointer is *not* stored.
The amortized cost may be better using the parent pointer.

For .adjust_key(), the topmost element is removed and set aside as
the return value.  A left-to-right pass is made over its children,
comparison-linking them pairwise.  For example, if there are 6
children labeled A, B, C, D, E, and F, A and B will be
comparison-linked, as will C and D, and E and F.  Next, the rightmost
node is repeatedly comparison-linked with its sibling until there is
only one node left.  This node is the new root.

For .delete(), first the node is cut from its parent.  This results in
a tree rooted at the node being deleted.  The .adjust_key()
procedure is invoked on the node, then the resulting tree is
comparison-linked with the original tree.

See [1] for more details.  It also includes figures which are helpful
for understanding.

"""

import logging

logger=logging.getLogger("pairing_heap")

class WrongHeap (Exception):
    """The user attempted to .delete() or .adjust() from the wrong heap.

    This can also occur if the user attempts to delete or adjust a
    node that has already been removed from the heap.
    """
    def __str__(self):
        return 'WrongHeap'

class Underflow (Exception):
    "The user attempted to .peek() or extract() from an empty heap."
    def __str__(self):
        return 'Underflow'

class IncompatibleHeaps (Exception):
    "The user attempted to meld heaps that have different compare operations."
    def __str__(self):
        return 'IncompatibleHeaps'

class InternalError (Exception):
    "This should never get raised in practice."
    def __str__(self):
        return 'InternalError'

class WrongAdjustKeyDirection (Exception):
    "The user attempted to adjust a key in the wrong direction."
    def __str__(self):
        return 'WrongAdjustKeyDirection'

class heap_node:
    """This class is used internally to implement the heap.

    Users will have be able to see references to this data structure,
    as this is neccessary so they can pass the reference to the
    .adjust_key() and .delete() methods.  However, the data structure
    should be considered opaque.  The only public methods are .value()
    and .values().

    Users should not instantiate heap_node on their own.

    In the C implementation, all of this will be enforced.

    """
    __slots__ = ('_item', '_child', '_sibling', '_parent')

    def __init__(self, item):
        self._item = item
        self._clean()

    def value(self):
        "Return the value associated with this node."
        return self._item

    def _clean(self):
        self._child = None
        self._sibling = None
        self._parent = None

    def _add_child(self, node):
        node._parent = self
        node._sibling = self._child
        self._child = node

    def _cut(self):
        "Remove this node and its children from this node's parent."
        n = self._parent
        if id(n._child) == id(self):
            n._child = self._sibling
        else:
            n = n._child
            try:
                while id(n._sibling) != id(self):
                    n = n._sibling
                n._sibling = self._sibling
            except AttributeError:
                raise InternalError
        self._parent = None
        self._sibling = None

    def _extract(self, heap):
        "Link this node's children as a prelude to this node's removal"

        #  Make the left-to-right pairing pass
        children = []
        n = self._child
        while n and n._sibling:
            next = n._sibling._sibling
            pair = n._sibling
            n._parent = None
            n._sibling = None
            pair._parent = None
            pair._sibling = None
            children.append(heap._link(n, pair))
            n = next
        if n:
            children.append(n)
        if not children:
            return None

        # Repeated comparison-link to the rightmost node
        root = children[-1]
        for i in range(len(children)-2,-1,-1):
            root = heap._link(root, children[i])
        root._parent = None

        # Return the new root
        return root

    def values(self):
        "Return an unsorted sequence of the values of this part of the tree."
        rv = [self._item,]
        if self._child:
            rv += self._child.values()
        if self._sibling:
            rv += self._sibling.values()
        return rv

class pairing_heap:
    "Implements a min-heap using the pairing-heap algorithm."

    def __init__(self, values=None, key=None):
        logger.debug("pairing heap key {0}".format(key))
        if key is not None:
            self._key = lambda x: key(x._item)
        else:
            self._key = lambda x: x._item
        self._root = None
        self._len = 0
        if values is not None:
            for i in values:
                self.insert(i)

    def empty(self):
        "True if there is nothing in the heap."
        return self._root is None

    def peek(self):
        "Return, but do not remove, the top of the heap."
        if self._root:
            return self._root._item
        else:
            raise Underflow

    def __len__(self):
        "Return the number of items in the heap."
        return self._len

    def _link(self, node1, node2):
        "Perform the comparison-link operation."
        if not node1: return node2
        if not node2: return node1

        node1, node2=sorted([node1, node2], key=self._key)
        node1._add_child(node2)
        return node1

    def insert(self, item):
        "Insert an item into the heap."
        node = heap_node(item)
        self._root = self._link(node, self._root)
        self._len += 1
        return node

    def meld(self, other):
        """Merge another heap into this one.

        other will be an empty heap after completion
        """

        if not isinstance(other, pairing_heap):
            raise TypeError

        # Check that these heaps are compatible
        if self._key != other._key:
            raise IncompatibleHeaps

        self._root = self._link(other._root, self._root)
        self._len += other._len

        # Destroy other
        other._len = 0
        other._root = None

    def adjust_key(self, node, item):
        """Adjust the value of node to item.

        This must be decrease for a min-heap or an increase for a max-heap.

        If node is not part of this heap, both heaps may become gibberish.
        """
        self._check_heap_node (node)
        if self._cmp(node._item, item) < 0:
            raise WrongAdjustKeyDirection
        node._item = item
        if id(self._root) == id(node):
            return

        # Cut this node out of the tree
        node._cut()

        # Link it with the root
        self._root = self._link(node, self._root)

    def delete(self, node):
        """Delete a node from the middle of the heap.

        If node is not part of this heap, both heaps may become gibberish.
        """
        self._check_heap_node (node)
        if id(self._root) == id(node):
            self.extract()
            return
        node._cut()
        self._root = self._link (node._extract(self), self._root)
        self._len -= 1
        node._clean()
        return node._item

    def extract(self, n=0):
        "Extract the top value of the heap."
        if n > 0:
            return [self.extract() for i in xrange(n)]
        if not self._root:
            raise Underflow
        old_root = self._root
        self._root = self._root._extract(self)
        self._len -= 1
        old_root._clean()
        return old_root._item

    def extract_all(self):
        "Empty the heap into a sorted list of all the values"
        return [self.extract() for i in xrange(self._len)]

    def values(self):
        "Return an unsorted list of all the values in the heap"
        return self._root.values()

    def _check_heap_node(self, node):
        """Raise an error if node will cause problems.

        This function is a simple check to catch certain user errors.
        It will raise TypeError if node isn't a heap_node object at
        all.  It *may* raise WrongHeap if node isn't part of this heap.
        It will return if node is part of this heap.

        However, no guarantees are made about the return value if node
        is a heap_node, but part of a different heap.

        """

        if not isinstance(node, heap_node):
            raise TypeError
        if self._root == node:
            return
        if not node._parent:
            raise WrongHeap
        return
