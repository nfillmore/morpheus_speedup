using System;

namespace Morpheus
{
    // This class is similar to C#'s System.Collections.Generic.List class,
    // except that it is designed to minimize allocations and garbage
    // collections when the contained elements are heap-allocated objects.
    //
    // This is useful for, and designed for, reusable thread-local storage.
    //
    // Specifically, this class manages the allocation of its contained objects
    // internally, reusing them when possible. Thus, (i) the Add method returns
    // a pointer to an object instead of being passed a pointer to an object.
    // (ii) The Clear method does not release any object references. (iii) The
    // [] accessor is read-only. See comments on each method for more details.
    // 
    // Caveat: References to objects contained in this list are invalidated by
    // the Clear method. For example:
    // 
    //   FastListOfBoxes<MyObj> l;
    //   MyObj m = l.Add();
    //   DoSomethingWith(m);
    //   l.Clear();
    //   MyObj n = l.Add();
    //
    // After the last line above, m and n point to the same object, which is
    // probably not what was expected.
    public class FastListOfBoxes<T> where T: new()
    {
        private T[] buf;

        // The number of elements logically in the list.
        public int Count { get; private set; }

        // Construct the buffer with initial capacity for the given number of
        // elements, but initial logical size 0.
        public FastListOfBoxes(int initial_capacity)
        {
            buf = new T[initial_capacity];
            Count = 0;
        }

        // Set the number of items in the list to 0. The underlying buffer is
        // unchanged and none of the object references are released for the
        // garbage collector.
        //
        // The fact that no references are released by this method is an
        // important difference from System.Collections.Generic.List.
        public void Clear()
        {
            Count = 0;
        }

        // Increase the logical size of the list and return a pointer to the
        // object at the end of the list.
        // 
        // If the underlying buffer is (i) big enough and (ii) already contains
        // a non-null pointer in the relevant location, this pointee is reused
        // and no allocation is performed. If this list is cleared and reused
        // many times, this should be the most likely situation.
        //
        // If the underlying buffer is big enough, but contains a null pointer
        // in the relevant location, then a new object of type T is allocated
        // and default-constructed, a pointer to it is placed in the relevant
        // buffer location, and this pointer is also returned.
        //
        // If the underlying buffer isn't big enough, then the buffer is
        // reallocated with approximately twice the size, and then the steps of
        // the previous paragraph are carried out.
        public T Add()
        {
            if(Count >= buf.Length)
            {
                // We add 1 here in case buf.Length == 0.
                Array.Resize<T>(ref buf, 2*buf.Length + 1);
            }
            if(buf[Count] == null)
            {
                buf[Count] = new T();
            }
            ++Count;
            return buf[Count - 1];
        }

        // Get the element at the given index. Setting the element is not
        // allowed, since that would release a reference to the old pointed-to
        // object, potentially leading to garbage collection.
        //
        // If the contained elements need to be re-initialized, you should add
        // a separate Init method to their class, in addition to its
        // constructor, and then call the Init method on the object returned by
        // the getter below.
        public T this[int index]
        {
            get
            {
#if DEBUG
                if(index >= Count)
                {
                    string msg = String.Format("Index {0} is >= the list " +
                                               "size {1}", index, Count);
                    throw new IndexOutOfRangeException(msg);
                }
#endif
                return buf[index];
            }
        }
    }
}
