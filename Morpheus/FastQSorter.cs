using System;

namespace Morpheus
{
    // These classes are copied almost wholesale from Mono's system.Array
    // implementation. The only difference is that the the FastQSortStack array
    // is allocated ahead of time instead of once every time Sort is called.
    //
    // This is useful when Sort is called (i) many times on (ii) small arrays
    // in (iii) multithreaded code. Under such circumstances, the heap
    // allocation of a separate QSortStack array for each call of Sort leads to
    // a lot of locking in the allocator.
    //
    // Note that a different FastQSorter object needs to be used in each
    // thread.
    public class FastQSorter
    {
        private struct FastQSortStack {
            public int high;
            public int low;
        }

        private FastQSortStack[] stack;

        public FastQSorter()
        {
            stack = new FastQSortStack[32];
        }

        private static void swap<T> (T [] array, int i, int j)
        {
            T tmp = array[i];
            array[i] = array[j];
            array[j] = tmp;
        }

		private static bool QSortArrange<T> (T [] array, int lo, int hi, Comparison<T> compare)
		{
			if (array[lo] != null) {
				if (array[hi] == null || compare (array[hi], array[lo]) < 0) {
					swap<T> (array, lo, hi);
					return true;
				}
			}
			
			return false;
		}

		public void Sort<T> (T [] array, int index, int length, Comparison<T> compare)
		{
			int low0 = index;
			int high0 = index + length - 1;

			//FastQSortStack[] stack = new QSortStack[32];
			const int QSORT_THRESHOLD = 7;
			int high, low, mid, i, k;
			int sp = 1;
			T key;
			
			// initialize our stack
			stack[0].high = high0;
			stack[0].low = low0;
			
			do {
				// pop the stack
				sp--;
				high = stack[sp].high;
				low = stack[sp].low;
				
				if ((low + QSORT_THRESHOLD) > high) {
					// switch to insertion sort
					for (i = low + 1; i <= high; i++) {
						for (k = i; k > low; k--) {
							// if keys[k] >= keys[k-1], break
							if (array[k-1] == null)
								break;
							
							if (array[k] != null && compare (array[k], array[k-1]) >= 0)
								break;
							
							swap<T> (array, k - 1, k);
						}
					}
					
					continue;
				}
				
				// calculate the middle element
				mid = low + ((high - low) / 2);
				
				// once we re-order the lo, mid, and hi elements to be in
				// ascending order, we'll use mid as our pivot.
				QSortArrange<T> (array, low, mid, compare);
				if (QSortArrange<T> (array, mid, high, compare))
					QSortArrange<T> (array, low, mid, compare);
				
				key = array[mid];
				
				// since we've already guaranteed that lo <= mid and mid <= hi,
				// we can skip comparing them again
				k = high - 1;
				i = low + 1;
				
				do {
					// Move the walls in
					if (key != null) {
						// find the first element with a value >= pivot value
						while (i < k && compare (key, array[i]) > 0)
							i++;
						
						// find the last element with a value <= pivot value
						while (k > i && compare (key, array[k]) < 0)
							k--;
					} else {
						while (i < k && array[i] == null)
							i++;
						
						while (k > i && array[k] != null)
							k--;
					}
					
					if (k <= i)
						break;
					
					swap<T> (array, i, k);
					
					i++;
					k--;
				} while (true);
				
				// push our partitions onto the stack, largest first
				// (to make sure we don't run out of stack space)
				if ((high - k) >= (k - low)) {
					if ((k + 1) < high) {
						stack[sp].high = high;
						stack[sp].low = k;
						sp++;
					}
					
					if ((k - 1) > low) {
						stack[sp].high = k;
						stack[sp].low = low;
						sp++;
					}
				} else {
					if ((k - 1) > low) {
						stack[sp].high = k;
						stack[sp].low = low;
						sp++;
					}
					
					if ((k + 1) < high) {
						stack[sp].high = high;
						stack[sp].low = k;
						sp++;
					}
				}
			} while (sp > 0);
		}

        // This comparison is faster than the built-in one, but it doesn't take
        // into account NaN and infinite values.
        public static int FastDoubleComparison(double x, double y)
        {
            if (x < y)
                return -1;
            else if (x > y)
                return 1;
            else if (x == y)
                return 0;
            else
                throw new Exception("Incomparable doubles");
        }

		public void Sort(double [] array, int low0, int high0)
		{
            Sort(array, low0, high0, FastDoubleComparison);
        }
    }
}
