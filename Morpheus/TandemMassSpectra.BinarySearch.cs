using System;
using System.Collections.Generic;

namespace Morpheus
{
    public partial class TandemMassSpectra
    {
        // This function looks for spectra in the specified range. The original
        // version of this function returned IEnumerable<TandemMassSpectrum>,
        // and this is still conceptually what we are trying to do. However,
        // the above appears to generate a lot of garbage, because this
        // function occurs inside an inner loop of DatabaseSearcher.DoSearch.
        //
        // Thus, instead, we require the caller to preallocate an array of
        // spectrum indices, "indices", and pass this to us by reference. This
        // function fills in this array (from the beginning) with the indices
        // of the spectra in the specified range. If the array isn't big
        // enough, it is reallocated with double the size. The number of
        // filled-in elements of the array, i.e., the number of spectra in the
        // specified range, is returned. (However, to emphasize, in general,
        // the indices array will be bigger than this.)
        //
        // If precursorMonoisotopicPeakCorrection is false, then pass 
        //   minimumMonoisotopicPeakOffset = 0 and
        //   maximumMonoisotopicPeakOffset = 0.
        public int GetTandemMassSpectraInMassRange(double precursorMass, MassTolerance precursorMassTolerance, 
            int minimumMonoisotopicPeakOffset, int maximumMonoisotopicPeakOffset, ref int[] indices)
        {
            int j = 0;
            for(int i = minimumMonoisotopicPeakOffset; i <= maximumMonoisotopicPeakOffset; i++)
            {
                double precursor_mass = precursorMass + i * Constants.C12_C13_MASS_DIFFERENCE;

                double minimum_precursor_mass = precursor_mass - precursorMassTolerance;
                double maximum_precursor_mass = precursor_mass + precursorMassTolerance;

                int index = BinarySearch(NextHigherDouble(maximum_precursor_mass));
                if(index == Count)
                {
                    index--;
                }
                while(index >= 0 && this[index].PrecursorMass >= minimum_precursor_mass)
                {
                    if(this[index].PrecursorMass <= maximum_precursor_mass)
                    {
                        // Record index as one of the spectra in the range. If
                        // the list of spectrum indices is now not big enough,
                        // resize it.
                        indices[j] = index;
                        j++;
                        if (j >= indices.Length)
                        {
                            System.Array.Resize<int>(ref indices, j*2);
                        }
                    }
                    index--;
                }
            }
            int number_of_indices = j;
            return number_of_indices;
        }

        private int BinarySearch(double precursorMass)
        {
            int low_index = 0;
            int high_index = Count - 1;
            while(low_index <= high_index)
            {
                int mid_index = low_index + ((high_index - low_index) / 2);
                int comparison = this[mid_index].PrecursorMass.CompareTo(precursorMass);
                if(comparison == 0)
                {
                    return mid_index;
                }
                if(comparison < 0)
                {
                    low_index = mid_index + 1;
                }
                else
                {
                    high_index = mid_index - 1;
                }
            }
            return low_index;
        }

        // only works for positive doubles; does not handle special cases
        private static double NextHigherDouble(double value)
        {
            long bits = BitConverter.DoubleToInt64Bits(value);
            long next_double_bits = bits + 1;
            return BitConverter.Int64BitsToDouble(next_double_bits);
        }
    }
}
