using System;
using System.Collections.Generic;

// This class implements comparison of FastSubstrings with 'L' and 'I' amino
// acides treated as equivalent. The purpose of doing so is to avoid allocation
// of a new string with 'I' explicitly replaced by 'L', and hence to avoid
// generating excessive garbage, in the main DatabaseSearcher.DoSearch method
// and elsewhere.
class LeucineSequenceEqualityComparer : IEqualityComparer<FastSubstring>
{
    public bool Equals(FastSubstring fs1, FastSubstring fs2)
    {
        if(fs1.Length != fs2.Length)
            return false;
        for(int i = 0; i < fs1.Length; ++i)
        {
            bool equal = (fs1[i] == fs2[i]) ||
                         (fs1[i] == 'I' && fs2[i] == 'L') ||
                         (fs1[i] == 'L' && fs2[i] == 'I');
            if(!equal)
                return false;
        }
        return true;
    }

    // Adapted from http://en.wikipedia.org/wiki/Jenkins_hash_function
    public int GetHashCode(FastSubstring fs)
    {
        Int32 hash, i;
        for(hash = i = 0; i < fs.Length; ++i)
        {
            char c = fs[i];
            if (c == 'I')
                c = 'L';
            hash += Convert.ToInt32(c);
            hash += (hash << 10);
            hash ^= (hash >> 6);
        }
        hash += (hash << 3);
        hash ^= (hash >> 11);
        hash += (hash << 15);
        return hash;
    }
}
