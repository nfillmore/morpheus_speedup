using System;

// A FastSubstring object represents the left-closed, right-open substring
// sequence[start, start+length). Initially the object is in an invalid state,
// and one of the Init methods needs to be called to initialize it.
//
// We define Init methods in addition to the constructors because we want to
// reuse buffers of FastSubstring objects. 
//
// We define FastSubstring as a struct instead of a class because it has hardly
// any data and it seems desirable to reduce allocator pressure and memory
// indirections.
public struct FastSubstring
{
    private string sequence;
    private int start;
    private int length;

    // Initialize the substring as a subsequence of a string. Should be
    // identical to the corresponding Init method below.
    public FastSubstring(string sequence, int start, int length)
    {
        this.sequence = sequence;
        this.start = start;
        this.length = length;
    }

    // Initialize the substring as a subsequence of another substring. Should
    // be identical to the corresponding Init method below.
    public FastSubstring(FastSubstring fs, int start, int length)
    {
        this.sequence = fs.sequence;
        this.start = fs.start + start;
        this.length = length;
    }

    // Initialize the substring as a subsequence of a string.
    public void Init(string sequence, int start, int length)
    {
        this.sequence = sequence;
        this.start = start;
        this.length = length;
    }

    // Initialize the substring as a subsequence of another substring.
    public void Init(FastSubstring fs, int start, int length)
    {
        this.sequence = fs.sequence;
        this.start = fs.start + start;
        this.length = length;
    }

    // Return a substring of this string.
    public FastSubstring Substring(int start, int length)
    {
        return new FastSubstring(this, start, length);
    }

    public char this[int index]
    {
        get
        {
            return sequence[index + start];
        }
    }

    public int Length
    {
        get
        {
            return length;
        }
    }

    // Is this subsequence equal to the other string?
    public bool EqualsIgnoringCase(string other)
    {
        if(length != other.Length)
            return false;
        for(int i = 0; i < length; ++i)
        {
            if(Char.ToUpper(this[i]) != Char.ToUpper(other[i]))
            {
                return false;
            }
        }
        return true;
    }

    // This converts the fast substring to an explicit newly allocated string
    // object. Avoid calling this in places where it is important to avoid
    // garbage generation.
    public string ToString()
    {
        return sequence.Substring(start, length);
    }
}
