using System;
using System.Collections.Generic;

namespace Morpheus
{
    // Storage for the thread-local state in the main loop of
    // DatabaseSearcher.DoSearch.
    public class DatabaseSearcherThreadLocalStorage
    {
        // Counts based on this thread's data. These will be combined at the
        // end of the loop.
        public int num_target_peptides;
        public int num_decoy_peptides;

        // Statistic to keep track of our progress.
        public int proteins;

        // Storage for the peptides encountered in this thread. Each peptide
        // sequence is mapped to true if it is a decoy, and false otherwise.
        public Dictionary<FastSubstring, bool> peptides_observed;

        // Temporary storage for Protein.Digest. Holds the digested peptides of
        // the current protein.
        public FastListOfBoxes<Peptide> digested_peptides_buf;

        // Temporary storage for Peptide.GetVariablyModifiedPeptides. Holds the
        // variably modified peptides generated based on the current peptide.
        public FastListOfBoxes<Peptide> modified_peptides_buf;

        // Temporary storage for the fixed and possibile modifications of the
        // current peptide. Used internally by SetFixedModifications and
        // GetVariablyModifiedPeptides.
        public Dictionary<int, List<Modification>> fixed_modifications_buf;
        public Dictionary<int, List<Modification>> possible_modifications_buf;

        // Temporary storage for GetTandemMassSpectraInMassRange. Holds indices
        // of current spectrum matches.
        public List<int> mass_spectra_indices_buf;

        // Temporary storage for Peptide.Init. Used internally by that method,
        // or rather, by further methods it calls.
        public double[] product_masses_buf;

        // Object with internal state used by Peptide.Init, or rather, by
        // further methods it calls.
        public FastQSorter fast_q_sorter;

        // The match currently under consideration.
        public PeptideSpectrumMatch psm;

        // Holds the best match for each spectrum found so far by this thread.
        // The matches are indexed by the spectrum's number.
        public PeptideSpectrumMatch[] psms;

        public DatabaseSearcherThreadLocalStorage(bool minimize_memory_usage, int psms_length)
        {
            this.num_target_peptides = 0;
            this.num_decoy_peptides = 0;
            this.proteins = 0;

            if (!minimize_memory_usage)
                this.peptides_observed = new Dictionary<FastSubstring, bool>(new LeucineSequenceEqualityComparer());
            else
                this.peptides_observed = null;

            digested_peptides_buf = new FastListOfBoxes<Peptide>(1000);
            modified_peptides_buf = new FastListOfBoxes<Peptide>(1000);
            fixed_modifications_buf = new Dictionary<int, List<Modification>>(1000);
            possible_modifications_buf = new Dictionary<int, List<Modification>>(1000);
            mass_spectra_indices_buf = new List<int>(1000);
            product_masses_buf = new double[1000];
            fast_q_sorter = new FastQSorter();
            psm = new PeptideSpectrumMatch();
            psms = new PeptideSpectrumMatch[psms_length];
        }
    }
}
