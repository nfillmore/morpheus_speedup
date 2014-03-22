using System.Collections.Generic;
using System.Text;
using System.Text.RegularExpressions;

namespace Morpheus
{
    public abstract class AminoAcidPolymer
    {
        private static MassType productMassType = MassType.Monoisotopic;

        public FastSubstring BaseSequence;

        public int Length
        {
            get
            {
                return BaseSequence.Length;
            }
        }

        public char this[int index]
        {
            get
            {
                return BaseSequence[index];
            }
        }

        public double MonoisotopicMass
        {
            get
            {
                double monoisotopic_mass = Constants.WATER_MONOISOTOPIC_MASS;

                for(int i = 0; i < BaseSequence.Length; ++i)
                {
                    monoisotopic_mass += AminoAcidMasses.GetMonoisotopicMass(BaseSequence[i]);
                }
                if(fixedModifications != null)
                {
                    foreach(List<Modification> fixed_modifications in fixedModifications.Values)
                    {
                        foreach(Modification fixed_modification in fixed_modifications)
                        {
                            monoisotopic_mass += fixed_modification.MonoisotopicMassShift;
                        }
                    }
                }
                if(variableModifications != null)
                {
                    foreach(Modification variable_modification in variableModifications.Values)
                    {
                        monoisotopic_mass += variable_modification.MonoisotopicMassShift;
                    }
                }

                return monoisotopic_mass;
            }
        }

        public double AverageMass
        {
            get
            {
                double average_mass = Constants.WATER_AVERAGE_MASS;

                for(int i = 0; i < BaseSequence.Length; ++i)
                {
                    average_mass += AminoAcidMasses.GetAverageMass(BaseSequence[i]);
                }
                if(fixedModifications != null)
                {
                    foreach(List<Modification> fixed_modifications in fixedModifications.Values)
                    {
                        foreach(Modification fixed_modification in fixed_modifications)
                        {
                            average_mass += fixed_modification.AverageMassShift;
                        }
                    }
                }
                if(variableModifications != null)
                {
                    foreach(Modification variable_modification in variableModifications.Values)
                    {
                        average_mass += variable_modification.AverageMassShift;
                    }
                }

                return average_mass;
            }
        }

        // Warning: Allocates two new strings. Avoid this in inner loops. Instead use LeucineSequenceEqualityComparer.
        public string BaseLeucineSequence
        {
            get { return BaseSequence.ToString().Replace('I', 'L'); }
        }

        // Warning: Allocates a lot of stuff. Avoid this in inner loops.
        public string Sequence
        {
            get
            {
                StringBuilder sequence = new StringBuilder();

                // fixed modifications on protein N-terminus
                if(fixedModifications != null)
                {
                    List<Modification> prot_n_term_fixed_mods;
                    if(fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                    {
                        foreach(Modification fixed_modification in prot_n_term_fixed_mods)
                        {
                            sequence.Append('[' + fixed_modification.Description + ']');
                        }
                    }
                }
                // variable modification on protein N-terminus
                if(variableModifications != null)
                {
                    Modification prot_n_term_variable_mod;
                    if(variableModifications.TryGetValue(0, out prot_n_term_variable_mod))
                    {
                        sequence.Append('(' + prot_n_term_variable_mod.Description + ')');
                    }
                }

                // fixed modifications on peptide N-terminus
                if(fixedModifications != null)
                {
                    List<Modification> pep_n_term_fixed_mods;
                    if(fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                    {
                        foreach(Modification fixed_modification in pep_n_term_fixed_mods)
                        {
                            sequence.Append('[' + fixed_modification.Description + ']');
                        }
                    }
                }
                // variable modification on peptide N-terminus
                if(variableModifications != null)
                {
                    Modification pep_n_term_variable_mod;
                    if(variableModifications.TryGetValue(1, out pep_n_term_variable_mod))
                    {
                        sequence.Append('(' + pep_n_term_variable_mod.Description + ')');
                    }
                }

                for(int r = 0; r < Length; r++)
                {
                    sequence.Append(this[r]);
                    // fixed modifications on this residue
                    if(fixedModifications != null)
                    {
                        List<Modification> residue_fixed_mods;
                        if(fixedModifications.TryGetValue(r + 2, out residue_fixed_mods))
                        {
                            foreach(Modification fixed_modification in residue_fixed_mods)
                            {
                                sequence.Append('[' + fixed_modification.Description + ']');
                            }
                        }
                    }
                    // variable modification on this residue
                    if(variableModifications != null)
                    {
                        Modification residue_variable_mod;
                        if(variableModifications.TryGetValue(r + 2, out residue_variable_mod))
                        {
                            sequence.Append('(' + residue_variable_mod.Description + ')');
                        }
                    }
                }

                // fixed modifications on peptide C-terminus
                if(fixedModifications != null)
                {
                    List<Modification> pep_c_term_fixed_mods;
                    if(fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                    {
                        foreach(Modification fixed_modification in pep_c_term_fixed_mods)
                        {
                            sequence.Append('[' + fixed_modification.Description + ']');
                        }
                    }
                }
                // variable modification on peptide C-terminus
                if(variableModifications != null)
                {
                    Modification pep_c_term_variable_mod;
                    if(variableModifications.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                    {
                        sequence.Append('(' + pep_c_term_variable_mod.Description + ')');
                    }
                }

                // fixed modifications on protein C-terminus
                if(fixedModifications != null)
                {
                    List<Modification> prot_c_term_fixed_mods;
                    if(fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                    {
                        foreach(Modification fixed_modification in prot_c_term_fixed_mods)
                        {
                            sequence.Append('[' + fixed_modification.Description + ']');
                        }
                    }
                }
                // variable modification on protein C-terminus
                if(variableModifications != null)
                {
                    Modification prot_c_term_variable_mod;
                    if(variableModifications.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                    {
                        sequence.Append('(' + prot_c_term_variable_mod.Description + ')');
                    }
                }

                return sequence.ToString();
            }
        }

        // Warning: Allocates a lot of stuff. Avoid this in inner loops. Instead, cache the sequence, and then use
        // LeucineSequenceEqualityComparer.
        public string LeucineSequence
        {
            get { return Sequence.Replace('I', 'L'); }
        }

        protected static readonly Regex INVALID_AMINO_ACIDS = new Regex("[^ACDEFGHIKLMNPQRSTVWY]", RegexOptions.Compiled);

        // protected AminoAcidPolymer(string baseSequence, bool prevalidated)
        // {
        //     if(prevalidated)
        //     {
        //         BaseSequence = baseSequence;
        //     }
        //     else
        //     {
        //         BaseSequence = INVALID_AMINO_ACIDS.Replace(baseSequence, string.Empty);
        //     }
        // }

        protected void BaseInit(FastSubstring baseSequence)
        {
            BaseSequence = baseSequence;
            initializeProductArrays = true;
            fixedModifications = null;
            variableModifications = null;
            //cumulativeNTerminalMass = null;
            //cumulativeCTerminalMass = null;
        }


        public override string ToString()
        {
            return Sequence;
        }

        public static void SetProductMassType(MassType productMassType)
        {
            AminoAcidPolymer.productMassType = productMassType;
        }

        private bool initializeProductArrays = true;

        private Dictionary<int, List<Modification>> fixedModifications;

        public Dictionary<int, List<Modification>> FixedModifications
        {
            get { return fixedModifications; }
            set
            {
                fixedModifications = value;
                initializeProductArrays = true;
            }
        }

        private Dictionary<int, Modification> variableModifications;

        public Dictionary<int, Modification> VariableModifications
        {
            get { return variableModifications; }
            set
            {
                variableModifications = value;
                initializeProductArrays = true;
            }
        }

        // This method sets the fixed modifications for this polymer, based on
        // the fixedModifications passed as the first argument.
        //
        // The second argument, fixedModificationsBuffer, will be cleared and
        // refilled by this method, and hence should be thread-local storage.
        // The point of doing this is to avoid allocation of a new dictionary
        // each time this method is called, since profiling showed that a large
        // number of dictionaries were allocated by this method. (However,
        // since each dictionary entry is itself a list of modifications, and
        // these are still newly allocated by this method, it isn't clear how
        // much is really saved, in practice, by reusing only the dictionary.
        // Further optimization might be beneficial here.)
        //
        // Note that the fixedModificationsBuffer (saved as
        // this.fixedModifications) is used in other methods in this class, so
        // you won't want to reuse the fixedModificationsBuffer elsewhere until
        // you are completely done processing this AminoAcidPolymer object.
        //
        // The fixedModificationsBuffer is not useful outside of this class
        // (e.g., by the caller of this method).
        public void SetFixedModifications(IEnumerable<Modification> fixedModifications, Dictionary<int, List<Modification>> fixedModificationsBuffer)
        {
            this.fixedModifications = fixedModificationsBuffer;
            this.fixedModifications.Clear();

            foreach(Modification fixed_modification in fixedModifications)
            {
                if(fixed_modification.Type == ModificationType.ProteinNTerminus && (this is Protein ||
                    (this is Peptide && (((Peptide)this).StartResidueNumber == 1 || (((Peptide)this).StartResidueNumber == 2 && ((Peptide)this).Parent[0] == 'M')))) 
                    && (fixed_modification.AminoAcid == char.MinValue || this[0] == fixed_modification.AminoAcid))
                {
                    List<Modification> prot_n_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                    {
                        prot_n_term_fixed_mods = new List<Modification>();
                        prot_n_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(0, prot_n_term_fixed_mods);
                    }
                    else
                    {
                        prot_n_term_fixed_mods.Add(fixed_modification);
                    }
                }

                if(fixed_modification.Type == ModificationType.PeptideNTerminus && (fixed_modification.AminoAcid == char.MinValue || this[0] == fixed_modification.AminoAcid))
                {
                    List<Modification> pep_n_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                    {
                        pep_n_term_fixed_mods = new List<Modification>();
                        pep_n_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(1, pep_n_term_fixed_mods);
                    }
                    else
                    {
                        pep_n_term_fixed_mods.Add(fixed_modification);
                    }
                }

                for(int r = 0; r < Length; r++)
                {
                    if(fixed_modification.Type == ModificationType.AminoAcidResidue && this[r] == fixed_modification.AminoAcid)
                    {
                        List<Modification> residue_fixed_mods;
                        if(!this.fixedModifications.TryGetValue(r + 2, out residue_fixed_mods))
                        {
                            residue_fixed_mods = new List<Modification>();
                            residue_fixed_mods.Add(fixed_modification);
                            this.fixedModifications.Add(r + 2, residue_fixed_mods);
                        }
                        else
                        {
                            residue_fixed_mods.Add(fixed_modification);
                        }
                    }
                }

                if(fixed_modification.Type == ModificationType.PeptideCTerminus && (fixed_modification.AminoAcid == char.MinValue || this[Length - 1] == fixed_modification.AminoAcid))
                {
                    List<Modification> pep_c_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                    {
                        pep_c_term_fixed_mods = new List<Modification>();
                        pep_c_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(Length + 2, pep_c_term_fixed_mods);
                    }
                    else
                    {
                        pep_c_term_fixed_mods.Add(fixed_modification);
                    }
                }

                if(fixed_modification.Type == ModificationType.ProteinCTerminus && (this is Protein || (this is Peptide && ((Peptide)this).EndResidueNumber == ((Peptide)this).Parent.Length - 1) 
                    && (fixed_modification.AminoAcid == char.MinValue || this[Length - 1] == fixed_modification.AminoAcid)))
                {
                    List<Modification> prot_c_term_fixed_mods;
                    if(!this.fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                    {
                        prot_c_term_fixed_mods = new List<Modification>();
                        prot_c_term_fixed_mods.Add(fixed_modification);
                        this.fixedModifications.Add(Length + 3, prot_c_term_fixed_mods);
                    }
                    else
                    {
                        prot_c_term_fixed_mods.Add(fixed_modification);
                    }
                }
            }

            if(this.fixedModifications.Count == 0)
            {
                this.fixedModifications = null;
            }

            initializeProductArrays = true;
        }

        private double[] cumulativeNTerminalMass;
        private double[] cumulativeCTerminalMass;

        private void InitializeProductArrays()
        {
            double mass_shift;

            //cumulativeNTerminalMass = new double[Length];
            if (cumulativeNTerminalMass == null || Length >= cumulativeNTerminalMass.Length)
            {
                System.Array.Resize<double>(ref cumulativeNTerminalMass, Length);
            }
            for(int i = 0; i < Length; ++i)
            {
                cumulativeNTerminalMass[i] = 0.0;
            }

            mass_shift = 0.0;
            // fixed modifications on protein N-terminus
            if(fixedModifications != null)
            {
                List<Modification> prot_n_term_fixed_mods;
                if(fixedModifications.TryGetValue(0, out prot_n_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in prot_n_term_fixed_mods)
                    {
                        mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on the protein N-terminus
            if(variableModifications != null)
            {
                Modification protein_n_term_variable_mod;
                if(variableModifications.TryGetValue(0, out protein_n_term_variable_mod))
                {
                    mass_shift += productMassType == MassType.Average ? protein_n_term_variable_mod.AverageMassShift : protein_n_term_variable_mod.MonoisotopicMassShift;
                }
            }
            // fixed modifications on peptide N-terminus
            if(fixedModifications != null)
            {
                List<Modification> pep_n_term_fixed_mods;
                if(fixedModifications.TryGetValue(1, out pep_n_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in pep_n_term_fixed_mods)
                    {
                        mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on peptide N-terminus
            if(variableModifications != null)
            {
                Modification pep_n_term_variable_mod;
                if(variableModifications.TryGetValue(1, out pep_n_term_variable_mod))
                {
                    mass_shift += productMassType == MassType.Average ? pep_n_term_variable_mod.AverageMassShift : pep_n_term_variable_mod.MonoisotopicMassShift;
                }
            }
            cumulativeNTerminalMass[0] = mass_shift;

            for(int r = 1; r < Length; r++)
            {
                mass_shift = 0.0;
                // fixed modifications on this residue
                if(fixedModifications != null)
                {
                    List<Modification> residue_fixed_mods;
                    if(fixedModifications.TryGetValue(r + 1, out residue_fixed_mods))
                    {
                        foreach(Modification fixed_modification in residue_fixed_mods)
                        {
                            mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                        }
                    }
                }
                // variable modification on this residue
                if(variableModifications != null)
                {
                    Modification residue_variable_mod;
                    if(variableModifications.TryGetValue(r + 1, out residue_variable_mod))
                    {
                        mass_shift += residue_variable_mod.MonoisotopicMassShift;
                    }
                }
                cumulativeNTerminalMass[r] = cumulativeNTerminalMass[r - 1] + (productMassType == MassType.Average ? AminoAcidMasses.GetAverageMass(this[r - 1]) : AminoAcidMasses.GetMonoisotopicMass(this[r - 1])) + mass_shift;
            }

            //cumulativeCTerminalMass = new double[Length];
            if (cumulativeCTerminalMass == null || Length >= cumulativeCTerminalMass.Length)
            {
                System.Array.Resize<double>(ref cumulativeCTerminalMass, Length);
            }
            for(int i = 0; i < Length; ++i)
            {
                cumulativeCTerminalMass[i] = 0.0;
            }

            mass_shift = 0.0;
            // fixed modifications on protein C-terminus
            if(fixedModifications != null)
            {
                List<Modification> prot_c_term_fixed_mods;
                if(fixedModifications.TryGetValue(Length + 3, out prot_c_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in prot_c_term_fixed_mods)
                    {
                        mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on protein C-terminus
            if(variableModifications != null)
            {
                Modification prot_c_term_variable_mod;
                if(variableModifications.TryGetValue(Length + 3, out prot_c_term_variable_mod))
                {
                    mass_shift += prot_c_term_variable_mod.MonoisotopicMassShift;
                }
            }
            // fixed modifications on peptide C-terminus
            if(fixedModifications != null)
            {
                List<Modification> pep_c_term_fixed_mods;
                if(fixedModifications.TryGetValue(Length + 2, out pep_c_term_fixed_mods))
                {
                    foreach(Modification fixed_modification in pep_c_term_fixed_mods)
                    {
                        mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                    }
                }
            }
            // variable modification on peptide C-terminus
            if(variableModifications != null)
            {
                Modification pep_c_term_variable_mod;
                if(variableModifications.TryGetValue(Length + 2, out pep_c_term_variable_mod))
                {
                    mass_shift += productMassType == MassType.Average ? pep_c_term_variable_mod.AverageMassShift : pep_c_term_variable_mod.MonoisotopicMassShift;
                }
            }
            cumulativeCTerminalMass[0] = mass_shift;

            for(int r = 1; r < Length; r++)
            {
                mass_shift = 0.0;
                // fixed modifications on this residue
                if(fixedModifications != null)
                {
                    List<Modification> residue_fixed_mods;
                    if(fixedModifications.TryGetValue(Length - r + 2, out residue_fixed_mods))
                    {
                        foreach(Modification fixed_modification in residue_fixed_mods)
                        {
                            mass_shift += productMassType == MassType.Average ? fixed_modification.AverageMassShift : fixed_modification.MonoisotopicMassShift;
                        }
                    }
                }
                // variable modification on this residue
                if(variableModifications != null)
                {
                    Modification residue_variable_mod;
                    if(variableModifications.TryGetValue(Length - r + 2, out residue_variable_mod))
                    {
                        mass_shift += productMassType == MassType.Average ? residue_variable_mod.AverageMassShift : residue_variable_mod.MonoisotopicMassShift;
                    }
                }

                cumulativeCTerminalMass[r] = cumulativeCTerminalMass[r - 1] + (productMassType == MassType.Average ? AminoAcidMasses.GetAverageMass(this[Length - r]) : AminoAcidMasses.GetMonoisotopicMass(this[Length - r])) + mass_shift;
            }

            initializeProductArrays = false;
        }

        private static readonly ProductCaps PRODUCT_CAPS = ProductCaps.Instance;

        public double CalculateProductMass(ProductType productType, int productNumber)
        {
            if(initializeProductArrays)
            {
                InitializeProductArrays();
            }

            switch(productType)
            {
                case ProductType.b:
                case ProductType.c:
                    return cumulativeNTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType];
                case ProductType.y:
                case ProductType.zdot:
                    return cumulativeCTerminalMass[productNumber] + PRODUCT_CAPS[productType, productMassType];
                default:
                    return double.NaN;
            }
        }

        // This method is called, via PeptideSpectrumMatch.ScoreMatch, in a
        // tight loop in DatabaseSearcher.DoSearch. Originally, (i) a
        // List<Double> containing the product masses was newly allocated and
        // returned each time this method was called, and then (ii)
        // PeptideSpectrumMatch.ScoreMatch converted this List<Double> to an
        // array of doubles. Thus, a lot of garbage was generated and needed to
        // be collected.
        //
        // In the current version of this function, we instead take a
        // preallocated array, product_masses_buf, by reference, and we fill in
        // the beginning of this array with the calculated product masses. If
        // product_masses_buf to small to store the calculated product masses,
        // we resize it so it is big enough. If product_masses_buf is bigger
        // than necessary, we only fill in the beginning of it and leave the
        // rest of the entries unchanged (in case the extra capacity is needed
        // on subsequent calls). We return the total number of entries in the
        // product_masses_buf that are actually used.
        //
        // It is important to emphasize that, in general,
        // product_masses_buf.Length != the number of calculated product
        // masses. That is, the number of elements in the array does not equal
        // the number of elements that are filled in with meaningful info.
        // 
        // Note also that the fast_q_sorter object's internal state will be
        // modified by this method.
        public int CalculateProductMasses(ProductType[] productTypes, ref double[] product_masses_buf, FastQSorter fast_q_sorter)
        {
            // If product_masses isn't big enough to store all the product
            // masses we might want to store, resize it.
            int max_products = 2 * (Length - 1);
            if(product_masses_buf.Length < max_products)
            {
                System.Array.Resize<double>(ref product_masses_buf, max_products);
            }

            // Calculate the product masses and store them in the prefix of
            // product_masses_buf.
            int i = 0;
            for(int r = 1; r < Length; r++)
            {
                for(int p = 0; p < productTypes.Length; p++)
                {
                    if(!(productTypes[p] == ProductType.c && r < Length && this[r] == 'P') &&
                       !(productTypes[p] == ProductType.zdot && Length - r < Length && this[Length - r] == 'P'))
                    {
                        double product_mass = CalculateProductMass(productTypes[p], r);
                        product_masses_buf[i] = product_mass;
                        ++i;
                    }
                }
            }
            int total_products = i;

            // Sort the product masses. Only the filled-in prefix of
            // product_masses_buf is sorted.
            fast_q_sorter.Sort(product_masses_buf, 0, total_products);

            // Return the number of products, i.e., the size of the filled-in
            // prefix of product_masses_buf.
            return total_products;
        }

        protected IEnumerable<Dictionary<int, Modification>> GetVariableModificationPatterns(Dictionary<int, List<Modification>> possibleVariableModifications)
        {
            if(possibleVariableModifications.Count == 0)
            {
                yield return null;
            }
            else
            {
                List<KeyValuePair<int, List<Modification>>> possible_variable_modifications = new List<KeyValuePair<int, List<Modification>>>(possibleVariableModifications);
                int[] base_variable_modification_pattern = new int[Length + 4];
                for(int variable_modifications = 0; variable_modifications <= possible_variable_modifications.Count; variable_modifications++)
                {
                    foreach(int[] variable_modification_pattern in GetVariableModificationPatterns(possible_variable_modifications, possible_variable_modifications.Count - variable_modifications, base_variable_modification_pattern, 0))
                    {
                        yield return GetVariableModificationPattern(variable_modification_pattern, possibleVariableModifications);
                    }
                }
            }
        }

        private static IEnumerable<int[]> GetVariableModificationPatterns(List<KeyValuePair<int, List<Modification>>> possibleVariableModifications, int unmodifiedResiduesDesired, int[] variableModificationPattern, int index)
        {
            if(index < possibleVariableModifications.Count - 1)
            {
                if(unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    foreach(int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired - 1, variableModificationPattern, index + 1))
                    {
                        yield return new_variable_modification_pattern;
                    }
                }
                if(unmodifiedResiduesDesired < possibleVariableModifications.Count - index)
                {
                    for(int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        foreach(int[] new_variable_modification_pattern in GetVariableModificationPatterns(possibleVariableModifications, unmodifiedResiduesDesired, variableModificationPattern, index + 1))
                        {
                            yield return new_variable_modification_pattern;
                        }
                    }
                }
            }
            else
            {
                if(unmodifiedResiduesDesired > 0)
                {
                    variableModificationPattern[possibleVariableModifications[index].Key] = 0;
                    yield return variableModificationPattern;
                }
                else
                {
                    for(int i = 1; i <= possibleVariableModifications[index].Value.Count; i++)
                    {
                        variableModificationPattern[possibleVariableModifications[index].Key] = i;
                        yield return variableModificationPattern;
                    }
                }
            }
        }

        private static Dictionary<int, Modification> GetVariableModificationPattern(int[] variableModificationArray, Dictionary<int, List<Modification>> possibleVariableModifications)
        {
            Dictionary<int, Modification> modification_pattern = new Dictionary<int, Modification>();

            foreach(KeyValuePair<int, List<Modification>> kvp in possibleVariableModifications)
            {
                if(variableModificationArray[kvp.Key] > 0)
                {
                    modification_pattern.Add(kvp.Key, kvp.Value[variableModificationArray[kvp.Key] - 1]);
                }
            }

            return modification_pattern;
        }
    }
}
