using System.Collections.Generic;

namespace Morpheus
{
    public class Peptide : AminoAcidPolymer
    {
        public Protein Parent { get; private set; }

        public int StartResidueNumber { get; private set; }

        public int EndResidueNumber { get; private set; }

        public int MissedCleavages { get; private set; }

        public char PreviousAminoAcid { get; private set; }

        public char NextAminoAcid { get; private set; }

        public bool Decoy
        {
            get { return Parent.Decoy; }
        }

        public bool Target
        {
            get { return !Decoy; }
        }

        // this doesn't need to be efficient
        public string ExtendedSequence
        {
            get
            {
                return PreviousAminoAcid.ToString() + '.' + Sequence + '.' + NextAminoAcid.ToString();
            }
        }

        // this isn't used anywhere
        public string ExtendedLeucineSequence
        {
            get
            {
                return ExtendedSequence.Replace('I', 'L');
            }
        }

        // public Peptide(Protein parent, int startResidueNumber, int endResidueNumber, int missedCleavages)
        //     : base(parent.BaseSequence.Substring(startResidueNumber - 1, endResidueNumber - startResidueNumber + 1), true)
        // {
        //     Parent = parent;
        //     StartResidueNumber = startResidueNumber;
        //     EndResidueNumber = endResidueNumber;
        //     MissedCleavages = missedCleavages;
        //     if(startResidueNumber - 1 - 1 >= 0)
        //     {
        //         PreviousAminoAcid = parent[startResidueNumber - 1 - 1];
        //     }
        //     else
        //     {
        //         PreviousAminoAcid = '-';
        //     }
        //     if(endResidueNumber - 1 + 1 < parent.Length)
        //     {
        //         NextAminoAcid = parent[endResidueNumber - 1 + 1];
        //     }
        //     else
        //     {
        //         NextAminoAcid = '-';
        //     }
        // }

        // private Peptide(Peptide peptide) : this(peptide.Parent, peptide.StartResidueNumber, peptide.EndResidueNumber, peptide.MissedCleavages) { }

        public void Init(Protein parent, int startResidueNumber, int endResidueNumber, int missedCleavages)
        {
            BaseInit(new FastSubstring(parent.BaseSequence, startResidueNumber - 1, endResidueNumber - startResidueNumber + 1));
            Parent = parent;
            StartResidueNumber = startResidueNumber;
            EndResidueNumber = endResidueNumber;
            MissedCleavages = missedCleavages;
            if(startResidueNumber - 1 - 1 >= 0)
            {
                PreviousAminoAcid = parent[startResidueNumber - 1 - 1];
            }
            else
            {
                PreviousAminoAcid = '-';
            }
            if(endResidueNumber - 1 + 1 < parent.Length)
            {
                NextAminoAcid = parent[endResidueNumber - 1 + 1];
            }
            else
            {
                NextAminoAcid = '-';
            }
        }

        public void CopyFrom(Peptide peptide)
        {
            Init(peptide.Parent, peptide.StartResidueNumber, peptide.EndResidueNumber, peptide.MissedCleavages);
        }

        public int GetVariablyModifiedPeptides(IEnumerable<Modification> variableModifications, int maximumVariableModificationIsoforms, ref Peptide[] peptides)
        {
            int p = 0;
            Dictionary<int, List<Modification>> possible_modifications = new Dictionary<int, List<Modification>>(Length + 4);

            foreach(Modification variable_modification in variableModifications)
            {
                if(variable_modification.Type == ModificationType.ProteinNTerminus && (StartResidueNumber == 1 || (StartResidueNumber == 2 && Parent[0] == 'M')) 
                    && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    List<Modification> prot_n_term_variable_mods;
                    if(!possible_modifications.TryGetValue(0, out prot_n_term_variable_mods))
                    {
                        prot_n_term_variable_mods = new List<Modification>();
                        prot_n_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(0, prot_n_term_variable_mods);
                    }
                    else
                    {
                        prot_n_term_variable_mods.Add(variable_modification);
                    }
                }

                if(variable_modification.Type == ModificationType.PeptideNTerminus && (variable_modification.AminoAcid == char.MinValue || this[0] == variable_modification.AminoAcid))
                {
                    List<Modification> pep_n_term_variable_mods;
                    if(!possible_modifications.TryGetValue(1, out pep_n_term_variable_mods))
                    {
                        pep_n_term_variable_mods = new List<Modification>();
                        pep_n_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(1, pep_n_term_variable_mods);
                    }
                    else
                    {
                        pep_n_term_variable_mods.Add(variable_modification);
                    }
                }

                for(int r = 0; r < Length; r++)
                {
                    if(variable_modification.Type == ModificationType.AminoAcidResidue && this[r] == variable_modification.AminoAcid)
                    {
                        List<Modification> residue_variable_mods;
                        if(!possible_modifications.TryGetValue(r + 2, out residue_variable_mods))
                        {
                            residue_variable_mods = new List<Modification>();
                            residue_variable_mods.Add(variable_modification);
                            possible_modifications.Add(r + 2, residue_variable_mods);
                        }
                        else
                        {
                            residue_variable_mods.Add(variable_modification);
                        }
                    }
                }

                if(variable_modification.Type == ModificationType.PeptideCTerminus && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    List<Modification> pep_c_term_variable_mods;
                    if(!possible_modifications.TryGetValue(Length + 2, out pep_c_term_variable_mods))
                    {
                        pep_c_term_variable_mods = new List<Modification>();
                        pep_c_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(Length + 2, pep_c_term_variable_mods);
                    }
                    else
                    {
                        pep_c_term_variable_mods.Add(variable_modification);
                    }
                }

                if(variable_modification.Type == ModificationType.ProteinCTerminus && (EndResidueNumber == Parent.Length - 1) 
                    && (variable_modification.AminoAcid == char.MinValue || this[Length - 1] == variable_modification.AminoAcid))
                {
                    List<Modification> prot_c_term_variable_mods;
                    if(!possible_modifications.TryGetValue(Length + 3, out prot_c_term_variable_mods))
                    {
                        prot_c_term_variable_mods = new List<Modification>();
                        prot_c_term_variable_mods.Add(variable_modification);
                        possible_modifications.Add(Length + 3, prot_c_term_variable_mods);
                    }
                    else
                    {
                        prot_c_term_variable_mods.Add(variable_modification);
                    }
                }
            }

            int variable_modification_isoforms = 0;
            foreach(Dictionary<int, Modification> kvp in GetVariableModificationPatterns(possible_modifications))
            {
                // Peptide peptide = new Peptide(this);
                // peptide.FixedModifications = FixedModifications;
                // peptide.VariableModifications = kvp;
                // yield return peptide;
                ReallocPeptideBuf(ref peptides, p);
                peptides[p].CopyFrom(this);
                peptides[p].FixedModifications = FixedModifications;
                peptides[p].VariableModifications = kvp;
                ++p;
                variable_modification_isoforms++;
                if(variable_modification_isoforms == maximumVariableModificationIsoforms)
                {
                    //yield break;
                    break;
                }
            }
            int num_peptides = p;
            return p;
        }

        public static void ReallocPeptideBuf(ref Peptide[] peptides, int p)
        {
            int old_length = (peptides == null) ? 0 : peptides.Length;
            if (p >= old_length)
            {
                int new_length = 2*p;
				System.Array.Resize<Peptide>(ref peptides, new_length);
                for(int q = old_length; q < new_length; ++q)
                {
                    peptides[q] = new Peptide();
                }
            }
        }
    }
}
