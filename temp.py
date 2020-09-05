
import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import Counter
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')

def main():
    """A Simple Bioinformatics App"""
    st.title('SeExplorer: Sequence Explorer')
    menus = ['Intro', 'DNA Sequence', 'Dotplot', 'About']
    choices = st.sidebar.selectbox("Select Activity", menus)

    if choices == 'Intro':
        st.subheader('Introduction to Bioinformatics')
    elif choices == 'DNA Sequence':
        st.subheader('DNA Sequence Analysis')
        seqfile = st.file_uploader("Upload a FASTA File", type=['fasta', 'fa', 'txt', ])
        if seqfile is not None:
            dna_record = SeqIO.read(seqfile, "fasta")
            # st.write(dna_record)
            seq = dna_record.seq

            details = st.radio('Details', ('Description', 'Sequence'))
            if details == 'Description':
                st.write(dna_record.description)
            elif details == 'Sequence':
                st.write(dna_record.seq)


            # Nucleotide Frequency
            st.subheader('Nucleotide Frequency')
            dna_freq = Counter(seq)
            st.write(dna_freq)

            # Color Picker
            adenine_color = st.beta_color_picker('Adenine Color')
            guanine_color = st.beta_color_picker('Guanine Color')
            thymine_color = st.beta_color_picker('Thymine Color')
            cytosil_color = st.beta_color_picker('Cytosil Color')

            # Frequency Plotting
            if st.button('Plot Frequency'):
                barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(guanine_color)
                barlist[2].set_color(thymine_color)
                barlist[3].set_color(cytosil_color)
                st.pyplot()

            st.subheader('DNA Composition')
            gc = GC(seq)
            st.write("gc)




    elif choices == 'DotPlot':
        st.subheader('Generate DotPlot for Two Sequences')
    elif choices == 'About':
        st.subheader('About App')









if __name__=='__main__':
        main()
