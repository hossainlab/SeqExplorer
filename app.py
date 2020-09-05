import streamlit as st
# Sequence Processing Packages
from collections import Counter
from skbio import Sequence
from skbio.sequence.distance import hamming
from Bio import SeqIO
from Bio.SeqUtils import GC
import pandas as pd
import numpy as np
from PIL import Image
# Data Visualization Packages
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('whitegrid')
# set figuresize and fontsize
plt.rcParams['figure.figsize'] = (10,6)
plt.rcParams['font.size'] = 18



def main():
    """A simple bioinformatics app"""
    st.title("SeqExplorer: Explore DNA and Protein Sequence")
    st.sidebar.title("App Options")

    # Analysis Functions
    def ntFrequency(seq, useall=False, calpc=False):
        """Count the frequencies of each bases in sequence including every letter"""
        length = len(seq)
        if calpc:
            # Make a dictionary "freqs" to contain the frequency(in % ) of each base.
            freqs = {}
        else:
            # Make a dictionary "base_counts" to contain the frequency(whole number) of each base.
            base_counts = {}
        if useall:
            # If we want to look at every letter that appears in the sequence.
            seqset = set(seq)
        else:
            # If we just want to look at the four bases A, T, C, G
            seqset = ("A", "T", "G", "C")

        for letter in seqset:
            num = seq.count(letter)
            if calpc:
                # The frequency is calculated out of the total sequence length, even though some bases are not A, T, C, G
                freq = round(num/length, 2)
                freqs[letter] = freq
            else:
                # Contain the actual number of bases.
                base_counts[letter] = num
        if calpc:
            return freqs
        else:
            return base_counts

    # Seqinfo Function
    def head(seq):
        """Returns frist 1000 bp of a sequence"""
        h = seq[0:1000]
        return h

    def tail(seq):
        """Returns last 1000 bp of a sequence"""
        t = seq[-1000:]
        return t


    # menus
    menu = ["Intro", "DNA Sequence", "Protein Sequence", "About"]
    choices = st.sidebar.selectbox("Activity", menu)

    if choices == 'Intro':
        st.subheader("Introduction to Bioinformatics")
        st.markdown("Bioinformatics is an interdisciplinary field that develops methods and software tools for understanding biological data, in particular when the data sets are large and complex.")
        # img = Image.open("./img/dna.jpg")
        # img.thumbnail((300,800))
        # st.image(img)
        st.subheader("Features")
        st.markdown("**1. DNA Sequence Statitstics**")
        st.markdown("**2. DNA Sequence Operations**")
        st.markdown("**3. Protein Sequence Analysis**")
        st.markdown("**4. Data Visualization**")

    # DNA Sequence Analysis
    elif choices == 'DNA Sequence':
        st.subheader("DNA Sequence Analysis")
        # sequence file upload
        seq_file = st.file_uploader("Upload a File", type=["fasta", "fa"])
        if seq_file is not None:
            file = SeqIO.read(seq_file, "fasta")
            seq = file.seq
            # seqinfo
            st.sidebar.subheader("Sequence Statistics")
            details = st.sidebar.radio("Seqinfo", ("Description","Head", "Tail", "Entire Sequence"))
            if details == 'Description':
                st.write(file.description)
            elif details == 'Head':
                st.write("Frist 1000 Base-pair(bp) of Sequence.\n", head(seq))
            elif details == 'Tail':
                st.write("Last 1000 Base-pair(bp) of Sequence.\n", tail(seq))
            elif details == 'Entire Sequence':
                st.write("Take a Look at Entire Eequence.\n", file.seq)

        # nucleotide frequency
            freq_analysis = st.sidebar.radio("Nucleotide Analysis", ("Frequency(Only ATGC)","Frequency(All)", "Percentage(Only ATGC)","Percentage(All)", "Frequency Plot(Bar)", "Frequency Plot(Pie)"))
            if freq_analysis == 'Frequency(Only ATGC)':
                freq = ntFrequency(seq)
                st.write(freq)
            if freq_analysis == 'Frequency(All)':
                freq = ntFrequency(seq, useall=True)
                st.write(freq)
            elif  freq_analysis == "Percentage(Only ATGC)":
                freq = ntFrequency(seq, calpc=True)
                st.write(freq)
            elif  freq_analysis == "Percentage(All)":
                freq = ntFrequency(seq, calpc=True, useall=True)
                st.write(freq)
            elif freq_analysis == 'Frequency Plot(Bar)':
                freq = ntFrequency(seq)
                plt.bar(freq.keys(), freq.values())
                plt.xlabel("Bases")
                plt.ylabel("Frequency")
                plt.title("Nucleotide Frquency Distribution")
                plt.tight_layout()
                st.pyplot()

            elif freq_analysis == 'Frequency Plot(Pie)':
                freq = ntFrequency(seq)
                plt.pie(freq.values(), labels=freq.keys(), autopct='%1.1f%%', shadow=True)
                plt.tight_layout()
                st.pyplot()



                        # DNA Composition
            def calculateGC(seq):
                """Take DNA sequence as input and calculate the AT content."""
                no_of_a = seq.count("G")
                no_of_t = seq.count("C")
                total = no_of_a + no_of_t
                gc = round(total/len(seq) * 100, 2)
                return gc

            def calculateAT(seq):
                """Take DNA sequence as input and calculate the AT content."""
                no_of_a = seq.count("A")
                no_of_t = seq.count("T")
                total = no_of_a + no_of_t
                at = round(total/len(seq) * 100, 2)
                return at

            dna_com = st.sidebar.radio("DNA Composition", ("Length", "GC Ratio", "AT Ratio"))
            if  dna_com == "Length":
                st.write("Length of DNA Sequence = ", len(seq))
            elif dna_com == 'GC Ratio':
                st.write("GC Content = ", calculateGC(seq))
            elif dna_com == 'AT Ratio':
                st.write("AT Content = ", calculateAT(seq))

            # Sequence Operations
            st.sidebar.subheader("Sequence Operations")
            if st.sidebar.checkbox("Complement"):
                st.write("Complementary Sequence\n", seq.complement())
            elif st.sidebar.checkbox("Reverse Complement"):
                st.write("Reverse Complement\n", seq.reverse_complement())
            elif st.sidebar.checkbox("Transcription"):
                st.write("Transcription\n", seq.transcribe())
            elif st.sidebar.checkbox("Back Transcription"):
                st.write("Back Transcription\n", seq.back_transcribe())
            elif st.sidebar.checkbox("Translation"):
                st.write("Translation\n", seq.translate())









    elif choices == 'Protein Sequence':
        st.subheader('Protein Sequence Analysis')
        st.markdown("### Comming Soon....")
        st.markdown("1. Protein Sequence Analysis")
        st.markdown("2. Protein Visualization")

    elif choices == 'About':
        st.subheader('About the Developer')
        st.markdown("I am a student of Microbiology at Jagannath University and Health Data Science Enthusiastic. My research interests are health data science, bioinformatics, microscopic image analysis, and machine learning. Recently, we have formed an organization called Health Data Research Organization for creating research opportunities in health data science and genomic data science in our country and I am working on it")

    else:
        print("Select an activity")














if __name__ == '__main__':
    main()
