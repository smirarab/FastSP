Overview
====
FastSP is a Java program for computing alignment error (SP-FN) quickly and using little memory. It compared two alignments, and computes a bunch of metrics (see below).

Older versions are available [here](http://www.cs.utexas.edu/~phylo/software/fastsp/older) and a change log is available [here](CHANGELOG). 
Dataset

Some test datasets are available [here](http://www.cs.utexas.edu/~phylo/software/fastsp/datasets.zip).

Usage
===

```java -jar FastSP.jar -r reference_alignment_file -e estimated_alignment_file```


FAQ
===
1. Is FastSP sensitive to case? Should I expect different results if I change the alignments from upper case to lower case or vice versa? 

    **Answer:** 
    * By default No. FastSP is by default not sensitive to case. In fact, it is not even sensitive to what characters you have in the alignment (and it doesn't need to). FastSP just cares about whether a certain position in the alignment is a residue or a gap. So, lower case letters are considered aligned as well as upper case case letters. Note that qscore is sensitive to case. qscore treats lower case letters as not aligned.

    * You can add a `-ml` option to make FastSP sensitive to case. `-ml` instructs FastSP that it should ignore any homologies in the estimated alignment where one of both of the characters are lower case. Similarly, `-mlr` instructs FastSP that it should ignore any homologies in the reference alignment.

2. What do I do if I get a `OutOFMemoryException`? 

    **A:** By default Java limits the memory available to programs. If you run out of memory, try increasing the maximum memory available to jvm using the `-Xmx` option. For example, to make 2GB available to jvm use: 

    ```java -Xmx2048m -jar FastSP.jar -r reference_alignment_file -e estimated_alignment_file``` 

    2GB has been more than enough on the largest alignments we have looked at so far (with more than 1,000,000,000 cells.) However, increasing available memory, if you have more memory available, could make FastSP run faster.

3. What is the output? 

    **A:** Run FastSP with a `-h` option to see the output format. The main output is:
	* SP-Score: number of shared homologies (aligned pairs) / total number of homologies in the reference alignment.
	* Modeler: number of shared homologies (aligned pairs) / total number of homologies in the estimated alignment.
	* SP-FN: 1 - SP-Score
	* SP-FP: 1- Modeler
	* TC: number of correctly aligned columns / total number of aligned columns.
	* Compression Factor: number of columns in the estimated alignment / number of columns in the reference alignment
	
	But FastSP also outputs (in standard error):
	* MaxLenNoGap: maximum number of non-gap characters
	* NumSeq: Number of sequences
	* LenRef: Length of reference alignment
	* LenEst: Length of estimated alignment
	* Cells: (LenEst+LenRef)*NumSeq
	* Number of shared homologies
	* Number of homologies in the reference alignment
	* Number of homologies in the estimated alignment
	* Number of correctly aligned columns
	* Number of aligned columns in reference alignment
	
	Make sure you capture and save standard error (using `2>somefilename` if you are interested in these quantities).
	

Publication
====

* FastSP: Linear time calculation of alignment accuracy 
by Siavash Mirarab and Tandy Warnow
Bioinformatics 2011; doi: [10.1093/bioinformatics/btr553](http://bioinformatics.oxfordjournals.org/content/27/23/3250.abstract)

