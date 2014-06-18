package phylolab.alg.sp;

/* 
 * This is FastS, developed by Siavash Mirarab and Tandy Warnow.
 */

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FastSP {

	public static String VERSION = "1.6.0";

	ReferenceAlignment referenceAlignment;
	int[][] s;
	int k1 = 0, n = 0, refColCount = 0, estColCount;

	long sharedHomologies = 0;
	long totalHomologies = 0;
	long totalHomologiesInEsitmated = 0;
	long correctColumns = 0;
	long effectiveRefColumns = 0;
	long effectiveEstColumns = 0;

	ArrayList<Integer> charsPerEstimatedColumn;

	private PrintStream out = System.out;

	boolean swapped = false;

	Set<Integer> dashChars = new HashSet<Integer>();
	{
		dashChars.add((int) '-');
		dashChars.add((int) '?');
	}

	boolean maskLower = false;
	boolean maskLowerReference = false;

	private long twoChoose(long n) {
		if (n < 2) {
			return 0;
		}
		return n * (n - 1) / 2;
	}

	private String readSequenceName(InputStream in) throws IOException {
		StringBuffer nameStringBuffer = new StringBuffer();
		int ch = 0;
		boolean append = true;
		while (ch != '\n') {
			ch = in.read();
			if (Character.isWhitespace(ch)) {
				append = false;
			}
			if (append) {
				nameStringBuffer.append((char) ch);
			}
			if (Character.isWhitespace(ch) && nameStringBuffer.length() == 0) {
				append = true;
			}
		}
		return nameStringBuffer.toString();
	}

	private int readFirstLine(InputStream in, StringBuffer stringBuffer,
			StringBuffer seqName, boolean toUpper) throws IOException {
		int ch = in.read();
		int colCount = 0, j = 0;
		// skip any headers
		while (ch != '>') {
			in.read();
		}
		seqName.append(readSequenceName(in));
		ch = 0;
		while (ch != '>') {
			ch = in.read();
			int CH = toUpper ? Character.toUpperCase(ch) : ch;
			if (dashChars.contains(CH)) {
				stringBuffer.append((char) ch);
				colCount++;
			} else if ((CH >= 'A' && CH <= 'Z') || (CH >= 'a' && CH <= 'z')) {
				stringBuffer.append((char) ch);
				colCount++;
				j++;
			}

		}
		k1 = Math.max(k1, j);
		return colCount;
	}

	private void readReferenceAlignment(InputStream in, char[] firstSequence)
			throws IOException {

		int ch = '>';
		int i = 0, j = 0;
		int nucInd = 0;
		char[] sequence = firstSequence;

		// Read the rest of sequences, creating a char []
		while (true) {
			if (ch == '>') {
				// Add previous sequence's char [] to the reference alignment
				referenceAlignment.addSequence(sequence);
				sequence = new char[refColCount];
				// update maximum number of (non-gap) sites if necessary.
				k1 = Math.max(k1, nucInd);
				// reset j and nucInd for the new sequence.
				i++;
				nucInd = 0;
				j = 0;
				referenceAlignment.setSequencePosition(readSequenceName(in), i);
			}
			int CH = maskLowerReference ? ch : Character.toUpperCase(ch);

			if (dashChars.contains(CH)) {
				sequence[j] = (char) ch;
				j++;
			} else if ((CH >= 'A' && CH <= 'Z')  || (CH >= 'a' && CH <= 'z')) {
				sequence[j] = (char) ch;
				j++;
				nucInd++;
			} else if (CH == -1) {
				referenceAlignment.addSequence(sequence);
				k1 = Math.max(k1, nucInd);
				break;
			}
			ch = in.read();
		}
		n = i + 1;
	}

	private void readEstimatedAlignment(InputStream in) throws IOException {

		s = new int[n][k1];
		charsPerEstimatedColumn = new ArrayList<Integer>(k1);

		int ch;
		int i = -1, j = 0;
		int nucInd = 0;
		int sequencePositionInReference = -1;
		String name = "";
		while (true) {
			try {
				ch = in.read();
				if (ch == '>') {
					i++;
					// reset j and nucInd for the new sequence.
					nucInd = 0;
					j = 0;
					// Find the position of this sequence in reference alignment
					name = readSequenceName(in);
					sequencePositionInReference = referenceAlignment
							.getSequencePosition(name);
				}
				if (i == 0) {
					charsPerEstimatedColumn.add(0);
				}
				int CH = maskLower ? ch : Character.toUpperCase(ch);
				if (dashChars.contains(CH)) {
					j++;
				} else if (CH >= 'A' && CH <= 'Z') {
					s[sequencePositionInReference][nucInd] = j;
					// increase the number of characters in column j
					charsPerEstimatedColumn.set(j,
							charsPerEstimatedColumn.get(j) + 1);
					j++;
					nucInd++;
				} else if (CH >= 'a' && CH <= 'z') {
					s[sequencePositionInReference][nucInd] = -1;
					j++;
					nucInd++;
				} else if (CH == -1) {
					break;
				}
			} catch (ArrayIndexOutOfBoundsException e) {
				System.err
						.println("Problem reading estimated alignment. Sequence lengths don't match? : "
								+ name);
				throw e;
			}
		}
		long cells = (j + refColCount);
		cells *= n;

		estColCount = j;

		System.err.println("MaxLenNoGap= " + k1 + ", NumSeq= " + n
				+ ", LenRef= " + (swapped ? estColCount : refColCount)
				+ ", LenEst= " + (swapped ? refColCount : estColCount)
				+ ", Cells= " + cells);

		if (k1 <= 0 || n < 2 || refColCount <= 0 || j <= 0 || cells <= 0) {
			System.err
					.println("Error: something is wrong with alignments. Checkout out the statistics above.");
			System.exit(1);
		}
		for (int x : charsPerEstimatedColumn) {
			totalHomologiesInEsitmated += twoChoose(x);
			if (x > 1)
				effectiveEstColumns++;
		}
		// TODO: perform some sanity checks to make sure the alignments are on
		// same data
	}

	private void computeSPFNAndTC() {

		// This array holds current j in N_{i,j}
		// (i.e. non-gap current index) for each row
		int[] refColInd = new int[n];
		// This hash map holds the number of times homologies of
		// the current reference alignment appear in different
		// locations of the estimated alignment (i.e. m(y))
		HashMap<Integer, Integer> estimatedSitesCount;
		// Visit each site in reference alignment
		for (int c = 0; c < refColCount; c++) {
			estimatedSitesCount = new HashMap<Integer, Integer>();
			long refCharCount = 0;
			// For each sequence, if current position is a nuclotide,
			// - find the associated position in the estimated alignment
			// - increment the size of the equivalence class of the associated position
			// - increment the total number of homologies in reference alignment
			for (int i = 0; i < n; i++) {
				int nucltide = referenceAlignment.getNucltide(i, c);
				if (dashChars.contains(nucltide))
					continue;
				if (maskLowerReference && nucltide >= 'a' && nucltide <= 'z') {
					
				} else {
					refCharCount++;
					int y = s[i][refColInd[i]];
					if (y != -1) { // This is the case where the residue was left in
						// an unaligned column in the estimated alignment
						int subjectColumnsCount = estimatedSitesCount
								.containsKey(y) ? estimatedSitesCount.get(y) : 0;
						estimatedSitesCount.put(y, subjectColumnsCount + 1);
					}
				}
				refColInd[i] = refColInd[i] + 1;
			}
			// Update the number of shared and total homologies.
			for (int y : estimatedSitesCount.keySet()) {
				sharedHomologies += twoChoose(estimatedSitesCount.get(y));
			}
			totalHomologies += twoChoose(refCharCount);
			// calculate correctly aligned sites.

			if (estimatedSitesCount.size() == 1 && refCharCount >= 2 && estimatedSitesCount.values().iterator().next() == refCharCount) {
				Integer estColumn = estimatedSitesCount.keySet().iterator().next();
				if (estimatedSitesCount.get(estColumn).equals(
						charsPerEstimatedColumn.get(estColumn))) {
					// System.out.println(estimatedSitesCount.size());
					correctColumns++;
				}
			}
			if (refCharCount >= 2) {
				effectiveRefColumns++;
			}
			// Clear the hash map for the next round
			estimatedSitesCount.clear();
		}
	}

	public double getSPScore() {
		return (sharedHomologies) / (.0 + totalHomologies);
	}

	public double getSPFN() {
		return 1 - getSPScore();
	}

	public double getSPFP() {
		return 1 - getModeler();
	}

	public double getCompressionFactor() {
		return swapped ? refColCount / (estColCount + .0) : estColCount
				/ (refColCount + .0);
	}

	public double getModeler() {
		return (sharedHomologies) / (.0 + totalHomologiesInEsitmated);
	}

	public double getTC() {
		return (correctColumns + 0.) / (effectiveRefColumns + 0.0);
	}

	public void run(String reference, String subject) {
		try {
			File referenceFile = new File(reference);
			System.err.println("Reference alignment: "
					+ referenceFile.getAbsolutePath() + " ...");
			InputStream refIn = new BufferedInputStream(new FileInputStream(
					referenceFile));

			File estimatedFile = new File(subject);
			InputStream subIn = new BufferedInputStream(new FileInputStream(
					estimatedFile));
			System.err.println("Estimated alignment: "
					+ estimatedFile.getAbsolutePath() + " ...");

			StringBuffer refSeq1Name = new StringBuffer();
			StringBuffer refSeq1 = new StringBuffer();
			int firstLineReference = readFirstLine(refIn, refSeq1, refSeq1Name,
					!maskLowerReference);

			StringBuffer estSeq1Name = new StringBuffer();
			StringBuffer estSeq1 = new StringBuffer();
			int firstLineEstimated = readFirstLine(subIn, estSeq1, estSeq1Name,
					!maskLower);

			if ((firstLineEstimated < firstLineReference) && !maskLower &&!maskLowerReference) {
				swapped = true;
				InputStream temp = subIn;
				subIn = refIn;
				refIn = temp;
				refSeq1Name = estSeq1Name;
				refSeq1 = estSeq1;
				refColCount = firstLineEstimated;
			} else {
				refColCount = firstLineReference;
			}

			referenceAlignment = new ReferenceAlignment();
			referenceAlignment.setSequencePosition(refSeq1Name.toString(), 0);

			readReferenceAlignment(refIn, refSeq1.toString().toCharArray());
			refIn.close();

			subIn.close();
			subIn = swapped ? new BufferedInputStream(new FileInputStream(
					referenceFile)) : new BufferedInputStream(
					new FileInputStream(estimatedFile));

			try {

				readEstimatedAlignment(subIn);
			} catch (ArrayIndexOutOfBoundsException e) {
				subIn.close();
				throw e;
			}

			subIn.close();

			System.err.println("computing ...");
			computeSPFNAndTC();

		} catch (FileNotFoundException e) {
			System.err.println(e.getLocalizedMessage());
		} catch (IOException e) {
			System.err.println(e.getLocalizedMessage());
		}

	}

	protected void runWithArguments(String[] args) {
		String estimated = null, reference = null;
		boolean help = false;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-r")) {
				reference = args[i + 1];
				i++;
			}
			if (args[i].equals("-e")) {
				estimated = args[i + 1];
				i++;
			}
			if (args[i].equals("-o")) {
				try {
					out = new PrintStream(args[i + 1]);
					i++;
				} catch (FileNotFoundException e) {
					out = null;
					System.err.println(e.getLocalizedMessage());
				}
			}
			if (args[i].equals("-c")) {
				char[] v = args[i + 1].toCharArray();
				i++;
				for (int j = 0; j < v.length; j++) {
					dashChars.add((int) v[j]);
				}
				System.err.println("Following characters are considered a dash: "
								+ dashChars);
			}
			if (args[i].equals("-ml")) {
				maskLower = true;
				System.err.println("Lower case characters in estimated alignment would not be considered homologs to other characters in the same column.");
			}
			if (args[i].equals("-mlr")) {
				maskLowerReference  = true;	
				System.err.println("Lower case characters in reference alignment would not be considered homologs to other characters in the same column.\n"
						+ "This means that Modeler Score and SPFP are probably meaningless");
			}
			if (args[i].equals("-h")) {
				help = true;
			}
		}

		if (estimated == null || reference == null || out == null || help) {
			System.err
					.println("FastSP version "
							+ this.VERSION
							+ " develped by Siavash Mirarab (smirarab@cs.utexas.edu) and Tandy Warnow at UT-Austin\n\n"
							+ "Usage: FastSP -r reference_alignment_file -e estimated_alignmet_file [-o output_file] [-c GAP_CHARS] [-ml mask out lower case characters in estimated]  [-mlr mask out lower case characters in reference]");
			System.err
					.println("Output: \n"
							+ "	SP-Score:\t number of shared homologies (aligned pairs) / total number of homologies in the reference alignment. \n"
							+ "	Modeler: \t number of shared homologies (aligned pairs) / total number of homologies in the estimated alignment. \n"
							+ "	SP-FN:   \t 1 - SP-Score\n"
							+ "	SP-FP:   \t 1 - Modeler\n"
							+ "	Compression:   \t number of columns in the estimated alignment / number of columns in the reference alignment \n"
							+ "	TC:      \t number of correctly aligned columns / total number of aligned columns in the reference alignment. \n");
			System.exit(1);
		}

		run(reference, estimated);

		if (swapped) {
			long tmp = totalHomologies;
			totalHomologies = totalHomologiesInEsitmated;
			totalHomologiesInEsitmated = tmp;

			tmp = effectiveRefColumns;
			effectiveRefColumns = effectiveEstColumns;
			effectiveEstColumns = tmp;
		}

		System.err.println("Number of shared homologies: " + sharedHomologies);
		System.err.println("Number of homologies in the reference alignment: "
				+ totalHomologies);
		System.err.println("Number of homologies in the estimated alignment: "
				+ totalHomologiesInEsitmated);
		System.err.println("Number of correctly aligned columns: "
				+ correctColumns);
		System.err.println("Number of aligned columns in ref. alignment: "
				+ effectiveRefColumns);
		out.println("SP-Score " + getSPScore());
		out.println("Modeler " + getModeler());
		out.println("SPFN " + getSPFN());
		out.println("SPFP " + getSPFP());
		out.println("Compression " + getCompressionFactor());
		out.println("TC " + getTC());
		out.flush();
	}

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();

		FastSP calc = new FastSP();

		calc.runWithArguments(args);

		System.err.println("Time to compute (seconds): "
				+ (System.currentTimeMillis() - startTime) / 1000.0);
	}

	private class ReferenceAlignment {

		List<char[]> sequencesOfRefAlign = new ArrayList<char[]>();
		HashMap<String, Integer> sequencePositions = new HashMap<String, Integer>();

		void addSequence(char[] s) {
			sequencesOfRefAlign.add(s);
		}

		int getNucltide(int i, int c) {
			return sequencesOfRefAlign.get(i)[c];
		}

		void setSequencePosition(String name, Integer position) {
			sequencePositions.put(name, position);
		}

		int getSequencePosition(String name) {
			if (sequencePositions.containsKey(name)) {
				return sequencePositions.get(name);
			} else {
				throw new RuntimeException("Sequence name Not Found: " + name);
			}
		}
	}
}
