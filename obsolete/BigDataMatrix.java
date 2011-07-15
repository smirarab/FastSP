
/**
 * BigDataMatrix.java
 */

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.math.BigInteger;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Hashtable;

public class BigDataMatrix {

	public String[] Data;
	public String[] Name;
	public int numTaxa;
	
	public int gapOpen, gapExt, mismatchCost;

    protected static final int DIVISION_SCALE = 20;

	BigDataMatrix() {
		
	}
	
	BigDataMatrix(String[] d) {
		Data = d;
		numTaxa = d.length;
	}
	
	public int length() {
		return numTaxa;
	}
	
	public void printMatrix() {
		System.out.println("\nThis is the data matrix:");
		for (int i=0; i<Data.length; i++)
			System.out.println(Name[i] + ": " + "Length= " + Data[i].length() + " Data= " + Data[i]);
		System.out.println();
	}
	
	public void readFasta(String filename) {
		int ch;
		int j;

		initializeMatrix(filename);
		
		try {
			File inputFile = new File(filename);
			BufferedReader in = new BufferedReader(new FileReader(inputFile));

			ch = in.read();
			for (int i=0; i<numTaxa; i++) {
				while (ch != '>') 
					ch = in.read();

				ch = in.read();
				Name[i] = new String();
				while (ch != '\n') {
					Name[i] = Name[i] + (char) ch;
					ch = in.read();
				}
				// now uppercase it
				Name[i] = Name[i].toUpperCase();
				
				//j = 0;
				Data[i] = "";
				StringBuffer Sequence = new StringBuffer();
				ch = in.read();
				while (ch != '>' && ch != -1) {
					if (ch != '\n' && ch != ' ') {
					    Sequence.append((char) ch);
					}
					ch = in.read();
				}
				Data[i] = new String(Sequence.toString().toUpperCase());
				//Data[i] = new String(Sequence, 0, j);

				//System.out.println ("hello " + i + "|" + Name[i] + "|" + Data[i] + "|");
			}
		}
		catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} 
		catch (IOException e) {
			System.err.println(e);
			System.exit(-1);
		}
	}
	
	
	public void initializeMatrix(String filename) {
	    int c;
		
		try {
			File inputFile = new File(filename);
			BufferedReader in = new BufferedReader(new FileReader(inputFile));
			
			numTaxa = 0;
			c = in.read();
			while (c != -1) {
				if (c == '>') numTaxa++;
				c = in.read();
			}
			Data = new String[numTaxa];
			Name = new String[numTaxa];
		}
		catch (FileNotFoundException e) {
			System.err.println("FileStreamsTest: " + e);
		} 
		catch (IOException e) {
			System.err.println("FileStreamsTest: " + e);
		}
	}
	
	public boolean checkEqualLengths() {
		boolean flag = true;
		
		for (int i=0; i<numTaxa; i++)
			for (int j=i+1; j<numTaxa; j++)
				if (Data[i].length() != Data[j].length())
					flag = false;
		return flag;
	}
	
	// Compute the mnhd of the data matrix
	public double mnhd() {
		double mnhd = 0, nhd;
		
		for (int i = 0; i < numTaxa; i++)
		  for (int j = i+1; j < numTaxa; j++) {
		    nhd = hammingDistance (Data[i], Data[j]);
		    if (nhd == -1)
		    	return nhd;
		    else if (nhd > mnhd) {
		      mnhd = nhd;
		    }
		  }
		
		return mnhd;
	}
	
	public double hammingDistance(String A, String B) {
		int k, l = A.length();
		int score = 0, length = 0;
		
		if (A.length() != B.length()) {
			System.err.println("Data file contains sequences of unequal length!");
			return -1;
		}

		for (k = 0; k < l; k++) {
		    if (A.charAt(k) != '-' && B.charAt(k) != '-') {
		    	length++;
		    	if (A.charAt(k) != B.charAt(k)) {
		    		score++;
		    	}
		    }
		}
		if (length == 0)
		    return 0;
		else return (1.0 * score / length);
	}
	
	public double gappiness() {
		if (!checkEqualLengths()) {
			System.err.println("Data file contains sequences of unequal length!");
			return -1;
		}
		int l = Data[0].length();
		int gaps = 0;
		
		for (int n = 0; n < numTaxa; n++)
			for (int k = 0; k < l; k++)
				if (Data[n].charAt(k) == '-')
					gaps++;
		
		double gappiness = 1.0 * gaps / (numTaxa * l); 
		
		return gappiness;
	}
	
	public static BigInteger[] pairwiseSPscore(char[] A1, char[] A2, char[] B1, char[] B2) {
		int i1=0, i2=0;
		BigInteger score = BigInteger.ZERO;
		int A1cnt=0, A2cnt=0, B2cnt=0;
		BigInteger length = BigInteger.ZERO;
		int j;
		
		do {
			A1cnt=1;
			while (A1[i1] == '-' || A2[i1] == '-') {
				if (A1[i1] != '-') A1cnt++;
				if (A2[i1] != '-') A2cnt++;
				i1++;
				if (i1 >= A1.length) break;
			}
			if (i1 >= A1.length) break;
			if (A2[i1] != '-') A2cnt++;
			j = 0;
			while (j < A1cnt) {
				if (B1[i2] != '-') j++;
				if (B2[i2] != '-') B2cnt++;
				i2++;
				if (i2 >= B1.length) break;
			}
			i2--;
			if (i2 >= B1.length) break;
			if (A1[i1] == B1[i2] && A2[i1] == B2[i2] && A2cnt == B2cnt) {
				score = score.add(BigInteger.ONE);
			}
			i1++;
			i2++;
			length = length.add(BigInteger.ONE);
		}
		while (i1 < A1.length && i2 < B1.length);

		BigInteger ans[] = new BigInteger[2];
		ans[0] = score;
		ans[1] = length;

		// testing
		//		System.err.println ("length: " + length + " score: " + score);

		return ans;
	}

	public void SPscore(BigDataMatrix Alignment) {
		int i, j;
		char[] A1, A2, B1, B2;
		int i1, j1;
		BigInteger[] ans = new BigInteger[2];
		BigInteger score = BigInteger.ZERO;
		BigInteger length = BigInteger.ZERO;

		for (i=0; i<numTaxa; i++) {
			for (j=i+1; j<numTaxa; j++) {
				A1 = Data[i].toCharArray();
				A2 = Data[j].toCharArray();
				i1 = Alignment.findName(Name[i]);
				if (i1 == -1)
					System.out.println(Name[i] + " not found in estimated alignment file.");
				j1 = Alignment.findName(Name[j]);
				if (j1 == -1)
					System.out.println(Name[j] + " not found in estimated alignment file.");
				B1 = Alignment.Data[i1].toCharArray();
				B2 = Alignment.Data[j1].toCharArray();
				ans = pairwiseSPscore(A1, A2, B1, B2);
				score = score.add(ans[0]);
				length = length.add(ans[1]);

			}
		}

		// testing
		//		System.err.println ("length: " + length + " score: " + score);

		BigDecimal result = new BigDecimal(length).subtract(new BigDecimal(score));

		result = result.divide(new BigDecimal(length), DIVISION_SCALE, RoundingMode.FLOOR);
		
		//		spscore = (1.0 * length - score) / length;
//		System.out.println("score = " + score + " length = " + length);
		System.out.println(result.toPlainString());
	}

    public static boolean verifyRawSequencesFromAlignments (BigDataMatrix d1, BigDataMatrix d2) {
	if (d1.numTaxa != d2.numTaxa) {
	    return (false);
	}

	// invert d2 names array
	Hashtable<String, Integer> invertMap = new Hashtable<String, Integer>();
	for (int i = 0; i < d2.numTaxa; i++) {
	    invertMap.put(d2.Name[i], new Integer(i));
	}

	for (int i = 0; i < d1.numTaxa; i++) {
	    String name1 = d1.Name[i];
	    String data1 = d1.Data[i];
	    Integer jInteger = invertMap.get(name1);
	    if (jInteger == null) {
		return (false);
	    }
	    int j = jInteger.intValue();
	    String data2 = d2.Data[j];

	    // testing
	    //System.out.println (data2);

	    String raw1 = data1.replaceAll("\\-", "");
	    String raw2 = data2.replaceAll("\\-", "");

	    // testing
	    //System.out.println (raw1 + " " + raw2);

	    // drop all indels, remaining raw sequences must match
	    if (!raw1.equals(raw2)) {
		return (false);
	    }
	}

	return (true);
    }


	public int findName(String name) {
		boolean found = false;
		int i=0;

		while (i < numTaxa && !found) {
			if (Name[i].equals(name))
				found = true;
			else i++;
		}
		if (found) return i;
		else return -1;
	}

	public static void main(String[] args) {
		boolean computeMNHD = false, computeGappiness = false, computeSPScore = false;
		String filename = "", filename2 = "";
		
/*		String a1 = "ACATA--";
		String a2 = "A--TACC";
		String b1 = "AC-ATA--";
		String b2 = "A-TA--CC";
		
		pairwiseSPscore(a1.toCharArray(), a2.toCharArray(), b1.toCharArray(), b2.toCharArray());
*/		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-f")) {
				filename = args[i+1];
			}
			if (args[i].equals("-v")) {
				filename2 = args[i+1];
			}
			if (args[i].equals("-mnhd"))
				computeMNHD = true;
			if (args[i].equals("-gappiness"))
				computeGappiness = true;
			if (args[i].equals("-sp"))
				computeSPScore = true;
		}
		
		if (computeMNHD) {
			BigDataMatrix DM = new BigDataMatrix();
			DM.readFasta(filename);			
			System.out.println(DM.mnhd());
		}
		if (computeGappiness) {
			BigDataMatrix DM = new BigDataMatrix();
			DM.readFasta(filename);			
			System.out.println(DM.gappiness());
		}
		if (computeSPScore) {
			BigDataMatrix DM1 = new BigDataMatrix(), DM2 = new BigDataMatrix();
			DM1.readFasta(filename);
			//System.out.println ("dm1");
			DM2.readFasta(filename2);
			//System.out.println ("dm2");
			if (!verifyRawSequencesFromAlignments(DM1, DM2)) {
			    System.err.println ("ERROR: could not verify that raw sequences match from the two input alignments! Aborting.\n");
			    System.exit(1);
			}
			DM2.SPscore(DM1);
		}

	}
	
}
