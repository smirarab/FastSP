package phylolab.alg.sp;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public class FastSP {
	
	ReferenceAlignment referenceAlign;
	int [][] s;	
	int k1 = 0, k2 = 0, n = 0, refColCount = 0;
	
	long sharedHomologies = 0;
	long totalHomologies = 0;
	
 
	private long twoChoose (long n) {
		if (n<2) {return 0;}
		return n*(n-1)/2;
	}
	
	private String readSequenceName(InputStream in) throws IOException {
		StringBuffer nameStringBuffer = new StringBuffer();
		int ch = 0;
		boolean append = true;
		while (ch != '\n') 
		{
			ch = in.read();		
			if (Character.isWhitespace(ch)) {
				append = false;
			}
			if (append) {
				nameStringBuffer.append((char)ch);
			}
			if (Character.isWhitespace(ch) && nameStringBuffer.length() == 0) {
				append = true;
			}
		}		
		return nameStringBuffer.toString();
	}

	private void readReferenceAlignment(InputStream in) throws IOException {
		
		referenceAlign = new ReferenceAlignment();
		
		int ch;
		int i = 0, j = 0;
		int nucInd = 0;	
		char[] sequence;
		
		//Read the first sequence
		{
			StringBuffer stringBuffer = new StringBuffer(); 			
			ch = in.read();
			//skip any headers
			while (ch != '>') { in.read();}
			referenceAlign.setSequencePosition(readSequenceName(in), i);
			stringBuffer = new StringBuffer();	
			ch = 0;
			while (ch != '>') 
			{
				ch = in.read();							
				int CH = Character.toUpperCase(ch);
				if (CH >= 'A' && CH <= 'Z') {
					k1 ++;
					stringBuffer.append((char)ch);
				} else if (CH == '-'){
					stringBuffer.append((char)ch);
				} 
			}	
			refColCount = stringBuffer.length();
			sequence = stringBuffer.toString().toCharArray();
		}
		// Read the rest of sequences, directly creating a char []		
		while (true) {
			if (ch == '>') {
				// Add previous sequence's char [] to the reference alignment 
				referenceAlign.addSequence(sequence);
				sequence = new char [refColCount];
				// update maximum number of (non-gap) sites if necessary. 
				k1 = Math.max(k1, nucInd);
				// reset j and nucInd for the new sequence.
				i++;
				nucInd = 0;
				j = 0;
				referenceAlign.setSequencePosition(readSequenceName(in), i);
			}			
			int CH = Character.toUpperCase(ch);

			if ((CH >= 'A' && CH <= 'Z')) {				
				sequence[j]=(char)ch;
				j++;
				nucInd ++;
			} else if (CH == '-'){
				sequence[j]=(char)ch;
				j++;
			} else if (CH == -1) {
				referenceAlign.addSequence(sequence);	
				k1 = Math.max(k1, nucInd);
				break;
			}
			ch = in.read();
		}	
		n = i + 1;		
	}

	private void readSubjectAlignment (InputStream in) throws IOException {
		
		s = new int [n][k1];
		
		int ch;
		int i = -1, j = 0;
		int nucInd = 0;		
		int sequencePositionInReference = -1;
		
		while (true) {
			ch = in.read();
			if (ch == '>') {				
				i++;
				// reset j and nucInd for the new sequence.
				nucInd = 0;
				j = 0;			
				// Find the position of this sequence in reference alignment
				sequencePositionInReference = referenceAlign.getSequencePosition(readSequenceName(in));
			}
			int CH = Character.toUpperCase(ch);
			if (CH >= 'A' && CH <= 'Z') {
				s[sequencePositionInReference][nucInd] = j;
				j++;	
				nucInd ++;
			} else if (CH == '-'){
				j++;
			} else if (CH == -1) {
				break;
			}
		}	
		long cells = (j+refColCount)*n;
		System.err.println("k= " + k1 + ", n= " + n + 
				", K1= " + refColCount + ", K2= " + j + 
				", cells= " + cells);
	}
	
	private void computeSPFN (){
		
		// This array holds current j in N_{i,j} 
		// (i.e. non-gap current index) for each row
		int [] refColInd = new int[n];		
		// This hash map holds the number of times homologies of 
		// the current reference alignment appear in different
		// locations of the estimated alignment (i.e. m(y))
		HashMap<Integer, Integer> estimatedSitesCount;	
		// Visit each site in reference alignment
		for (int c = 0; c < refColCount; c++) {		
			estimatedSitesCount = new HashMap<Integer, Integer>();	
			long refHomogiesCount = 0;
			// For each sequence, if current position is a nuclotide, 
			// - find the associated position in the estimated alignment
			// - increment the size of the equivalence class of the associated position
			// - increment the total number of homologies in reference alignment
			for (int i = 0; i < n; i++) {
				char nucltide = referenceAlign.getNucltide(i , c);
				if (nucltide == '-')
					continue;
				refHomogiesCount++;
				int y = s[i][refColInd[i]];
				refColInd[i] = refColInd[i] + 1;
				int subjectColumnsCount = 0;
				if (estimatedSitesCount.containsKey(y)) {					
					subjectColumnsCount = estimatedSitesCount.get(y); 
				}
				estimatedSitesCount.put(y, subjectColumnsCount + 1);
			}
			// Update the number of shared and total homologies.
			for (int y: estimatedSitesCount.keySet()) {
				sharedHomologies += twoChoose(estimatedSitesCount.get(y));
			}
			totalHomologies += twoChoose(refHomogiesCount);	
			// Clert the hash map for the next round
			estimatedSitesCount.clear();
		}		
	}
	
	public double getSPFN () {
		return (totalHomologies - sharedHomologies)/(.0+totalHomologies);
	}
	
	public void run(String reference, String subject) {
		try {
			File referenceFile = new File(reference);
			System.err.println("reading reference alignment " + referenceFile.getAbsolutePath() + " ...");
			InputStream refIn = new BufferedInputStream(new FileInputStream(referenceFile));
			readReferenceAlignment(refIn);			
			refIn.close();
						
			File subjectFile = new File(subject);
			System.err.println("reading subject alignment " + subjectFile.getAbsolutePath() + " ...");
			InputStream subIn = new BufferedInputStream(new FileInputStream(subjectFile));
			readSubjectAlignment(subIn);
			subIn.close();
						
			System.err.println("computing SPFN ...");
			computeSPFN();		
									
		} catch (FileNotFoundException e) {
			System.err.println(e.getLocalizedMessage());
		} catch (IOException e) {
			System.err.println(e.getLocalizedMessage());
		}
	}
	
	protected void runWithArguments(String[] args) {
		String estimated = null, reference = null;		
		for (int i=0; i<args.length; i++) {
			if (args[i].equals("-r")) {
				reference = args[i+1];
			}
			if (args[i].equals("-e")) {
				estimated = args[i+1];
			}
		}
		if (estimated == null || reference == null) {
			System.err.println("Usage: FastSP -r reference_alignment_file -e estimated_alignmet_file");
			System.exit(1);
		}
		
		run(reference, estimated);
		
		System.err.println("Number of shared homologies: " + sharedHomologies);
		System.err.println("Number of homologies in the reference alignment: " + 
				totalHomologies);
		System.out.println("SPFN: " + getSPFN());
	}
	
	public static void main (String [] args) {		
		long startTime = System.currentTimeMillis();
		
		FastSP calc = new FastSP();
				
		calc.runWithArguments(args);
		
		System.err.println("Time to compute (seconds): " + 
				(System.currentTimeMillis() - startTime)/1000.0); 
	}
	 
	private class ReferenceAlignment{
		
		List<char[]> sequencesOfRefAlign = new ArrayList<char[]>();
		HashMap<String, Integer> sequencePositions = new HashMap<String, Integer>();		
		
		void addSequence(char[] s) {
			sequencesOfRefAlign.add(s);
		}
		char getNucltide(int i, int c) {
			return sequencesOfRefAlign.get(i)[c];
		} 
		void setSequencePosition(String name, Integer position){
			sequencePositions.put(name, position);
		}
		int getSequencePosition(String name){
			return sequencePositions.get(name);
		}
	}	
}
