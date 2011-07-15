import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;


public class CopyOfFastSP {
	
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
			setSequencePosition(readSequenceName(in), i);
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
				addSequence(sequence);
				sequence = new char [refColCount];
				// update maximum number of (non-gap) sites if necessary. 
				k1 = Math.max(k1, nucInd);
				// reset j and nucInd for the new sequence.
				i++;
				nucInd = 0;
				j = 0;
				setSequencePosition(readSequenceName(in), i);
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
				addSequence(sequence);	
				k1 = Math.max(k1, nucInd);
				break;
			}
			ch = in.read();
		}	
		n = i + 1;		
	}

	
	private void computeSPFN (){
		
		// This array holds current j in N_{i,j} 
		// (i.e. non-gap current index) for each row
		int [] refColInd = new int[n];		
		// This hash map holds the number of times homologies of 
		// the current reference alignment appear in different
		// locations of the estimated alignment (i.e. m(y))
		HashSet <Integer> gapcols = new HashSet<Integer>();	
		// Visit each site in reference alignment
		for (int c = 0; c < refColCount; c++) {		
			long refHomogiesCount = 0;
			// For each sequence, if current position is a nuclotide, 
			// - find the associated position in the estimated alignment
			// - increment the size of the equivalence class of the associated position
			// - increment the total number of homologies in reference alignment
			boolean allgaps = true;
			for (int i = 0; i < n; i++) {
				char nucltide = getNucltide(i , c);
				if (nucltide == '-') {
					continue;
				}
				allgaps = false;				
			}
			if (allgaps) {
				gapcols.add(c);
			}
		}
		for (int i = 0; i < n; i++) {
			String nucltides = getNucltides(i , gapcols);
			System.out.println(">"+getSequencePosition(i));
			System.out.println(nucltides);
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
												
			System.err.println("computing ...");
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
		if (reference == null) {
			System.err.println("Usage: FastSPFN -r reference_alignment_file");
			System.exit(1);
		}
		
		run(reference, estimated);
		
		//System.out.println("SPFN: " + getSPFN());
	}
	
	public static void main (String [] args) {		
		long startTime = System.currentTimeMillis();
		
		CopyOfFastSP calc = new CopyOfFastSP();
				
		calc.runWithArguments(args);
		
		System.err.println("Time to compute (seconds): " + 
				(System.currentTimeMillis() - startTime)/1000.0); 
	}
	 
	List<char[]> sequencesOfRefAlign = new ArrayList<char[]>();
	HashMap<Integer,String> sequencePositions = new HashMap< Integer,String>();		
	
	void addSequence(char[] s) {
		sequencesOfRefAlign.add(s);
	}
	public String getNucltides(int i, HashSet<Integer> gapcols) {
		ArrayList<String> e;
		StringBuffer buffer = new StringBuffer();
		for (int j = 0; j < sequencesOfRefAlign.get(i).length; j++) {
			if (gapcols.contains(j)) {
				continue;
			}
			buffer.append(sequencesOfRefAlign.get(i)[j]);
		}
		                   
		return buffer.toString();
	}
	char getNucltide(int i, int c) {
		return sequencesOfRefAlign.get(i)[c];
	}
	void setSequencePosition(String name, Integer position){
		sequencePositions.put( position, name);
	}
	String getSequencePosition(Integer pos){
		return sequencePositions.get(pos);
	}
	
}
