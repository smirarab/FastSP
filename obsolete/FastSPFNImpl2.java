import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;


public class FastSPFNImpl2 {
	
	ReferenceAlignment referenceAlign;
	int [][] s;	
	int k = 0, n = 0;
	
	long sharedHomologies = 0;
	long totalHomologies = 0;
	
 
	private long twoChoose (long n) {
		if (n<2) {return 0;}
		return n*(n-1)/2;
	}
	
	private void readReferenceAlignment(InputStream in) throws IOException {
		
		referenceAlign = new ReferenceAlignment();
		
		int ch;
		int i = -1, j = -1;
		int nucInd = 0;	
		
		while (true) {
			ch = in.read();
			if (ch == '>') {				
				// update maximum number of (non-gap) sites if necessary. 
				if (k < nucInd) {
					k = nucInd;
				}
				// reset j and nucInd for the new sequence.
				i++;
				nucInd = 0;
				j = -1;
				// Read off the name. 
				StringBuffer nameBuffer = new StringBuffer();
				while (ch != '\n') 
				{
					ch = in.read();					
					nameBuffer.append((char)ch);
				}
				referenceAlign.setSequencePosition(nameBuffer.toString(), i);
				
				continue;
			}			
			int CH = Character.toUpperCase(ch);
			if (i == 0 && ( (CH >= 'A' && CH <= 'Z') || (CH == '-') ) ){
				referenceAlign.addColumn(new ArrayList<N>());
			}
			if (CH >= 'A' && CH <= 'Z') {
				j++;
				List<N> column = referenceAlign.getColumn(j);
				column.add(new N (i, nucInd));
				nucInd ++;
			} else if (CH == '-'){
				j++;
			} else if (CH == -1) {
				break;
			}
		}	
		if (k < nucInd) {
			k = nucInd;
		}
		n = i + 1;
		System.err.println("k= " + k + ", n= " + n + ", K= " + referenceAlign.getColumns().size());
	}
	
	private void readSubjectAlignment (InputStream in) throws IOException {
		
		s = new int [n][k];
		
		int ch;
		int i = -1, j = -1;
		int nucInd = 0;		
		int sequencePositionInReference = -1;
		
		while (true) {
			ch = in.read();
			if (ch == '>') {				
				i++;
				// reset j and nucInd for the new sequence.
				nucInd = 0;
				j = -1;
				// Read off the name. 
				StringBuffer nameBuffer = new StringBuffer();
				while (ch != '\n') 
				{
					ch = in.read();					
					nameBuffer.append((char)ch);
				}				
				sequencePositionInReference = referenceAlign.getSequencePosition(nameBuffer.toString());
			}
			int CH = Character.toUpperCase(ch);
			if (CH >= 'A' && CH <= 'Z') {
				j++;		
				s[sequencePositionInReference][nucInd] = j;
				nucInd ++;
			} else if (CH == '-'){
				j++;
			} else if (CH == -1) {
				break;
			}
		}	
	}
	
	private void computeSPFN (){
		for (List<N> refColumn : referenceAlign.getColumns()) {
			HashMap<Integer, Integer> subjectCounts = new HashMap<Integer, Integer>();
			for (N cell: refColumn) {
				int y = s[cell.i][cell.j];
				Integer subjectColumnsCount = 0;
				if (subjectCounts.containsKey(y)) {					
					subjectColumnsCount = subjectCounts.get(y); 
				}
				subjectCounts.put(y, subjectColumnsCount + 1);
			}
			for (int y: subjectCounts.keySet()) {
				int count = subjectCounts.get(y);
				sharedHomologies += twoChoose(count);
			}
			totalHomologies += twoChoose(refColumn.size());	
		}		
	}
	
	public double getSPFN () {
		return (totalHomologies - sharedHomologies)/(.0+totalHomologies);
	}
	
	public void run(String reference, String subject) {
		try {
			MemoryMXBean memorymbean = ManagementFactory.getMemoryMXBean();
			System.err.println("Heap: " + memorymbean.getHeapMemoryUsage().getUsed()/1024);
			
			System.err.println("reading reference alignment " + reference + " ...");
			InputStream refIn = new BufferedInputStream(new FileInputStream(reference));
			readReferenceAlignment(refIn);			
			refIn.close();
			
			System.err.println("Heap: " + memorymbean.getHeapMemoryUsage().getUsed()/1024);
			
			System.err.println("reading subject alignment " + subject + " ...");
			InputStream subIn = new BufferedInputStream(new FileInputStream(subject));
			readSubjectAlignment(subIn);
			subIn.close();
			
			System.err.println("Heap: " + memorymbean.getHeapMemoryUsage().getUsed()/1024);
			
			System.err.println("computing SPFN ...");
			computeSPFN();		
						
			System.err.println("Heap: " + memorymbean.getHeapMemoryUsage().getUsed()/1024);
		} catch (FileNotFoundException e) {
			System.err.println(e.getLocalizedMessage());
		} catch (IOException e) {
			System.err.println(e.getLocalizedMessage());
		}
	}
	
	protected void readArgumentsAndRun(String[] args) {
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
			System.err.println("Usage: FastSPFN -r reference_tree_file -e estimated_tree_file");
			System.exit(1);
		}
		
		run(reference, estimated);
	}
	
	public static void main (String [] args) {
		
		long startTime = System.currentTimeMillis();
		
		FastSPFNImpl2 calc = new FastSPFNImpl2();
				
		calc.readArgumentsAndRun(args);
		
		System.out.println("Time to compute (s): " + 
				(System.currentTimeMillis() - startTime)/1000.0); 
		System.out.println("Number of shared homologies: " + 
				calc.sharedHomologies);
		System.out.println("Total number of homologies in reference alignment:" + 
				calc.totalHomologies);
		System.out.println("SPFN: " + calc.getSPFN());				
	}
	
	private class ReferenceAlignment{
		List<List<N>> sitesOfRefAlign = new ArrayList<List<N>>();
		HashMap<String, Integer> sequencePositions = new HashMap<String, Integer>();		
		List<N> getColumn(int j){
			return sitesOfRefAlign.get(j); 	
		}
		void addColumn(List<N> col) {
			sitesOfRefAlign.add(col);
		}
		List<List<N>> getColumns() {
			return sitesOfRefAlign;
		}
		void setSequencePosition(String name, Integer position){
			sequencePositions.put(name, position);
		}
		int getSequencePosition(String name){
			return sequencePositions.get(name);
		}
	}
	
	private class N {
		int i,j;
		N (int i, int j) {this.i=i; this.j=j;}
	}
}
