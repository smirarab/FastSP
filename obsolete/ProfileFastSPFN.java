import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;
import java.util.List;


public class ProfileFastSPFN {
		
	
	public static void main (String [] args) {

		BigDataMatrix.main(args);
				
		List<MemoryPoolMXBean> memoryPoolMXBeans = ManagementFactory.getMemoryPoolMXBeans();
	
		long totalMemory = 0;
		for (MemoryPoolMXBean memoryPoolMXBean : memoryPoolMXBeans) {
			totalMemory += memoryPoolMXBean.getPeakUsage().getUsed();			
		}
		
		System.err.print("Peak Memory (Kb): " + totalMemory/1024);
		System.err.println(", Running time (s): " + ManagementFactory.getRuntimeMXBean().getUptime()/1000.);
		
	}

}
