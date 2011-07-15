import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.util.List;


public class ProfileBigDataMatrix {
		
	
	public static void main (String [] args) {

		BigDataMatrix.main(args);
		
		List<MemoryPoolMXBean> memoryPoolMXBeans = ManagementFactory.getMemoryPoolMXBeans();
	
		long totalMemory = 0;
		for (MemoryPoolMXBean memoryPoolMXBean : memoryPoolMXBeans) {
			totalMemory += memoryPoolMXBean.getPeakUsage().getUsed();			
		}
		
		System.err.print("Peak (k): " + totalMemory/1024);
		System.err.println(", Time (milisec): " + ManagementFactory.getRuntimeMXBean().getUptime());
		
	}

}
