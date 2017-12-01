package genomicUtils;

public class DNAUtils {

	public String comp(String strand) {
		strand = strand.toUpperCase();
		char[] Strand = strand.toCharArray();
		String comp = "";
		for(int i = 0; i < Strand.length; i++) {
			switch(Strand[i]) {
				case 'A': comp += "T";
				break;
				case 'T': comp += "A";
				break;
				case 'G': comp += "C";
				break;
				case 'C': comp += "G";
				break;
				case 'N': comp += "N";
				break;
			}
		}
		return comp;
	}
	
	public String revcomp(String strand) {
		String comp = comp(strand);
		char[] revcompa = comp.toCharArray();
		String revcomp = "";

		for(int i = revcompa.length-1; i >= 0; --i)
		{
			revcomp += revcompa[i];
		}	
		
		return revcomp;
	}
	
	
}
